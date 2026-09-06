function [IMG, info] = load_and_reconstruct_common_oct(source_file, bscan_index)
%LOAD_AND_RECONSTRUCT_COMMON_OCT Rebuild one common repeated-frame OCT B-scan.
%   [IMG, INFO] = LOAD_AND_RECONSTRUCT_COMMON_OCT(SOURCE_FILE, BSCAN_INDEX)
%   follows the shared upstream path in VCSEL_OMAGV4 and
%   VCSEL_SpeckleVarianceV4:
%
%       .oct frame read -> FFT -> crop 50:400 -> subpixel registration
%       -> phase compensation
%
%   IMG is a complex depth-by-A-line-by-repeat array. No OMAG or speckle
%   variance calculation, masking, filtering, log conversion, display
%   mapping, or DICOM conversion is performed here.

if nargin < 1 || isempty(source_file)
    error('OCTA:MissingSourceFile', 'A source .oct file is required.');
end
if nargin < 2
    bscan_index = [];
end
source_file = char(source_file);

% These settings are shared by the current OMAG V4 and SV V4 programs.
IMcropRg = 50:400;
boffset = 0;
do_reg = true;
Nsub = 10;
Colshift = false;
do_PhComp = true;
do_medianshift = true;
subtract_reference = false;
useDVV = false;
useKES = false;
usedispcoef = false;

loader_dir = fileparts(mfilename('fullpath'));
phase_config_file = fullfile(loader_dir, 'config.ini');
if ~exist(phase_config_file, 'file')
    error('OCTA:ThresholdConfigUnavailable', ...
        'The frozen phase configuration is missing: %s', phase_config_file);
end
ini = iniset(phase_config_file);
threshold_db = str2double(ini.input.thr);
if ~(isscalar(threshold_db) && isfinite(threshold_db))
    error('OCTA:InvalidPhaseThreshold', ...
        'input.thr in matlab/config.ini must be a finite scalar.');
end
Thr = 10^(threshold_db / 20);

fid = fopen(source_file, 'r', 'ieee-le');
if fid == -1
    error('OCTA:CannotOpenSourceFile', 'Cannot open source OCT file: %s', source_file);
end
file_cleanup = onCleanup(@() fclose(fid));

header = struct();
header.bob = local_read_scalar(fid, 'uint32=>double', 'bob');
header.SPL = local_read_scalar(fid, 'double=>double', 'SPL');
header.nX = local_read_scalar(fid, 'uint32=>double', 'nX');
header.nY_raw = local_read_scalar(fid, 'uint32=>double', 'nY');
header.Boffset = local_read_scalar(fid, 'uint32=>double', 'Boffset');
header.Blength = local_read_scalar(fid, 'uint32=>double', 'Blength') + 1;
header.Xcenter = local_read_scalar(fid, 'double=>double', 'Xcenter');
header.Xspan = local_read_scalar(fid, 'double=>double', 'Xspan');
header.Ycenter = local_read_scalar(fid, 'double=>double', 'Ycenter');
header.Yspan = local_read_scalar(fid, 'double=>double', 'Yspan');
header.frame_per_pos = local_read_scalar( ...
    fid, 'uint32=>double', 'frame_per_pos');
header.reserved = local_read_scalar(fid, 'uint32=>double', 'reserved');
header.sizeBck = local_read_scalar(fid, 'uint32=>double', 'sizeBck');
[~, background_count] = fread(fid, header.sizeBck, 'int16=>double');
if background_count ~= header.sizeBck
    error('OCTA:TruncatedHeader', ...
        'The OCT background spectrum is truncated (%d of %d values read).', ...
        background_count, header.sizeBck);
end

local_validate_positive_integer(header.SPL, 'SPL');
local_validate_positive_integer(header.nX, 'nX');
local_validate_positive_integer(header.nY_raw, 'nY');
local_validate_positive_integer(header.Blength, 'Blength');
local_validate_positive_integer(header.frame_per_pos, 'frame_per_pos');
if header.Boffset < 0 || header.Boffset ~= round(header.Boffset)
    error('OCTA:InvalidHeader', 'Boffset must be a nonnegative integer.');
end
if header.Boffset + header.Blength > header.SPL
    error('OCTA:InvalidHeader', ...
        'Boffset + Blength exceeds the recorded spectrum length SPL.');
end
if max(IMcropRg) > header.SPL
    error('OCTA:CropRangeOutOfBounds', ...
        'The verified crop 50:400 exceeds SPL=%d.', header.SPL);
end

nR = header.frame_per_pos;
nY_adjusted = header.nY_raw;
if boffset > 0
    nY_adjusted = nY_adjusted - 1;
end
n_bscans = floor(nY_adjusted / nR);
if n_bscans < 1
    error('OCTA:NoCompleteRepeatedBscan', ...
        'The file contains no complete repeated-frame B-scan group.');
end
if mod(nY_adjusted, nR) ~= 0
    warning('OCTA:TrailingFramesIgnored', ...
        ['nY=%d is not divisible by frame_per_pos=%d; the final %d frame(s) ' ...
         'are ignored, matching the current floor(nY/nR) organization.'], ...
        nY_adjusted, nR, mod(nY_adjusted, nR));
end

if isempty(bscan_index)
    if n_bscans == 1
        bscan_index = 1;
    else
        fprintf('The OCT file contains %d repeated-frame B-scan positions.\n', ...
            n_bscans);
        bscan_index = input(sprintf('B-scan index [1-%d] = ', n_bscans));
    end
end
local_validate_positive_integer(bscan_index, 'B-scan index');
if bscan_index > n_bscans
    error('OCTA:BscanIndexOutOfRange', ...
        'B-scan index must satisfy 1 <= index <= %d; received %g.', ...
        n_bscans, bscan_index);
end
bscan_index = double(bscan_index);

bytes_per_frame = (header.SPL * header.nX + 2) * 2;
bscan_offset = header.bob + bytes_per_frame * nR * (bscan_index - 1) + ...
    bytes_per_frame * boffset;
if fseek(fid, bscan_offset, 'bof') ~= 0
    error('OCTA:DataSeekFailed', ...
        'Cannot seek to B-scan %d at byte offset %.0f.', ...
        bscan_index, bscan_offset);
end

Bs = zeros(header.Blength, header.nX, nR, 'double');
samples_per_frame = header.SPL * header.nX;
for i_repeat = 1:nR
    if fseek(fid, 4, 'cof') ~= 0
        error('OCTA:FrameHeaderSeekFailed', ...
            'Cannot skip the frame header for repeat %d.', i_repeat);
    end
    [raw_frame, sample_count] = fread( ...
        fid, samples_per_frame, 'int16=>double');
    if sample_count ~= samples_per_frame
        error('OCTA:TruncatedFrame', ...
            ['B-scan %d repeat %d is truncated (%d of %.0f spectral ' ...
             'samples read).'], ...
            bscan_index, i_repeat, sample_count, samples_per_frame);
    end
    raw_frame = reshape(raw_frame, [header.SPL, header.nX]);
    Bs(:, :, i_repeat) = raw_frame( ...
        header.Boffset + 1:header.Boffset + header.Blength, :);
end

if subtract_reference || useDVV || useKES || usedispcoef
    error('OCTA:UnsupportedCommonSetting', ...
        'The current verified common path requires all optional branches to remain disabled.');
end

Bimg = fft(Bs, header.SPL, 1);
IMG = Bimg(IMcropRg, :, :);

[IMG, registration_shift_px] = OCTA_F_SubPixReg(IMG, Nsub, Colshift);
IMG = SSOCT_F_PhCompV3(IMG, 1, Thr, do_medianshift);

expected_size = [numel(IMcropRg), header.nX, nR];
if ~isequal(size(IMG), expected_size)
    error('OCTA:UnexpectedIMGSize', ...
        'Common reconstruction produced IMG size [%s], expected [%s].', ...
        num2str(size(IMG)), num2str(expected_size));
end

info = struct();
info.source_file = source_file;
info.bscan_index = bscan_index;
info.n_bscan_positions = n_bscans;
info.IMG_size = size(IMG);
info.dimension_order = 'depth x A-line x repeated frames';
info.IMcropRg = IMcropRg;
info.phase_threshold_db = threshold_db;
info.phase_threshold_linear = Thr;
info.phase_threshold_config = 'matlab/config.ini';
info.do_registration = do_reg;
info.registration_Nsub = Nsub;
info.registration_Colshift = Colshift;
info.registration_shift_px = registration_shift_px;
info.do_phase_compensation = do_PhComp;
info.do_median_bulk_phase_shift = do_medianshift;
info.subtract_reference = subtract_reference;
info.useDVV = useDVV;
info.useKES = useKES;
info.use_dispersion_coefficient = usedispcoef;
info.header = header;
info.axial_calibration_available = false;
end

function value = local_read_scalar(fid, precision, field_name)
[value, count] = fread(fid, 1, precision);
if count ~= 1
    error('OCTA:TruncatedHeader', ...
        'Cannot read required OCT header field: %s.', field_name);
end
end

function local_validate_positive_integer(value, label)
if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
        isfinite(value) && value >= 1 && value == round(value))
    error('OCTA:InvalidPositiveInteger', ...
        '%s must be a finite positive integer; received %g.', label, value);
end
end
