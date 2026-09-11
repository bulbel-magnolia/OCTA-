function meta = export_p_bld_ed_volume(rawPath, outPath, opts)
%EXPORT_P_BLD_ED_VOLUME Reproduce VCSEL_OMAGV4 up to unfiltered p_bld_ed.
%
% This is deliberately a read-only raw-data reader.  It does not call the
% interactive VCSEL_OMAGV4 entry point, does not create/replace DICOM files,
% and never writes beside the source .oct file.  Output is a float32
% MATLAB-v7.3 file written through a same-directory temporary file and only
% renamed to outPath after all requested frames and metadata are complete.
%
% Full-volume example:
%   export_p_bld_ed_volume(raw, out, DatasetId="d185_s1_scan32")
%
% Small smoke-test example (the third output axis follows this order):
%   export_p_bld_ed_volume(raw, out, FrameIndices=[1 250 500])

arguments
    rawPath (1,1) string
    outPath (1,1) string
    opts.DatasetId (1,1) string = ""
    opts.CodeDir (1,1) string = "D:\OCT_Project\matlabwork\matlab_20260107"
    opts.ConfigIni (1,1) string = ""
    opts.FrameIndices (1,:) double {mustBeInteger,mustBePositive} = []
    opts.SavePTis (1,1) logical = false
    opts.SaveRegistrationShifts (1,1) logical = true
    opts.Overwrite (1,1) logical = false
    opts.ExpectedFrames (1,1) double {mustBeInteger,mustBePositive} = 500
    opts.ExpectedRawBytes (1,1) double {mustBeInteger,mustBePositive} = 1536021062
    opts.LogPath (1,1) string = ""
    opts.ExtraMetadata (1,1) struct = struct()
end

exporterVersion = "1.0.0";
rawPath = absolute_existing_file(rawPath, "rawPath");
codeDir = absolute_existing_folder(opts.CodeDir, "CodeDir");
outPath = string(char(java.io.File(char(outPath)).getCanonicalPath()));

if strlength(opts.ConfigIni) == 0
    configIni = fullfile(codeDir, "config.ini");
else
    configIni = opts.ConfigIni;
end
configIni = absolute_existing_file(configIni, "ConfigIni");

if isfile(outPath) && ~opts.Overwrite
    error("tailq:OutputExists", ...
        "Output already exists; refusing to overwrite: %s", outPath);
end
if startsWith(lower(outPath), lower(fileparts(rawPath) + filesep))
    error("tailq:UnsafeOutput", ...
        "Output must not be placed in the raw-data directory: %s", outPath);
end

outDir = string(fileparts(outPath));
if ~isfolder(outDir), mkdir(outDir); end
if strlength(opts.LogPath) == 0
    logPath = outPath + ".export.log";
else
    logPath = string(char(java.io.File(char(opts.LogPath)).getCanonicalPath()));
end
logDir = string(fileparts(logPath));
if ~isfolder(logDir), mkdir(logDir); end

partialPath = outPath + ".partial-" + string(feature("getpid")) + ".mat";
partialMetaPath = outPath + ".metadata.json.partial-" + string(feature("getpid"));
finalMetaPath = outPath + ".metadata.json";
if isfile(partialPath), delete(partialPath); end
if isfile(partialMetaPath), delete(partialMetaPath); end

logFid = fopen(logPath, "a", "n", "UTF-8");
assert(logFid >= 0, "Cannot open export log: %s", logPath);
logCleanup = onCleanup(@() fclose(logFid)); %#ok<NASGU>
log_line(logFid, "START exporter=%s dataset=%s raw=%s output=%s", ...
    exporterVersion, opts.DatasetId, rawPath, outPath);

try
    addpath(codeDir);
    dependencies = ["iniset.m", "OCTA_F_SubPixReg.m", ...
        "SSOCT_F_PhCompV3.m", "OCTA_F_ED_Clutter_EigFeed.m"];
    for i = 1:numel(dependencies)
        assert(isfile(fullfile(codeDir, dependencies(i))), ...
            "Missing MATLAB dependency: %s", fullfile(codeDir, dependencies(i)));
    end

    rawInfoBefore = dir(rawPath);
    assert(~isempty(rawInfoBefore), "Cannot stat raw source: %s", rawPath);
    if opts.ExpectedRawBytes > 0
        assert(rawInfoBefore.bytes == opts.ExpectedRawBytes, ...
            "Unexpected raw size for %s: got %d, expected %d bytes", ...
            rawPath, rawInfoBefore.bytes, opts.ExpectedRawBytes);
    end

    ini = iniset(char(configIni));
    thrDb = str2double(ini.input.thr);
    assert(isfinite(thrDb), "Invalid Thr in %s", configIni);
    Thr = 10^(thrDb/20);
    IMcropRg = 50:400;
    Nsub = 10;
    Colshift = false;
    NumEV = 2;
    doMedianShift = true;

    fid = fopen(rawPath, "r");
    assert(fid >= 0, "Cannot open raw file: %s", rawPath);
    rawCleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    header = read_raw_header(fid);
    nY = floor(double(header.nYraw) / double(header.frame_per_pos));
    assert(header.SPL == 1024 && header.nX == 500 && header.Blength == 1024, ...
        "Unexpected raw dimensions SPL=%g nX=%d Blength=%d", ...
        header.SPL, header.nX, header.Blength);
    assert(header.frame_per_pos == 3, ...
        "Unexpected frame_per_pos=%d (expected 3)", header.frame_per_pos);
    assert(nY == opts.ExpectedFrames, ...
        "Unexpected spatial frame count=%d (expected %d)", nY, opts.ExpectedFrames);

    if isempty(opts.FrameIndices)
        frameIdx = 1:nY;
    else
        frameIdx = opts.FrameIndices;
    end
    assert(all(frameIdx <= nY), "Requested frame exceeds nY=%d", nY);
    assert(numel(unique(frameIdx)) == numel(frameIdx), ...
        "FrameIndices contains duplicates; output axis would be ambiguous");

    nF = numel(frameIdx);
    nZ = numel(IMcropRg);
    nX = double(header.nX);
    nR = double(header.frame_per_pos);

    m = matfile(partialPath, "Writable", true);
    m.p_bld_ed(nZ, nX, nF) = single(0);
    if opts.SavePTis
        m.p_tis(nZ, nX, nF) = single(0);
    end
    if opts.SaveRegistrationShifts
        m.registration_Dx(nR, nF) = single(0);
    end

    t0 = tic;
    frameSeconds = nan(1, nF);
    for jj = 1:nF
        frameTic = tic;
        iY = frameIdx(jj);
        Bs = zeros(double(header.Blength), nX, nR);
        offset = double(header.bob) + ...
            (double(header.SPL) * nX + 2) * nR * (double(iY)-1) * 2;
        status = fseek(fid, offset, "bof");
        assert(status == 0, "fseek failed for frame %d", iY);
        for ic = 1:nR
            status = fseek(fid, 4, "cof");
            assert(status == 0, "repeat-header fseek failed at frame %d repeat %d", iY, ic);
            vals = fread(fid, double(header.SPL)*nX, "int16=>double");
            assert(numel(vals) == double(header.SPL)*nX, ...
                "Short raw read at frame %d repeat %d", iY, ic);
            B = reshape(vals, [double(header.SPL), nX]);
            first = double(header.Boffset) + 1;
            last = double(header.Boffset) + double(header.Blength);
            Bs(:,:,ic) = B(first:last, :);
        end

        Bimg = fft(Bs, double(header.SPL), 1);
        IMG = Bimg(IMcropRg, :, :);
        [IMG, Dx] = OCTA_F_SubPixReg(IMG, Nsub, Colshift);
        IMG = SSOCT_F_PhCompV3(IMG, 1, Thr, doMedianShift);
        [~, pBldEd] = OCTA_F_ED_Clutter_EigFeed(IMG, NumEV);
        assert(isequal(size(pBldEd), [nZ nX]) && all(isfinite(pBldEd), "all"), ...
            "Invalid p_bld_ed at frame %d", iY);
        m.p_bld_ed(:,:,jj) = single(pBldEd);

        if opts.SavePTis
            pTis = mean(abs(IMG), 3);
            assert(all(isfinite(pTis), "all"), "Invalid p_tis at frame %d", iY);
            m.p_tis(:,:,jj) = single(pTis);
        end
        if opts.SaveRegistrationShifts
            m.registration_Dx(:,jj) = single(Dx(:));
        end

        frameSeconds(jj) = toc(frameTic);
        if jj == 1 || mod(jj, 10) == 0 || jj == nF
            elapsed = toc(t0);
            eta = median(frameSeconds(1:jj), "omitnan") * (nF-jj);
            log_line(logFid, "PROGRESS %d/%d source_frame=%d elapsed_sec=%.1f eta_sec=%.1f", ...
                jj, nF, iY, elapsed, eta);
        end
    end
    elapsedSec = toc(t0);

    rawInfoAfter = dir(rawPath);
    assert(rawInfoAfter.bytes == rawInfoBefore.bytes && ...
           rawInfoAfter.datenum == rawInfoBefore.datenum, ...
        "Raw source changed while it was being read; export rejected: %s", rawPath);

    meta = build_metadata(exporterVersion, opts, rawPath, outPath, ...
        rawInfoBefore, codeDir, configIni, dependencies, header, nY, ...
        frameIdx, IMcropRg, thrDb, Thr, Nsub, Colshift, NumEV, ...
        doMedianShift, elapsedSec, frameSeconds);
    metaJson = jsonencode(meta, PrettyPrint=true);

    % These variables make the logical axis order explicit both to MATLAB
    % readers and to h5py readers of the underlying v7.3 HDF5 container.
    m.meta_json = metaJson;
    m.axis_order_matlab = "z,x,frame";
    m.axis_order_h5py = "frame,x,z";
    m.logical_shape_z_x_frame = uint32([nZ nX nF]);
    m.source_frame_indices = uint32(frameIdx);
    clear m

    sidecarFid = fopen(partialMetaPath, "w", "n", "UTF-8");
    assert(sidecarFid >= 0, "Cannot create metadata sidecar: %s", partialMetaPath);
    fwrite(sidecarFid, metaJson, "char");
    fclose(sidecarFid);

    if isfile(outPath) && opts.Overwrite, delete(outPath); end
    if isfile(finalMetaPath) && opts.Overwrite, delete(finalMetaPath); end
    [ok, msg] = movefile(partialPath, outPath, "f");
    assert(ok, "Failed to finalize output: %s", msg);
    [ok, msg] = movefile(partialMetaPath, finalMetaPath, "f");
    assert(ok, "Failed to finalize metadata sidecar: %s", msg);
    log_line(logFid, "DONE frames=%d elapsed_sec=%.3f bytes=%d", ...
        nF, elapsedSec, dir(outPath).bytes);
catch ME
    log_line(logFid, "ERROR id=%s message=%s", string(ME.identifier), string(ME.message));
    % A .partial-* file is intentionally never promoted.  It is left for
    % forensic diagnosis and can never be mistaken for a valid .mat output.
    rethrow(ME);
end
end

function header = read_raw_header(fid)
status = fseek(fid, 0, "bof");
assert(status == 0, "Cannot seek to raw header");
header = struct();
header.bob = fread(fid, 1, "uint32");
header.SPL = fread(fid, 1, "double");
header.nX = fread(fid, 1, "uint32");
header.nYraw = fread(fid, 1, "uint32");
header.Boffset = fread(fid, 1, "uint32");
header.Blength_stored = fread(fid, 1, "uint32");
header.Blength = header.Blength_stored + 1;
header.Xcenter = fread(fid, 1, "double");
header.Xspan = fread(fid, 1, "double");
header.Ycenter = fread(fid, 1, "double");
header.Yspan = fread(fid, 1, "double");
header.frame_per_pos = fread(fid, 1, "uint32");
fields = fieldnames(header);
for i = 1:numel(fields)
    assert(~isempty(header.(fields{i})), "Truncated raw header at %s", fields{i});
end
end

function meta = build_metadata(versionText, opts, rawPath, outPath, rawInfo, ...
        codeDir, configIni, dependencies, header, nY, frameIdx, IMcropRg, ...
        thrDb, Thr, Nsub, Colshift, NumEV, doMedianShift, elapsedSec, frameSeconds)
meta = struct();
meta.schema_version = "1.0.0";
meta.exporter_name = "export_p_bld_ed_volume";
meta.exporter_version = versionText;
meta.dataset_id = opts.DatasetId;
meta.created_utc = string(datetime("now", TimeZone="UTC", ...
    Format="yyyy-MM-dd'T'HH:mm:ss.SSSXXX"));
meta.raw_path = rawPath;
meta.raw_bytes = rawInfo.bytes;
meta.raw_last_write_local = string(datetime(rawInfo.datenum, ...
    ConvertFrom="datenum", Format="yyyy-MM-dd'T'HH:mm:ss.SSS"));
meta.output_path = outPath;
meta.output_format = "MAT-v7.3";
meta.output_variable = "p_bld_ed";
meta.output_dtype = "single";
meta.axis_order_matlab = "z,x,frame";
meta.axis_order_h5py = "frame,x,z";
meta.logical_shape_z_x_frame = [numel(IMcropRg), double(header.nX), numel(frameIdx)];
meta.source_frame_indices = frameIdx;
meta.is_full_volume = numel(frameIdx) == nY && isequal(frameIdx, 1:nY);
meta.save_p_tis = opts.SavePTis;
meta.save_registration_shifts = opts.SaveRegistrationShifts;
meta.raw_header = header;
meta.raw_header.nY_spatial = nY;
meta.processing = struct( ...
    "IMcropRg_matlab", IMcropRg, ...
    "threshold_db", thrDb, ...
    "threshold_linear", Thr, ...
    "subpixel_registration", true, ...
    "Nsub", Nsub, ...
    "Colshift", Colshift, ...
    "phase_compensation", true, ...
    "median_phase_shift", doMedianShift, ...
    "eigen_clutter_NumEV", NumEV, ...
    "filter_omag_bscan_applied", false, ...
    "dicom_normalization_applied", false, ...
    "threshold_clipping_applied", false);
meta.matlab_version = string(version);
meta.matlab_release = string(version("-release"));
meta.computer = string(computer);
meta.code_dir = codeDir;
meta.config_ini = configIni;
meta.dependencies = file_inventory(codeDir, dependencies);
meta.elapsed_sec = elapsedSec;
meta.frame_seconds_median = median(frameSeconds, "omitnan");
meta.frame_seconds_iqr = iqr(frameSeconds);
meta.extra = opts.ExtraMetadata;
end

function inventory = file_inventory(codeDir, names)
inventory = repmat(struct("name", "", "path", "", "bytes", 0, ...
    "last_write_local", ""), 1, numel(names)+2);
allNames = ["VCSEL_OMAGV4.m", "config.ini", names];
for i = 1:numel(allNames)
    p = fullfile(codeDir, allNames(i));
    d = dir(p);
    inventory(i).name = allNames(i);
    inventory(i).path = string(p);
    inventory(i).bytes = d.bytes;
    inventory(i).last_write_local = string(datetime(d.datenum, ...
        ConvertFrom="datenum", Format="yyyy-MM-dd'T'HH:mm:ss.SSS"));
end
end

function value = absolute_existing_file(value, label)
assert(isfile(value), "%s does not exist: %s", label, value);
value = string(char(java.io.File(char(value)).getCanonicalPath()));
end

function value = absolute_existing_folder(value, label)
assert(isfolder(value), "%s does not exist: %s", label, value);
value = string(char(java.io.File(char(value)).getCanonicalPath()));
end

function log_line(fid, formatText, varargin)
stamp = string(datetime("now", Format="yyyy-MM-dd HH:mm:ss.SSS"));
message = sprintf(formatText, varargin{:});
line = sprintf("[%s] %s\n", stamp, message);
fprintf(1, "%s", line);
fprintf(fid, "%s", line);
end
