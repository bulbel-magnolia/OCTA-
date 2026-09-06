function export_sv_omag_frame(source_file, bscan_index, output_file, omag_num_ev)
%EXPORT_SV_OMAG_FRAME Reconstruct and export one auditable interim frame.
%   BSCAN_INDEX is MATLAB 1-based. OUTPUT_FILE is a classic MAT file read
%   directly by svrecttail. The formal signal is never normalized, clipped,
%   log-transformed, or display-mapped.

if nargin < 4 || isempty(omag_num_ev)
    omag_num_ev = 2;
end
if nargin < 3 || isempty(output_file)
    error('SVRectTail:MissingOutput', 'An output MAT path is required.');
end
if exist(output_file, 'file')
    error('SVRectTail:OutputExists', ...
        'Refusing to overwrite existing output: %s', output_file);
end
if ~(isscalar(omag_num_ev) && isfinite(omag_num_ev) && ...
        omag_num_ev >= 0 && omag_num_ev == round(omag_num_ev))
    error('SVRectTail:InvalidOMAGNumEV', ...
        'omag_num_ev must be a nonnegative integer.');
end

output_parent = fileparts(output_file);
if ~isempty(output_parent) && ~exist(output_parent, 'dir')
    mkdir(output_parent);
end
[IMG, reconstruction] = load_and_reconstruct_common_oct(source_file, bscan_index);
if omag_num_ev >= size(IMG, 3)
    error('SVRectTail:OMAGNumEVOutOfRange', ...
        'omag_num_ev must be smaller than the repeat count (%d).', size(IMG, 3));
end
[sv_raw, stru_amp] = compute_sv_maps(IMG);
cv2_epsilon = eps('double');
sv_cv2 = sv_raw ./ (stru_amp .^ 2 + cv2_epsilon);
[~, omag_raw] = OCTA_F_ED_Clutter_EigFeed(IMG, omag_num_ev);
if ~isequal(size(sv_raw), size(omag_raw), size(stru_amp))
    error('SVRectTail:MapSizeMismatch', ...
        'Export maps must share depth-by-A-line dimensions.');
end

metadata = struct();
metadata.schema_version = 'SV_Rectangle_v1_interim';
metadata.source_file = char(source_file);
metadata.bscan_index_matlab_1based = double(bscan_index);
metadata.formal_signal_name = 'sv_raw';
metadata.formal_signal_definition = 'var(abs(E), 1, 3)';
metadata.variance_denominator = 'N';
metadata.sv_cv2_definition = 'sv_raw / (stru_amp^2 + cv2_epsilon)';
metadata.cv2_epsilon = cv2_epsilon;
metadata.dimension_order = 'depth x A-line';
metadata.anchor_coordinate_convention = ...
    'downstream manifest uses zero-based pixel-centre coordinates';
metadata.omag_num_ev = double(omag_num_ev);
metadata.reconstruction = reconstruction;
metadata.exported_at_utc = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX'));

save(output_file, 'sv_raw', 'sv_cv2', 'omag_raw', 'stru_amp', ...
    'metadata', '-v7');
end
