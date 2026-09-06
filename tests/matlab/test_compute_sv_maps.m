function test_compute_sv_maps()
%TEST_COMPUTE_SV_MAPS Verify the frozen population-variance definition.

test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(test_dir));
matlab_dir = fullfile(repo_root, 'matlab');
path_entries = strsplit(path, pathsep);
if ~any(strcmp(path_entries, matlab_dir))
    addpath(matlab_dir);
end

IMG = reshape(complex([1, 2, 3, 4], 0), [1, 1, 4]);
[sv_raw, stru_amp] = compute_sv_maps(IMG);
assert(abs(sv_raw - 1.25) < 1e-12, ...
    'sv_raw must use denominator N, not N-1.');
assert(abs(stru_amp - 2.5) < 1e-12, ...
    'stru_amp must be the repeat-wise amplitude mean.');
fprintf('test_compute_sv_maps PASS\n');
end
