function test_reconstruction_dependencies
%TEST_RECONSTRUCTION_DEPENDENCIES Exercise the copied registration/OMAG closure.

rng(7, 'twister');
test_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(test_dir));
ini = iniset(fullfile(repo_root, 'matlab', 'config.ini'));
assert(str2double(ini.input.thr) == 85);
base = randn(32, 20) + 1i * randn(32, 20);
IMG = cat(3, base, base .* exp(1i * 0.15), circshift(base, [1, 0]));

[registered, shift_px] = OCTA_F_SubPixReg(IMG, 10, false);
assert(isequal(size(registered), size(IMG)));
assert(isequal(size(shift_px), [1, size(IMG, 3)]));
assert(all(isfinite(real(registered(:)))));
assert(all(isfinite(imag(registered(:)))));

compensated = SSOCT_F_PhCompV3(registered, 1, 0, true);
assert(isequal(size(compensated), size(IMG)));
assert(all(isfinite(real(compensated(:)))));
assert(all(isfinite(imag(compensated(:)))));

[omag_tissue, omag_flow] = OCTA_F_ED_Clutter_EigFeed(compensated, 1);
assert(isequal(size(omag_tissue), size(base)));
assert(isequal(size(omag_flow), size(base)));
assert(all(isfinite(omag_tissue(:))));
assert(all(isfinite(omag_flow(:))));

fprintf('test_reconstruction_dependencies PASS\n');
end
