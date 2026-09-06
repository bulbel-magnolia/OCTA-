function [sv_raw, stru_amp] = compute_sv_maps(IMG)
%COMPUTE_SV_MAPS Formal linear speckle-variance and structure maps.
%   IMG must be depth x A-line x repeat complex data. The variance flag 1
%   selects denominator N, exactly matching var(abs(E), 1, 3).

if ndims(IMG) ~= 3 || size(IMG, 3) < 1
    error('SVRectTail:InvalidComplexVolume', ...
        'IMG must be a nonempty depth-by-A-line-by-repeat array.');
end
amplitude = abs(IMG);
sv_raw = var(amplitude, 1, 3);
stru_amp = mean(amplitude, 3);
end
