function export_sv_collection(jobs_file, selected_scan)
% Batch wrapper only; the frozen per-frame reconstruction/export is unchanged.
if nargin < 2, selected_scan = ''; end
addpath(fileparts(mfilename('fullpath')));
jobs = jsondecode(fileread(jobs_file));
for j = 1:numel(jobs)
    job = jobs(j);
    if ~isempty(selected_scan) && ~strcmp(job.scan_id, selected_scan), continue; end
    if ~exist(job.interim_dir, 'dir'), mkdir(job.interim_dir); end
    for frame = reshape(job.frame_indices_0based, 1, [])
        output = fullfile(job.interim_dir, sprintf('frame_%03d.mat', frame));
        if exist(output, 'file')
            prior = load(output, 'metadata');
            assert(strcmp(prior.metadata.source_file, job.source_file) && ...
                prior.metadata.bscan_index_matlab_1based == frame + 1, ...
                'SVCollection:ResumeMismatch', 'Existing export belongs to another input.');
            continue;
        end
        temporary = [tempname(job.interim_dir) '.mat'];
        export_sv_omag_frame(job.source_file, frame + 1, temporary, 2);
        movefile(temporary, output);
        if mod(frame, 25) == 0 || frame == job.frame_indices_0based(end)
            fprintf('EXPORTED %s frame %d\n', job.scan_id, frame);
        end
    end
end
end
