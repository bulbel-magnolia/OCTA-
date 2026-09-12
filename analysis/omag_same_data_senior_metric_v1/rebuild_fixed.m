function rebuild_fixed()
% Fixed raw-reader transcription of senior exporter, plus independent frozen reader.
out=fileparts(mfilename('fullpath'));
repo=fileparts(fileparts(out));
addpath(fullfile(repo,'matlab'));
addpath(fullfile(out,'reconstruction_source'),'-begin');
jobs=jsondecode(fileread(fullfile(out,'intermediate_local','jobs.json')));
env=struct('matlab_version',version,'release',version('-release'),'computer',computer,...
 'iniset',which('iniset'),'subpix',which('OCTA_F_SubPixReg'),...
 'phase',which('SSOCT_F_PhCompV3'),'eigfeed',which('OCTA_F_ED_Clutter_EigFeed'));
fid=fopen(fullfile(out,'intermediate_local','matlab_environment.json'),'w','n','UTF-8');
fwrite(fid,jsonencode(env),'char');fclose(fid);
ini=iniset(fullfile(out,'reconstruction_source','config.ini'));
Thr=10^(str2double(ini.input.thr)/20);
assert(str2double(ini.input.thr)==85);
for j=1:numel(jobs)
 job=jobs(j);
 IMG=senior_reader(job.raw_path,job.frame_index+1,Thr);
 sv_rebuilt=var(abs(IMG),1,3);
 [~,omag_rebuilt]=OCTA_F_ED_Clutter_EigFeed(IMG,2);
 [IMG_independent,info]=load_and_reconstruct_common_oct(job.raw_path,job.frame_index+1);
 [~,omag_independent]=OCTA_F_ED_Clutter_EigFeed(IMG_independent,2);
 complex_max_error=max(abs(IMG-IMG_independent),[],'all');
 assert(isequal(size(IMG),[351 500 3]) && all(isfinite(IMG),'all'));
 save(job.output_path,'sv_rebuilt','omag_rebuilt','omag_independent','complex_max_error','info','-v7');
 fprintf('FIXED %d/%d %s frame=%d complex_error=%g\n',j,numel(jobs),job.scan_id,job.frame_index,complex_max_error);
end
end

function IMG=senior_reader(rawPath,iY,Thr)
fid=fopen(rawPath,'r','ieee-le');assert(fid>=0);cleanup=onCleanup(@() fclose(fid));
bob=fread(fid,1,'uint32');SPL=fread(fid,1,'double');nX=fread(fid,1,'uint32');
nYraw=fread(fid,1,'uint32');Boffset=fread(fid,1,'uint32');Blength=fread(fid,1,'uint32')+1;
fread(fid,4,'double');nR=fread(fid,1,'uint32');
assert(SPL==1024 && nX==500 && nYraw==1500 && nR==3 && Blength==1024);
assert(fseek(fid,bob+(SPL*nX+2)*nR*(iY-1)*2,'bof')==0);
Bs=zeros(Blength,nX,nR);
for ic=1:nR
 assert(fseek(fid,4,'cof')==0);
 vals=fread(fid,SPL*nX,'int16=>double');assert(numel(vals)==SPL*nX);
 B=reshape(vals,[SPL,nX]);Bs(:,:,ic)=B(Boffset+1:Boffset+Blength,:);
end
Bimg=fft(Bs,SPL,1);IMG=Bimg(50:400,:,:);
[IMG,~]=OCTA_F_SubPixReg(IMG,10,false);
IMG=SSOCT_F_PhCompV3(IMG,1,Thr,true);
end
