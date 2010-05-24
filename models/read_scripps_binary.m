function [info]=read_scripps_binary(file)
%READ_SCRIPPS_BINARY    Reads a Scripps mantle model unformatted binary
%
%    Usage:    info=read_scripps_binary(binfile)

% open the file
fid=fopen(file,'r','ieee-be');

% first read nblks/phibar
fread(fid,1,'*int32'); % unformatted header
info.nblks=fread(fid,1,'*int32');
info.phibar1=fread(fid,1,'*float32');
fread(fid,1,'*int32'); % unformatted tailer

% now read in model
fread(fid,1,'*int32'); % unformatted header
info.iter=fread(fid,1,'*int32');
info.model=fread(fid,double(info.nblks),'*float32');
info.phibar2=fread(fid,1,'*float32');
info.r=fread(fid,1,'*float32');
fread(fid,1,'*int32'); % unformatted tailer
fclose(fid);

end
