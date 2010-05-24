function [s]=scripps2matlab(dfile,bfile)
%SCRIPPS2MATLAB    Creates a Matlab struct from a Scripps mantle model
%
%    Usage:    s=scripps2matlab(descfile,binfile)

% read binary info
s=read_scripps_binary(bfile);

% read in model description file
words=getwords(readtxt(dfile));
words=str2double(words);
s.blocksize=words(1);
s.nlayers=words(2);
s.deplimits=6371-reshape(words(3:end),2,[]).';

% now get lat/lon block limits (lon is not regular)
s.latlimits=[90:-4:-86; 86:-4:-90].';
s.lonblks_at_lat=max(round(360./s.blocksize...
    *sin((90-mean(s.latlimits,2))*pi/180)),1);
s.lonwidth=360./s.lonblks_at_lat;
s.lat_block_idx=cumsum(s.lonblks_at_lat);
s.lat_block_idx=[[1; s.lat_block_idx(1:end-1)+1] s.lat_block_idx];
s.blocks_per_layer=s.lat_block_idx(end);
s.lonlimits=nan(s.blocks_per_layer,2);
for i=1:numel(s.lonblks_at_lat)
    s.lonlimits(s.lat_block_idx(i,1):s.lat_block_idx(i,2),:)=...
        [0:s.lonwidth(i):360-s.lonwidth(i);
        s.lonwidth(i):s.lonwidth(i):360].';
end

end
