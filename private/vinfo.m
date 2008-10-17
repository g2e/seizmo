function [h,idx]=vinfo(data)
%VINFO    Returns version info for SAClab data records
%
%    Description: [H,IDX]=VINFO(DATA) returns all necessary version
%     definitions pertaining to the records in the SAClab structure DATA in
%     the struct array H.  IDX has one entry per record in DATA that
%     gives the index in H that corresponds to that record's version
%     definition. This is an internal function to reduce rampant code
%     repetition.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    [h,idx]=vinfo(data)
%
%    Examples:
%
%    See also: seisdef

%     Version History:
%        Sep. 25, 2008 - initial version
%        Sep. 26, 2008 - check internal header version
%        Oct. 16, 2008 - removed header version consistency check
%        Oct. 17, 2008 - filetype support, removed several pointless
%                        outputs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 17, 2008 at 01:45 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data))

% get filetypes
ft={data.filetype}.';
uft=unique(ft);
nft=numel(uft);

% get versions
v=[data.version].';

% loop through each filetype
count=0; idx=nan(size(data));
for i=1:nft
    % who has this filetype
    ift=strcmpi(uft{i},ft);
    
    % get local versions
    uv=unique(v(ift));
    nv=numel(uv);
    
    % grab definitions
    for j=1:nv
        count=count+1;
        h(count)=seisdef(uft{i},uv(j));
        idx(v==uv(j) & ift)=count;
    end
end

end
