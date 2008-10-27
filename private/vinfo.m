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
%     - currently assumes each definition struct has the same top tier
%       layout (ok for now) and also does not preallocate the struct which
%       may be a speed bottle neck as this is a low level function
%
%    Tested on: Matlab r2007b
%
%    Header changes: NONE
%
%    Usage:    [h,idx]=vinfo(data)
%
%    Examples:
%     Get the undefined value for your data:
%      h=vinfo(data);
%      h.undef
%
%    See also: seisdef, vvseis

%     Version History:
%        Sep. 25, 2008 - initial version
%        Sep. 26, 2008 - check internal header version
%        Oct. 16, 2008 - removed header version consistency check
%        Oct. 17, 2008 - filetype support, removed several pointless output
%        Oct. 25, 2008 - doc update, block mlint hint for now
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 25, 2008 at 21:30 GMT

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
        h(count)=seisdef(uft{i},uv(j)); %#ok<AGROW>
        idx(v==uv(j) & ift)=count;
    end
end

end
