function [h,idx]=versioninfo(data)
%VERSIONINFO    Returns version info for SEIZMO data records
%
%    Usage:    [h,idx]=versioninfo(data)
%
%    Description: [H,IDX]=VERSIONINFO(DATA) returns all necessary version
%     definitions pertaining to the records in the SEIZMO structure DATA in
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
%    Examples:
%     Get the undefined value for your data:
%      h=versioninfo(data);
%      h.undef
%
%    See also: seizmodef, validseizmo

%     Version History:
%        Sep. 25, 2008 - initial version
%        Sep. 26, 2008 - check internal header version
%        Oct. 16, 2008 - removed header version consistency check
%        Oct. 17, 2008 - filetype support, removed several pointless output
%        Oct. 25, 2008 - doc update, block mlint hint for now
%        Nov. 13, 2008 - renamed from VINFO to VERSIONINFO
%        Nov. 15, 2008 - update for new name schema
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 21:15 GMT

% todo:

% check number of inputs
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

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
        h(count)=seizmodef(uft{i},uv(j)); %#ok<AGROW>
        idx(v==uv(j) & ift)=count;
    end
end

end
