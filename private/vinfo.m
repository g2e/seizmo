function [h,vi,v,uv,nv]=vinfo(data)
%VINFO    Returns version info for SAClab data records
%
%    Description: VINFO(DATA) returns an abundance of version related info
%     pertaining to the records in the SAClab structure DATA.  This is an
%     internal function to reduce rampant code repetition.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: DATA must be a valid SAClab structure
%
%    Header changes: N/A
%
%    Usage:    [h,vi,v,uv,nv]=vinfo(data)
%
%    Examples:
%     NONE
%
%    See also: get_n_check

%     Version History:
%        Sep. 25, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 25, 2008 at 06:45 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data))

% get versions
v=[data.version].';
uv=unique(v);
nv=length(uv);

% grab definitions
vi=nan(size(v));
h(nv)=seisdef(uv(nv));
vi(v==uv(nv))=nv;
for i=1:nv-1
    h(i)=seisdef(uv(i));
    vi(v==uv(i))=i;
end

end