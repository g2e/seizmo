function [h,idx]=versioninfo(data,varargin)
%VERSIONINFO    Returns version info for SEIZMO data records
%
%    Usage:    [h,idx]=versioninfo(data)
%              [h,idx]=versioninfo(data,'reqfield1',...,'reqfieldN')
%
%    Description:
%     [H,IDX]=VERSIONINFO(DATA) returns all necessary version definitions
%     pertaining to the records in the SEIZMO structure DATA in the struct
%     array H.  IDX has one entry per record in DATA that gives the index
%     in H that corresponds to that record's version definition. This is an
%     internal function to reduce rampant code repetition.
%
%     [H,IDX]=VERSIONINFO(DATA,'REQFIELD1',...,'REQFIELDN') passes
%     additional required SEIZMO struct fields on to SEIZMOCHECK.
%
%    Notes:
%
%    Examples:
%     % Get the undefined value for your data:
%     h=versioninfo(data);
%     h.undef
%
%    See also: SEIZMODEF, VALIDSEIZMO, SEIZMOCHECK

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
%        Sep.  7, 2009 - fixed multi-filetype bottleneck
%        Nov. 18, 2009 - allow extra input into seizmocheck
%        Jan. 29, 2010 - add export to SEIZMO and caching
%        Aug. 16, 2010 - fix error usages
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 00:05 GMT

% todo:

% check number of inputs
error(nargchk(1,inf,nargin));

% get global info
global SEIZMO

% use cache?
try
    if(SEIZMO.VERSIONINFO.USECACHE)
        h=SEIZMO.VERSIONINFO.H;
        idx=SEIZMO.VERSIONINFO.IDX;
        return;
    end
catch
end

% check data structure
error(seizmocheck(data,varargin{:}));

% get filetype and version
ft={data.filetype}.';
v=[data.version].';

% get unique filetype/versions
[u,a,idx]=unique([double(char(ft)) v],'rows');
nu=numel(a);

% loop through unique filetype/versions
h(1:nu,1)=deal(seizmodef(ft{a(1)},v(a(1))));
for i=2:nu
    h(i)=seizmodef(ft{a(i)},v(a(i)));
end

% store this info in global SEIZMO
% - very temporal (overwritten the next call)
SEIZMO.VERSIONINFO.H=h;
SEIZMO.VERSIONINFO.IDX=idx;

end
