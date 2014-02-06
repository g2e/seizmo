function [s]=ww3cat(s,varargin)
%WW3CAT    Concatenates WW3 hindcast data from WW3STRUCT
%
%    Usage:    s=ww3cat(s)
%              s=ww3cat(s1,s2,...)
%
%    Description:
%     S=WW3CAT(S) concatenates the data in each element of S into one
%     element.  So if the input S is 1x3 the output is 1x1 with the data
%     and time fields concatenated and sorted by time.  This requires that
%     the elements all have the same data type & lat/lon grid.
%
%     S=WW3CAT(S1,S2,...) concatenates across the same indices of S1, S2,
%     etc.  This operation could be done in a loop with the 1st usage form
%     like this:
%      for i=1:numel(S1)
%          S(i)=WW3CAT([S1(i) S2(i) ...]);
%      end
%     All structs S1, S2, etc should be the same size and should be
%     concatenatible across indices.
%
%    Notes:
%
%    Examples:
%     % Combine 2 months of wave height data:
%     s=ww3struct({'*hs*200501*grb' '*hs*200502*grb'});
%     s=ww3cat(s);
%
%    See also: WW3STRUCT, WW3REC, PLOTWW3, PLOTWW3TS, WW3MOV, WW3MAP,
%              WW3MAPMOV, WW3UV2SA, WW3BAZ2AZ

%     Version History:
%        May  11, 2012 - initial version
%        Jan. 15, 2014 - updated See also list
%        Feb.  5, 2014 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  5, 2014 at 00:40 GMT

% todo:

% check number of inputs
error(nargchk(1,inf,nargin));

% concatenate struct elements or separate structs
if(nargin==1)
    % elements
    error(chkww3(s));
    if(isequal(s.lat) && isequal(s.lon) ...
            && isequal(s.description) && isequal(s.units))
        % concatenate
        s(1).time=cat(2,s.time);
        s(1).data=catnindex(s.data);
        s(2:end)=[];
        
        % sort by time
        [s.time,idx]=sort(s.time);
        for i=1:numel(s.data); s.data{i}=s.data{i}(:,:,idx); end
    else
        error('seizmo:ww3cat:badInput',...
            'S elements must have the same [LAT LON] grid & datatype!');
    end
else % nargin>1
    % structs
    error(chkww3(s));
    for i=2:nargin; error(chkww3(varargin{i-1})); end
    [s,varargin{:}]=expandscalars(s,varargin{:});
    sz=size(s);
    s=cat(numel(sz)+1,s,varargin{:});
    s=permute(s,[numel(sz)+1 1:numel(sz)]);
    for i=1:prod(sz)
        s(1,i)=ww3cat(s(:,i));
    end
    s=ipermute(s,[numel(sz)+1 1:numel(sz)]);
    c={':'}; c=c(ones(1,numel(sz)));
    s(c{:},2:nargin)=[];
end

end

function [d]=catnindex(varargin)
    d=cat(3,varargin{:});
    for i=1:size(d,1)
        for j=1:size(d,2)
            d{i,j,1}=cat(3,d{i,j,:});
            [d{i,j,2:end}]=deal([]);
        end
    end
    d(:,:,2:end)=[];
end

function [report]=chkww3(s)
report=[];
if(~isstruct(s) || any(~isfield(s,...
        {'description' 'units' 'lat' 'lon' 'data' 'time'})))
    report.identifier='seizmo:chkww3:badStruct';
    report.message='S must be a struct from WW3STRUCT!';
end
end
