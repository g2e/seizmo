function [s]=ssidx(s,idx)
%SSIDX    Scalar struct database indexing
%
%    Usage:    s=ssidx(s,idx)
%
%    Description:
%     S=SSIDX(S,IDX) indexes into a scalar struct db S using IDX.  S is
%     expected to be a struct where every field has the same number of
%     rows (can be of any type, just needs to be NROWSx?).  S does not need
%     to be scalar (works on each struct element separately).
%
%    Notes:
%
%    Examples:
%     % Extract 8 CMT solutions from the cmt database:
%     cmt=ssidx(findcmts,33:40);
%
%    See also: SSCAT, READNDK, PARSE_ISC_ORIGIN, READSODEVENTCSV

%     Version History:
%        July 30, 2010 - initial version
%        Aug.  2, 2010 - logical indexing support
%        Mar.  7, 2011 - minor doc formatting
%        Feb. 29, 2012 - extract rows rather than elements, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 29, 2012 at 15:40 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check inputs
if(~isstruct(s))
    error('seizmo:ssidx:badInput',...
        'S must be a struct!');
end
fields=fieldnames(s);

% loop over each struct element
for i=1:numel(s)
    nrows=size(s(i).(fields{1}),1);
    for j=2:numel(fields)
        if(size(s(i).(fields{j}),1)~=nrows)
            error('seizmo:ssidx:badInput',...
                ['S(%d).%s must have the same number ' ...
                'of rows as other fields!'],...
                i,fields{j});
        end
    end
    if(~isempty(idx))
        if((islogical(idx) && isequal(numel(idx),nrows)) ...
                || (isreal(idx) && all(idx==fix(idx)) ...
                && all(idx>0) && all(idx<=nrows)))
            % good
        else
            error('seizmo:ssidx:badInput',...
                'IDX does not have valid indices!');
        end
    end
    
    % extract rows
    for j=1:numel(fields)
        s(i).(fields{j})=s(i).(fields{j})(idx,:);
    end
end

end
