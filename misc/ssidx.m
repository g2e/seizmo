function [s]=ssidx(s,idx)
%SSIDX    Scalar struct database indexing
%
%    Usage:    s=ssidx(s,idx)
%
%    Description:
%     S=SSIDX(S,IDX) indexes into a scalar struct db S using IDX.  S is
%     expected to be a struct where every field has the same number of
%     elements (can be of any type).  This will likely require char array
%     fields to be converted to cellstr arrays.  S does not need to be
%     scalar.
%
%    Notes:
%
%    Examples:
%     %Extract 8 CMT solutions from the cmt database:
%     cmt=load('globalcmt_full.mat');
%     cmt=ssidx(cmt,33:40);
%
%    See also: READNDK, PARSE_ISC_ORIGIN

%     Version History:
%        July 30, 2010 - initial version
%        Aug.  2, 2010 - logical indexing support
%        Mar.  7, 2011 - minor doc formatting
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  7, 2011 at 15:40 GMT

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
    maxidx=numel(s(i).(fields{1}));
    for j=2:numel(fields)
        if(numel(s(i).(fields{j}))~=maxidx)
            error('seizmo:ssidx:badInput',...
                'S(%d).%s must have the same numel as other fields!',...
                i,fields{j});
        end
    end
    if(~isempty(idx))
        if((islogical(idx) && isequal(numel(idx),maxidx)) ...
                || (isreal(idx) && all(idx==fix(idx)) ...
                && all(idx>0) && all(idx<=maxidx)))
            % good
        else
            error('seizmo:ssidx:badInput',...
                'IDX does not have valid indices!');
        end
    end
    
    % extract
    for j=1:numel(fields)
        s(i).(fields{j})=s(i).(fields{j})(idx);
    end
end

end
