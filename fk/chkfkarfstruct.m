function [report]=chkfkarfstruct(fk)
%CHKFKARFSTRUCT    Validate if a struct is as defined by FKARF
%
%    Usage:    msg=chkfkarfstruct(arf)
%
%    Description:
%     MSG=CHKFKARFSTRUCT(ARF) checks that ARF is a struct as output from
%     FKARF.  See FKARF for details on the struct layout.  MSG is an error
%     structure if a problem is found (otherwise it is empty).
%
%    Notes:
%
%    Examples:
%     % Check that output from FKARF is ok:
%     sarf=fkarf(stla,stlo,50,201,0,0,[1/30 1/20]);
%     error(chkfkarfstruct(sarf));
%
%    See also: CHKFKSTRUCT, PLOTFKARF, FKARF

%     Version History:
%        May  11, 2010 - initial version
%        May  13, 2010 - minor bug fix
%        May  24, 2010 - minor doc touch (don't forget to update Contents)
%        May  27, 2010 - fixed an error message
%        June 16, 2010 - minor code update
%        July  6, 2010 - major update for new struct
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 17:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check map is proper struct
report=[];
fields={'beam' 'nsta' 'stla' 'stlo' 'x' 'y' 'npw' 's' 'baz' 'f' ...
    'polar' 'npairs' 'method' 'center' 'normdb'};
if(~isstruct(fk) || ~all(ismember(fields,fieldnames(fk))))
    report.identifier='seizmo:chkfkarfstruct:dataNotStruct';
    report.message=['FK data must be a struct with fields:' ...
        sprintf('\n') sprintf('%s ',fields{:})];
    return;
end

% valid method strings
valid.METHOD={'center' 'coarray' 'full' 'user'};

% loop over each frame/volume/map
for i=1:numel(fk)
    % check consistency
    sfk=size(fk(i).beam);
    sx=numel(fk(i).x);
    sy=numel(fk(i).y);
    if(~isreal(fk(i).nsta) || ~isscalar(fk(i).nsta) ...
            || fk(i).nsta~=fix(fk(i).nsta) || fk(i).nsta<2 ...
            || ~isequal(size(fk(i).stla),[fk(i).nsta 1]) ...
            || ~isequal(size(fk(i).stlo),[fk(i).nsta 1]) ...
            || ~isreal(fk(i).stla) || ~isreal(fk(i).stlo))
        report.identifier='seizmo:chkfkarfstruct:fkStnInfoCorrupt';
        report.message='FK station info fields appear corrupt!';
        return;
    elseif(~isreal(fk(i).npw) || ~isscalar(fk(i).npw) ...
            || fk(i).npw~=fix(fk(i).npw) || fk(i).npw<1 ...
            || ~isequal(size(fk(i).s),[fk(i).npw 1]) ...
            || ~isequal(size(fk(i).baz),[fk(i).npw 1]) ...
            || ~isequal(size(fk(i).f),[fk(i).npw 1]) ...
            || ~isreal(fk(i).s) || ~isreal(fk(i).baz) ...
            || ~isreal(fk(i).f) || any(fk(i).f<=0))
        report.identifier='seizmo:chkfkarfstruct:fkPWInfoCorrupt';
        report.message='FK plane wave info fields appear corrupt!';
        return;
    elseif(~islogical(fk(i).polar))
        report.identifier='seizmo:chkfkstruct:polarInvalid';
        report.message='FK polar field must be TRUE/FALSE!';
        return;
    elseif(~isreal(fk(i).npairs) || fk(i).npairs~=fix(fk(i).npairs) ...
            || fk(i).npairs<0)
        report.identifier='seizmo:chkfkarfstruct:npairsInvalid';
        report.message='FK npairs field must be a positive integer!';
        return;
    elseif(~any(strcmpi(fk(i).method,valid.METHOD)))
        report.identifier='seizmo:chkfkarfstruct:methodInvalid';
        report.message=['FK method field must be ''USER'', ' ...
            '''CENTER'', ''COARRAY'' or ''FULL'''];
        return;
    elseif(isnumeric(fk(i).center) && (~isreal(fk(i).center) ...
            || ~numel(fk(i).center)==2))
        report.identifier='seizmo:chkfkarfstruct:centerInvalid';
        report.message='FK center field must be [LAT LON]!';
        return;
    elseif(~isreal(fk(i).normdb) || ~isscalar(fk(i).normdb))
        report.identifier='seizmo:chkfkarfstruct:normdbInvalid';
        report.message='FK normdb field must be a real-valued scalar!';
        return;
    elseif(~isreal(fk(i).x) || ~isreal(fk(i).y))
        report.identifier='seizmo:chkfkarfstruct:xyzCorrupt';
        report.message='FK x/y fields appear corrupt!';
        return;
    elseif(~isreal(fk(i).beam) || ~isequal(sfk,[sy sx]))
        report.identifier='seizmo:chkfkarfstruct:fkResponseCorrupt';
        report.message='FK beam field has wrong size or invalid data!';
        return;
    end
end

end
