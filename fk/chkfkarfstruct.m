function [report]=chkfkarfstruct(fk)
%CHKFKARFSTRUCT    True if is a struct as defined by FKARF
%
%    Usage:    msg=chkfkarfstruct(fk)
%
%    Description: MSG=CHKFKARFSTRUCT(FK) check that FK is a struct as
%     output from FKARF.  See FKARF for details on the struct layout.  MSG
%     is an error structure if a problem is found (otherwise it is empty).
%
%    Notes:
%
%    Examples:
%     Check that output from FKARF is ok:
%      sarf=fkarf(data,50,201,0,0,[1/30 1/20]);
%      error(chkfkarfstruct(sarf));
%
%    See also: CHKFKSTRUCT, PLOTFKARF, FKARF, FKMAP, PLOTFKMAP

%     Version History:
%        May  11, 2010 - initial version
%        May  13, 2010 - minor bug fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  13, 2010 at 01:30 GMT

% todo:

% check nargin
report=[];
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check map is proper struct
fields={'response' 'nsta' 'stla' 'stlo' 'x' 'y' 'npw' 's' 'baz' 'f' ...
    'polar' 'center' 'normdb'};
if(~isstruct(fk) || ~all(ismember(fields,fieldnames(fk))))
    report.identifier='seizmo:chkfkarfstruct:dataNotStruct';
    report.message=['FK data must be a struct with fields:\n' ...
        sprintf('%s ',fields{:})];
    return;
end

% valid center strings
valid.CENTER={'center' 'coarray' 'full'};

% loop over each frame/volume/map
for i=1:numel(fk)
    % check consistency
    sfk=size(fk(i).response);
    sx=numel(fk(i).x);
    sy=numel(fk(i).y);
    if(~isreal(fk(i).nsta) || ~isscalar(fk(i).nsta) ...
            || fk(i).nsta~=fix(fk(i).nsta) || fk(i).nsta<2 ...
            || ~isequal(size(fk(i).stla),[fk(i).nsta 1]) ...
            || ~isequal(size(fk(i).stlo),[fk(i).nsta 1]) ...
            || ~isreal(fk(i).stla) || ~isreal(fk(i).stlo))
        report.identifier='seizmo:chkfkstruct:fkStnInfoCorrupt';
        report.message='FK station info appears corrupt!';
        return;
    elseif(~isreal(fk(i).npw) || ~isscalar(fk(i).npw) ...
            || fk(i).npw~=fix(fk(i).npw) || fk(i).npw<1 ...
            || ~isequal(size(fk(i).s),[fk(i).npw 1]) ...
            || ~isequal(size(fk(i).baz),[fk(i).npw 1]) ...
            || ~isequal(size(fk(i).f),[fk(i).npw 1]) ...
            || ~isreal(fk(i).s) || ~isreal(fk(i).baz) ...
            || ~isreal(fk(i).f) || any(fk(i).f<=0))
        report.identifier='seizmo:chkfkstruct:fkPWInfoCorrupt';
        report.message='FK plane wave info appears corrupt!';
        return;
    elseif(~islogical(fk(i).polar))
        report.identifier='seizmo:chkfkstruct:polarInvalid';
        report.message='FK polar field must be TRUE/FALSE!';
        return;
    elseif((isnumeric(fk(i).center) && (~isreal(fk(i).center) ...
            || ~numel(fk(i).center)==2)) || (ischar(fk(i).center) ...
            && ~any(strcmpi(fk(i).center,valid.CENTER))))
        report.identifier='seizmo:chkfkstruct:centerInvalid';
        report.message=['CENTER must be [LAT LON], ''CENTER'', ' ...
            '''COARRAY'' or ''FULL'''];
        return;
    elseif(~isreal(fk(i).normdb) || ~isscalar(fk(i).normdb))
        report.identifier='seizmo:chkfkstruct:normdbInvalid';
        report.message='FK normdb field must be a real-valued scalar!';
        return;
    elseif(~isreal(fk(i).x) || ~isreal(fk(i).y))
        report.identifier='seizmo:chkfkstruct:xyzCorrupt';
        report.message='FK xy info appears corrupt!';
        return;
    elseif(~isreal(fk(i).response) || ~isequal(sfk,[sy sx]))
        report.identifier='seizmo:chkfkstruct:fkResponseCorrupt';
        report.message='FK response data size wrong or data invalid!';
        return;
    end
end

end
