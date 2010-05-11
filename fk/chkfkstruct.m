function [report]=chkfkstruct(fk)
%CHKFKSTRUCT    True if is a struct as defined by FKMAP/VOLUME/4D
%
%    Usage:    msg=chkfkstruct(fk)
%
%    Description: MSG=CHKFKSTRUCT(FK) check that FK is a struct as output
%     from FKMAP/FKVOLUME/FK4D.  See those functions for details on the
%     struct layout.  MSG is an error structure if a problem is found
%     (otherwise it is empty).
%
%    Notes:
%
%    Examples:
%     Check that output from FKMAP is compatible:
%      smap=fkmap(data,50,201,[1/30 1/20]);
%      error(chkfkstruct(smap));
%
%    See also: PLOTFKMAP, FKMAP, FKVOLUME, FK4D

%     Version History:
%        May  11, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  11, 2010 at 01:00 GMT

% todo:

% check nargin
report=[];
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check map is proper struct
fields={'response' 'nsta' 'stla' 'stlo' 'stel' 'stdp' 'butc' 'eutc' ...
    'npts' 'delta' 'x' 'y' 'z' 'polar' 'center' 'normdb' 'volume'};
if(~isstruct(fk) || ~all(ismember(fields,fieldnames(fk))))
    report.identifier='seizmo:chkfkstruct:dataNotStruct';
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
    sz=numel(fk(i).z);
    if(~isreal(fk(i).nsta) || ~isscalar(fk(i).nsta) ...
            || fk(i).nsta~=fix(fk(i).nsta) || fk(i).nsta<2 ...
            || ~isequal(size(fk(i).stla),[fk(i).nsta 1]) ...
            || ~isequal(size(fk(i).stlo),[fk(i).nsta 1]) ...
            || ~isequal(size(fk(i).stel),[fk(i).nsta 1]) ...
            || ~isequal(size(fk(i).stdp),[fk(i).nsta 1]) ...
            || ~isreal(fk(i).stla) || ~isreal(fk(i).stdp) ...
            || ~isreal(fk(i).stlo) || ~isreal(fk(i).stel))
        report.identifier='seizmo:chkfkstruct:fkStnInfoCorrupt';
        report.message='FK station info appears corrupt!';
        return;
    elseif(~isequal(size(fk(i).butc),[1 5]) || ~isreal(fk(i).butc) ...
            || ~isequal(size(fk(i).eutc),[1 5]) || ~isreal(fk(i).eutc))
        report.identifier='seizmo:chkfkstruct:fkTimingInfoCorrupt';
        report.message='FK timing info appears corrupt!';
        return;
    elseif(~isreal(fk(i).npts) || ~isscalar(fk(i).npts) ...
            || fk(i).npts~=fix(fk(i).npts) || fk(i).npts<0 ...
            || ~isreal(fk(i).delta) || ~isscalar(fk(i).delta) ...
            || fk(i).delta<=0)
        report.identifier='seizmo:chkfkstruct:fkWaveformInfoCorrupt';
        report.message='FK waveform info appears corrupt!';
        return;
    elseif(~islogical(fk(i).polar))
        report.identifier='seizmo:chkfkstruct:polarInvalid';
        report.message='FK polar field must be TRUE/FALSE!';
        return;
    elseif(~islogical(fk(i).volume))
        report.identifier='seizmo:chkfkstruct:volumeInvalid';
        report.message='FK volume field must be TRUE/FALSE!';
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
    elseif(~isreal(fk(i).x) || ~isreal(fk(i).y) || ~isreal(fk(i).z) ...
            || any(fk(i).z<0))
        report.identifier='seizmo:chkfkstruct:xyzCorrupt';
        report.message='FK xyz info appears corrupt!';
        return;
    elseif(~isreal(fk(i).response) ...
            || (fk(i).volume && ~isequal(sfk,[sy sx sz])) ...
            || (~fk(i).volume && ~isequal(sfk,[sy sx])))
        report.identifier='seizmo:chkfkstruct:fkResponseCorrupt';
        report.message='FK response data size wrong or data invalid!';
        return;
    end
end

end
