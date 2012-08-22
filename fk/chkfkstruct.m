function [report]=chkfkstruct(fk)
%CHKFKSTRUCT    Validate if a struct is as defined by FKMAP/VOLUME/4D
%
%    Usage:    msg=chkfkstruct(fk)
%
%    Description:
%     MSG=CHKFKSTRUCT(FK) checks that FK is a struct as output from FKMAP,
%     FKVOLUME, or FK4D.  See those functions for details on the struct
%     layout.  MSG is an error structure if a problem is found (otherwise
%     it is empty).
%
%    Notes:
%
%    Examples:
%     % Check that output from FKMAP is compatible:
%     smap=fkmap(data,50,201,[1/30 1/20]);
%     error(chkfkstruct(smap));
%
%    See also: CHKFKARFSTRUCT, PLOTFKMAP, FKMAP, FKVOLUME, FK4D

%     Version History:
%        May  11, 2010 - initial version
%        May  12, 2010 - allow single frequency volume
%        May  24, 2010 - minor doc update (don't forget to update Contents)
%        May  27, 2010 - fixed an error message
%        June 16, 2010 - minor code update
%        July  6, 2010 - major update for new struct
%        Apr.  4, 2012 - minor doc update
%        June 13, 2012 - add capon method support
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 13, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check map is proper struct
report=[];
fields={'beam' 'nsta' 'stla' 'stlo' 'stel' 'stdp' 'butc' 'eutc' 'npts' ...
    'delta' 'x' 'y' 'freq' 'polar' 'npairs' 'method' 'center' 'normdb' ...
    'volume'};
if(~isstruct(fk) || ~all(ismember(fields,fieldnames(fk))))
    report.identifier='seizmo:chkfkstruct:dataNotStruct';
    report.message=['FK data must be a struct with fields:' ...
        sprintf('\n') sprintf('%s ',fields{:})];
    return;
end

% valid method strings
valid.METHOD={'center' 'coarray' 'full' 'user' 'capon'};

% loop over each frame/volume/map
for i=1:numel(fk)
    % check consistency
    sfk=size(fk(i).beam);
    sx=numel(fk(i).x);
    sy=numel(fk(i).y);
    sf=numel(fk(i).freq);
    if(~isreal(fk(i).nsta) || ~isscalar(fk(i).nsta) ...
            || fk(i).nsta~=fix(fk(i).nsta) || fk(i).nsta<2 ...
            || ~isequal(size(fk(i).stla),[fk(i).nsta 1]) ...
            || ~isequal(size(fk(i).stlo),[fk(i).nsta 1]) ...
            || ~isequal(size(fk(i).stel),[fk(i).nsta 1]) ...
            || ~isequal(size(fk(i).stdp),[fk(i).nsta 1]) ...
            || ~isreal(fk(i).stla) || ~isreal(fk(i).stdp) ...
            || ~isreal(fk(i).stlo) || ~isreal(fk(i).stel))
        report.identifier='seizmo:chkfkstruct:fkStnInfoCorrupt';
        report.message='FK station info fields appear corrupt!';
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
        report.message='FK npts/delta fields appear corrupt!';
        return;
    elseif(~islogical(fk(i).polar))
        report.identifier='seizmo:chkfkstruct:polarInvalid';
        report.message='FK polar field must be TRUE/FALSE!';
        return;
    elseif(~isreal(fk(i).npairs) || fk(i).npairs~=fix(fk(i).npairs) ...
            || fk(i).npairs<0)
        report.identifier='seizmo:chkfkstruct:npairsInvalid';
        report.message='FK npairs field must be a positive integer!';
        return;
    elseif(~any(strcmpi(fk(i).method,valid.METHOD)))
        report.identifier='seizmo:chkfkstruct:methodInvalid';
        report.message=['FK method field must be ''USER'', ' ...
            '''CENTER'', ''COARRAY'' or ''FULL'''];
        return;
    elseif(isnumeric(fk(i).center) && (~isreal(fk(i).center) ...
            || ~numel(fk(i).center)==2))
        report.identifier='seizmo:chkfkstruct:centerInvalid';
        report.message='FK center field must be [LAT LON]!';
        return;
    elseif(~islogical(fk(i).volume))
        report.identifier='seizmo:chkfkstruct:volumeInvalid';
        report.message='FK volume field must be TRUE/FALSE!';
        return;
    elseif(~isreal(fk(i).normdb) || ~isscalar(fk(i).normdb))
        report.identifier='seizmo:chkfkstruct:normdbInvalid';
        report.message='FK normdb field must be a real-valued scalar!';
        return;
    elseif(~isreal(fk(i).x) || ~isreal(fk(i).y) || ~isreal(fk(i).freq) ...
            || any(fk(i).freq<0))
        report.identifier='seizmo:chkfkstruct:xyzCorrupt';
        report.message='FK x/y/freq fields appear corrupt!';
        return;
    elseif(~isreal(fk(i).beam) ...
            || (fk(i).volume && sf~=1 && ~isequal(sfk,[sy sx sf])) ...
            || (sf==1 && ~isequal(sfk,[sy sx])))
        report.identifier='seizmo:chkfkstruct:fkResponseCorrupt';
        report.message='FK beam field has wrong size or invalid data!';
        return;
    end
end

end
