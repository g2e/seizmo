function [report]=chkgeofkstruct(geofk)
%CHKGEOFKSTRUCT    Validate if a struct is as defined by GEOFK functions
%
%    Usage:    msg=chkgeofkstruct(geofk)
%
%    Description: MSG=CHKGEOFKSTRUCT(GEOFK) checks that GEOFK is a struct
%     as output from GEOFK functions.  See those functions for details on
%     the struct layout.  MSG is an error structure if a problem is found
%     (otherwise it is empty).
%
%    Notes:
%
%    Examples:
%     Check that output from GEOFKXCVOLUME is compatible:
%      sgeo=geofkxcvolume(xcdata,latlon,horzslow,freqrng);
%      error(chkgeofkstruct(sgeo));
%
%    See also: CHKFKSTRUCT, CHKGEOFKSTRUCT, GEOFKXCVOLUME

%     Version History:
%        June 22, 2010 - initial version
%        June 24, 2010 - bugfix
%        June 25, 2010 - more lenient to allow a variety of volume types
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 25, 2010 at 11:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check map is proper struct
report=[];
fields={'response' 'nsta' 'stla' 'stlo' 'stel' 'stdp' 'butc' 'eutc' ...
    'npts' 'delta' 'latlon' 'horzslow' 'freq' 'normdb' 'volume'};
if(~isstruct(geofk) || ~all(ismember(fields,fieldnames(geofk))))
    report.identifier='seizmo:chkgeofkstruct:dataNotStruct';
    report.message=['FK data must be a struct with fields:' ...
        sprintf('\n') sprintf('%s ',fields{:})];
    return;
end

% loop over each frame/volume/map
for i=1:numel(geofk)
    % check consistency
    sfk=size(geofk(i).response);
    sx=size(geofk(i).latlon,1);
    sy=numel(geofk(i).horzslow);
    sz=numel(geofk(i).freq);
    if(~isreal(geofk(i).nsta) || ~isscalar(geofk(i).nsta) ...
            || geofk(i).nsta~=fix(geofk(i).nsta) || geofk(i).nsta<2 ...
            || ~isequal(size(geofk(i).stla),[geofk(i).nsta 1]) ...
            || ~isequal(size(geofk(i).stlo),[geofk(i).nsta 1]) ...
            || ~isequal(size(geofk(i).stel),[geofk(i).nsta 1]) ...
            || ~isequal(size(geofk(i).stdp),[geofk(i).nsta 1]) ...
            || ~isreal(geofk(i).stla) || ~isreal(geofk(i).stdp) ...
            || ~isreal(geofk(i).stlo) || ~isreal(geofk(i).stel))
        report.identifier='seizmo:chkgeofkstruct:fkStnInfoCorrupt';
        report.message='GEOFK station info appears corrupt!';
        return;
    elseif(~isequal(size(geofk(i).butc),[1 5]) ...
            || ~isreal(geofk(i).butc) ...
            || ~isequal(size(geofk(i).eutc),[1 5]) ...
            || ~isreal(geofk(i).eutc))
        report.identifier='seizmo:chkgeofkstruct:fkTimingInfoCorrupt';
        report.message='GEOFK timing info appears corrupt!';
        return;
    elseif(~isreal(geofk(i).npts) || ~isscalar(geofk(i).npts) ...
            || geofk(i).npts~=fix(geofk(i).npts) || geofk(i).npts<0 ...
            || ~isreal(geofk(i).delta) || ~isscalar(geofk(i).delta) ...
            || geofk(i).delta<=0)
        report.identifier='seizmo:chkgeofkstruct:fkWaveformInfoCorrupt';
        report.message='GEOFK waveform info appears corrupt!';
        return;
    elseif(~islogical(geofk(i).volume) ...
            || ~isequal(size(geofk(i).volume),[1 2]))
        report.identifier='seizmo:chkgeofkstruct:volumeInvalid';
        report.message='GEOFK volume field must be TRUE/FALSE!';
        return;
    elseif(~isreal(geofk(i).normdb) || ~isscalar(geofk(i).normdb))
        report.identifier='seizmo:chkgeofkstruct:normdbInvalid';
        report.message='GEOFK normdb field must be a real-valued scalar!';
        return;
    elseif(~isreal(geofk(i).latlon) || size(geofk(i).latlon,2)~=2 ...
            || ndims(geofk(i).latlon)~=2 || ~isreal(geofk(i).horzslow) ...
            || any(geofk(i).horzslow<=0) ...
            || ~isreal(geofk(i).freq) || any(geofk(i).freq<0))
        report.identifier='seizmo:chkgeofkstruct:xyzCorrupt';
        report.message='GEOFK position/horzslow/freq info is corrupt!';
        return;
    elseif(~isreal(geofk(i).response) ...
            || (all(geofk(i).volume) ...
            && ((sz~=1 && ~isequal(sfk,[sx sy sz])) ...
            || (sz==1 && ~isequal(sfk,[sx sy])))) ...
            || (isequal(geofk(i).volume,[true false]) ...
            && (~isequal(sfk,[sx sy]))) ...
            || (isequal(geofk(i).volume,[false true]) ...
            && ((sz~=1 && ~isequal(sfk,[sx 1 sz])) ...
            || (sz==1 && ~isequal(sfk,[sx 1])))) ...
            || (isequal(geofk(i).volume,[false false]) ...
            && ~isequal(sfk,[sx 1])))
        report.identifier='seizmo:chkgeofkstruct:fkResponseCorrupt';
        report.message='GEOFK response data size wrong or data invalid!';
        return;
    end
end

end
