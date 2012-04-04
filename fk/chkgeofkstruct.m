function [report]=chkgeofkstruct(geofk)
%CHKGEOFKSTRUCT    Validate if a struct is as defined by GEOFK functions
%
%    Usage:    msg=chkgeofkstruct(geofk)
%
%    Description:
%     MSG=CHKGEOFKSTRUCT(GEOFK) checks that GEOFK is a struct as output
%     from GEOFK functions.  See those functions for details on the struct
%     layout.  MSG is an error structure if a problem is found (otherwise
%     it is empty).
%
%    Notes:
%
%    Examples:
%     % Check that output from GEOFKXCVOLUME is compatible:
%     sgeo=geofkxcvolume(xcdata,latlon,horzslow,freqrng);
%     error(chkgeofkstruct(sgeo));
%
%    See also: CHKFKSTRUCT, GEOFKXCVOLUME, GEOFKXCHORZVOLUME, PLOTGEOFKMAP

%     Version History:
%        June 22, 2010 - initial version
%        June 24, 2010 - bugfix
%        June 25, 2010 - more lenient to allow a variety of volume types
%        July  6, 2010 - major update for new struct
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 18:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check map is proper struct
report=[];
fields={'beam' 'nsta' 'stla' 'stlo' 'stel' 'stdp' 'butc' 'eutc' 'npts' ...
        'delta' 'latlon' 'horzslow' 'freq' 'npairs' 'method' 'center' ...
        'normdb' 'volume'};
if(~isstruct(geofk) || ~all(ismember(fields,fieldnames(geofk))))
    report.identifier='seizmo:chkgeofkstruct:dataNotStruct';
    report.message=['FK data must be a struct with fields:' ...
        sprintf('\n') sprintf('%s ',fields{:})];
    return;
end

% valid method strings
valid.METHOD={'center' 'coarray' 'full' 'user'};

% loop over each frame/volume/map
for i=1:numel(geofk)
    % check consistency
    sfk=size(geofk(i).beam);
    sx=size(geofk(i).latlon,1);
    sy=numel(geofk(i).horzslow);
    sf=numel(geofk(i).freq);
    if(~isreal(geofk(i).nsta) || ~isscalar(geofk(i).nsta) ...
            || geofk(i).nsta~=fix(geofk(i).nsta) || geofk(i).nsta<2 ...
            || ~isequal(size(geofk(i).stla),[geofk(i).nsta 1]) ...
            || ~isequal(size(geofk(i).stlo),[geofk(i).nsta 1]) ...
            || ~isequal(size(geofk(i).stel),[geofk(i).nsta 1]) ...
            || ~isequal(size(geofk(i).stdp),[geofk(i).nsta 1]) ...
            || ~isreal(geofk(i).stla) || ~isreal(geofk(i).stdp) ...
            || ~isreal(geofk(i).stlo) || ~isreal(geofk(i).stel))
        report.identifier='seizmo:chkgeofkstruct:fkStnInfoCorrupt';
        report.message='GEOFK station info fields appear corrupt!';
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
        report.message='GEOFK npts/delta fields appear corrupt!';
        return;
    elseif(~isreal(geofk(i).npairs) ...
            || geofk(i).npairs~=fix(geofk(i).npairs) || geofk(i).npairs<0)
        report.identifier='seizmo:chkgeofkstruct:npairsInvalid';
        report.message='GEOFK npairs field must be a positive integer!';
        return;
    elseif(~any(strcmpi(geofk(i).method,valid.METHOD)))
        report.identifier='seizmo:chkgeofkstruct:methodInvalid';
        report.message=['GEOFK method field must be ''USER'', ' ...
            '''CENTER'', ''COARRAY'' or ''FULL'''];
        return;
    elseif(isnumeric(geofk(i).center) && (~isreal(geofk(i).center) ...
            || ~numel(geofk(i).center)==2))
        report.identifier='seizmo:chkgeofkstruct:centerInvalid';
        report.message='GEOFK center field must be [LAT LON]!';
    elseif(~islogical(geofk(i).volume) ...
            || ~isequal(size(geofk(i).volume),[1 2]))
        report.identifier='seizmo:chkgeofkstruct:volumeInvalid';
        report.message='GEOFK volume field must be 1x2 logical array!';
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
    elseif(~isreal(geofk(i).beam) ...
            || (all(geofk(i).volume) ...
            && ((sf~=1 && ~isequal(sfk,[sx sy sf])) ...
            || (sf==1 && ~isequal(sfk,[sx sy])))) ...
            || (isequal(geofk(i).volume,[true false]) ...
            && (~isequal(sfk,[sx sy]))) ...
            || (isequal(geofk(i).volume,[false true]) ...
            && ((sf~=1 && ~isequal(sfk,[sx 1 sf])) ...
            || (sf==1 && ~isequal(sfk,[sx 1])))) ...
            || (isequal(geofk(i).volume,[false false]) ...
            && ~isequal(sfk,[sx 1])))
        report.identifier='seizmo:chkgeofkstruct:fkResponseCorrupt';
        report.message='GEOFK beam field has wrong size or invalid data!';
        return;
    end
end

end
