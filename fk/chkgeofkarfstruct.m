function [report]=chkgeofkarfstruct(geofk)
%CHKGEOFKARFSTRUCT    Validate if a struct is as defined by GEOFKARF
%
%    Usage:    msg=chkgeofkarfstruct(arf)
%
%    Description:
%     MSG=CHKGEOFKARFSTRUCT(ARF) checks that ARF is a struct as output from
%     GEOFKARF.  See GEOFKARF for details on the struct layout.  MSG is an
%     error structure if a problem is found (otherwise it is empty).
%
%    Notes:
%
%    Examples:
%     % Check that output from GEOFKARF is ok:
%     lat0=-20:20;
%     lon0=-20:20;
%     [lat,lon]=meshgrid(lat0,lon0);
%     arf=geofkarf([stla stlo],[lat(:) lon(:)],27.5:0.5:32.5,...
%                  [5 10],30,1/26.3,'center');
%     error(chkgeofkarfstruct(arf));
%
%    See also: CHKGEOFKSTRUCT, GEOFKARF, PLOTGEOFKARF, GEOFKSUBARF,
%              UPDATEGEOFKARF, GEOFKARFSLOWSLIDE, GEOFKARF2MAP

%     Version History:
%        July  7, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 18:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check map is proper struct
report=[];
fields={'beam' 'nsta' 'stla' 'stlo' 'latlon' 'horzslow' 'nsw' 'latlon0' ...
        'horzslow0' 'freq0' 'npairs' 'method' 'center' 'normdb' 'volume'};
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
    sf=numel(geofk(i).freq0); % not right but meh ...
    if(~isreal(geofk(i).nsta) || ~isscalar(geofk(i).nsta) ...
            || geofk(i).nsta~=fix(geofk(i).nsta) || geofk(i).nsta<2 ...
            || ~isequal(size(geofk(i).stla),[geofk(i).nsta 1]) ...
            || ~isequal(size(geofk(i).stlo),[geofk(i).nsta 1]) ...
            || ~isreal(geofk(i).stla) || ~isreal(geofk(i).stlo))
        report.identifier='seizmo:chkgeofkarfstruct:fkStnInfoCorrupt';
        report.message='ARF station info fields appear corrupt!';
        return;
    elseif(~isreal(geofk(i).nsw) || ~isscalar(geofk(i).nsw) ...
            || geofk(i).nsw~=fix(geofk(i).nsw) || geofk(i).nsw<1 ...
            || ~isequal(size(geofk(i).horzslow0),[geofk(i).nsw 1]) ...
            || ~isequal(size(geofk(i).latlon0),[geofk(i).nsw 2]) ...
            || ~isequal(size(geofk(i).freq0),[geofk(i).nsw 1]) ...
            || ~isreal(geofk(i).horzslow0) || ~isreal(geofk(i).latlon0) ...
            || any(geofk(i).horzslow0<=0) ...
            || ~isreal(geofk(i).freq0) || any(geofk(i).freq0<=0))
        report.identifier='seizmo:chkfkarfstruct:fkPWInfoCorrupt';
        report.message='ARF spherical wave info fields appear corrupt!';
        return;
    elseif(~isreal(geofk(i).npairs) ...
            || geofk(i).npairs~=fix(geofk(i).npairs) || geofk(i).npairs<0)
        report.identifier='seizmo:chkgeofkstruct:npairsInvalid';
        report.message='ARF npairs field must be a positive integer!';
        return;
    elseif(~any(strcmpi(geofk(i).method,valid.METHOD)))
        report.identifier='seizmo:chkgeofkstruct:methodInvalid';
        report.message=['ARF method field must be ''USER'', ' ...
            '''CENTER'', ''COARRAY'' or ''FULL'''];
        return;
    elseif(isnumeric(geofk(i).center) && (~isreal(geofk(i).center) ...
            || ~numel(geofk(i).center)==2))
        report.identifier='seizmo:chkgeofkstruct:centerInvalid';
        report.message='ARF center field must be [LAT LON]!';
    elseif(~islogical(geofk(i).volume) ...
            || ~isequal(size(geofk(i).volume),[1 2]))
        report.identifier='seizmo:chkgeofkstruct:volumeInvalid';
        report.message='ARF volume field must be 1x2 logical array!';
        return;
    elseif(~isreal(geofk(i).normdb) || ~isscalar(geofk(i).normdb))
        report.identifier='seizmo:chkgeofkstruct:normdbInvalid';
        report.message='ARF normdb field must be a real-valued scalar!';
        return;
    elseif(~isreal(geofk(i).latlon) || size(geofk(i).latlon,2)~=2 ...
            || ndims(geofk(i).latlon)~=2 || ~isreal(geofk(i).horzslow) ...
            || any(geofk(i).horzslow<=0))
        report.identifier='seizmo:chkgeofkstruct:xyzCorrupt';
        report.message='ARF latlon/horzslow fields are corrupt!';
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
        report.message='ARF beam field has wrong size or invalid data!';
        return;
    end
end

end
