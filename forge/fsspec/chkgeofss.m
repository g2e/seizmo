function [report]=chkgeofss(s)
%CHKGEOFSS    Validate if a struct is as defined by GEOFSS functions
%
%    Usage:    msg=chkgeofss(s)
%
%    Description:
%     MSG=CHKGEOFSS(S) checks that S is a struct as output from GEOFSS
%     functions.  See those functions for details on the struct layout.
%     MSG is an error structure if a problem is found (otherwise it is
%     empty).
%
%    Notes:
%
%    Examples:
%     % Check that output from GEOFSSXC is compatible:
%     s=geofssxc(xcdata,latlon,slow,frng);
%     error(chkgeofss(s));
%
%    See also: CHKFSS, GEOFSSXC, GEOFSSHORZXC, PLOTGEOFSS, ISGEOFSS

%     Version History:
%        June 22, 2010 - initial version
%        June 24, 2010 - bugfix
%        June 25, 2010 - more lenient to allow a variety of volume types
%        July  6, 2010 - major update for new struct
%        Apr.  4, 2012 - minor doc update
%        June  4, 2012 - adapted from chkgeofkstruct
%        June 10, 2012 - handle full method
%        June 13, 2012 - allow capon method
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 13, 2012 at 18:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check map is proper struct
report=[];
fields={'spectra' 'nsta' 'stla' 'stlo' 'stel' 'stdp' 'butc' 'eutc' ...
        'npts' 'delta' 'latlon' 'slow' 'freq' 'npairs' 'method' ...
        'center' 'vector'};
if(~isstruct(s) || ~all(ismember(fields,fieldnames(s))))
    report.identifier='seizmo:chkgeofss:dataNotStruct';
    report.message=['GEOFSS data must be a struct with fields:' ...
        sprintf('\n') sprintf('%s ',fields{:})];
    return;
end

% valid method strings
valid.METHOD={'center' 'coarray' 'user' 'full' 'capon'};

% loop over each frame/volume/map
for i=1:numel(s)
    % check consistency
    ss=size(s(i).spectra);
    sx=size(s(i).latlon,1);
    sy=numel(s(i).slow);
    sf=numel(s(i).freq);
    if(~isreal(s(i).nsta) || ~isscalar(s(i).nsta) ...
            || s(i).nsta~=fix(s(i).nsta) || s(i).nsta<2 ...
            || ~isequal(size(s(i).stla),[s(i).nsta 1]) ...
            || ~isequal(size(s(i).stlo),[s(i).nsta 1]) ...
            || ~isequal(size(s(i).stel),[s(i).nsta 1]) ...
            || ~isequal(size(s(i).stdp),[s(i).nsta 1]) ...
            || ~isreal(s(i).stla) || ~isreal(s(i).stdp) ...
            || ~isreal(s(i).stlo) || ~isreal(s(i).stel))
        report.identifier='seizmo:chkgeofss:fkStnInfoCorrupt';
        report.message='GEOFSS station info fields appear corrupt!';
        return;
    elseif(~isequal(size(s(i).butc),[1 5]) || ~isreal(s(i).butc) ...
            || ~isequal(size(s(i).eutc),[1 5]) || ~isreal(s(i).eutc))
        report.identifier='seizmo:chkgeofss:fkTimingInfoCorrupt';
        report.message='GEOFSS timing info appears corrupt!';
        return;
    elseif(~isreal(s(i).npts) || ~isscalar(s(i).npts) ...
            || s(i).npts~=fix(s(i).npts) || s(i).npts<0 ...
            || ~isreal(s(i).delta) || ~isscalar(s(i).delta) ...
            || s(i).delta<=0)
        report.identifier='seizmo:chkgeofss:fkWaveformInfoCorrupt';
        report.message='GEOFSS npts/delta fields appear corrupt!';
        return;
    elseif(~isreal(s(i).npairs) ...
            || s(i).npairs~=fix(s(i).npairs) || s(i).npairs<0)
        report.identifier='seizmo:chkgeofss:npairsInvalid';
        report.message='GEOFSS npairs field must be a positive integer!';
        return;
    elseif(~any(strcmpi(s(i).method,valid.METHOD)))
        report.identifier='seizmo:chkgeofss:methodInvalid';
        report.message=['GEOFSS method field must be ''USER'', ' ...
            '''CENTER'', or ''COARRAY'''];
        return;
    elseif(~isnumeric(s(i).center) || (~isreal(s(i).center) ...
            || ~numel(s(i).center)==2))
        report.identifier='seizmo:chkgeofss:centerInvalid';
        report.message='GEOFSS center field must be [LAT LON]!';
    elseif(~islogical(s(i).vector) || ~isequal(size(s(i).vector),[1 2]))
        report.identifier='seizmo:chkgeofss:vectorInvalid';
        report.message='GEOFSS vector field must be 1x2 logical array!';
        return;
    elseif(~isreal(s(i).latlon) || size(s(i).latlon,2)~=2 ...
            || ndims(s(i).latlon)~=2 || ~isreal(s(i).slow) ...
            || any(s(i).slow<=0) ...
            || ~isreal(s(i).freq) || any(s(i).freq<0))
        report.identifier='seizmo:chkgeofss:xyzCorrupt';
        report.message='GEOFSS position/slow/freq info is corrupt!';
        return;
    elseif(~isnumeric(s(i).spectra) || (all(s(i).vector) ...
            && ((sf~=1 && ~isequal(ss,[sx sy sf])) ...
            || (sf==1 && ~isequal(ss,[sx sy])))) ...
            || (isequal(s(i).vector,[false true]) ...
            && (~isequal(ss,[sx sy]))) ...
            || (isequal(s(i).vector,[true false]) ...
            && ((sf~=1 && ~isequal(ss,[sx 1 sf])) ...
            || (sf==1 && ~isequal(ss,[sx 1])))) ...
            || (isequal(s(i).vector,[false false]) ...
            && ~isequal(ss,[sx 1])))
        report.identifier='seizmo:chkgeofss:spectraCorrupt';
        report.message='GEOFSS spectra has wrong size or invalid data!';
        return;
    end
end

end
