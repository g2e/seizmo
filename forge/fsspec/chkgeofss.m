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
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     s=geofssxc(xcdata,[lat(:) lon(:)],30,[.01 .0125]);
%     error(chkgeofss(s));
%
%    See also: CHKFSS, GEOFSSXC, GEOFSSHORZ, GEOFSSHORZXC, ISGEOFSS

%     Version History:
%        June 22, 2010 - initial version
%        June 24, 2010 - bugfix
%        June 25, 2010 - more lenient to allow a variety of volume types
%        July  6, 2010 - major update for new struct
%        Apr.  4, 2012 - minor doc update
%        June  4, 2012 - adapted from chkgeofkstruct
%        June 10, 2012 - handle full method
%        June 13, 2012 - allow capon method
%        Sep. 12, 2012 - minor doc update
%        Sep. 29, 2012 - added new methods, drop vector field, add whiten
%                        & weights, combine st fields
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 18:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check map is proper struct
report=[];
fields={'spectra' 'nsta' 'st' 'butc' 'eutc' 'npts' 'delta' 'latlon' ...
        'slow' 'freq' 'npairs' 'method' 'center' 'whiten' 'weights'};
if(~isstruct(s) || ~all(ismember(fields,fieldnames(s))))
    report.identifier='seizmo:chkgeofss:dataNotStruct';
    report.message=['GEOFSS data must be a struct with fields:' ...
        sprintf('\n') sprintf('%s ',fields{:})];
    return;
end

% valid method strings
valid.METHOD={'center' 'coarray' 'user' 'full' ...
    'capon' 'xc' 'caponxc' 'tdxc'};

% loop over each spectra
for i=1:numel(s)
    % check consistency
    ss=size(s(i).spectra);
    if(numel(ss)==2); ss(3)=1; end
    sx=size(s(i).latlon,1);
    sy=numel(s(i).slow);
    sf=numel(s(i).freq);
    if(~isreal(s(i).nsta) || ~isscalar(s(i).nsta) ...
            || s(i).nsta~=fix(s(i).nsta) || s(i).nsta<2 ...
            || size(s(i).st,1)~=s(i).nsta || size(s(i).st,2)<2 ...
            || ~isreal(s(i).st))
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
        report.message='GEOFSS method field contains an invalid method!';
        return;
    elseif(~isnumeric(s(i).center) || (~isreal(s(i).center) ...
            || numel(s(i).center)<2))
        report.identifier='seizmo:chkgeofss:centerInvalid';
        report.message='GEOFSS center field must be [LAT LON]!';
        return;
    elseif(~isreal(s(i).latlon) || size(s(i).latlon,2)<2 ...
            || ndims(s(i).latlon)~=2 ...
            || (~isreal(s(i).slow) ...
            && ~isa(s(i).slow,'function_handle')) ...
            || (isreal(s(i).slow) && any(s(i).slow<=0)) ...
            || ~isreal(s(i).freq) || any(s(i).freq<0))
        report.identifier='seizmo:chkgeofss:xyzCorrupt';
        report.message='GEOFSS position/slow/freq info is corrupt!';
        return;
    elseif(~isreal(s(i).spectra) ...
            || (ss(3)>1 && ss(2)>1 && ~isequal(ss,[sx sy sf])) ...
            || (ss(3)>1 && ss(2)==1 && ~isequal(ss,[sx 1 sf])) ...
            || (ss(3)==1 && ss(2)>1 && ~isequal(ss,[sx sy 1])) ...
            || (ss(3)==1 && ss(2)==1 && ~isequal(ss,[sx 1 1])))
        report.identifier='seizmo:chkgeofss:spectraCorrupt';
        report.message='GEOFSS spectra has wrong size or invalid data!';
        return;
    end
end

end
