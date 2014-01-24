function [paths]=extract_upswing_raypaths(paths,depth)
%EXTRACT_UPSWING_RAYPATHS    Returns the "upswing" portion of raypaths
%
%    Usage:    paths=extract_upswing_raypaths(paths,depth)
%
%    Description:
%     PATHS=EXTRACT_UPSWING_RAYPATHS(PATHS,DEPTH) returns the "upswing"
%     portion of raypaths in PATHS.  The "upswing" is just the last section
%     of the path on the receiver side that is above the cutoff depth
%     DEPTH.  This is primarily for getting core-diffracted mantle
%     corrections for azimuthal profiles.  DEPTH must be in kilometers &
%     PATHS is expected to conform to the TAUPPATH format.
%
%    Notes:
%
%    Examples:
%     % Construct a large group of diffracted upswing paths and plot them:
%     evla=31.5; evlo=140.07; evdp=35;
%     [stla,stlo]=sphericalfwd(evla,evlo,...
%                              90+70*rand(100,1),360*rand(100,1));
%     paths=getraypaths('P,Pdiff','prem',evla,evlo,evdp,stla,stlo);
%     paths=extract_upswing_raypaths(paths,2891-500);
%     plotraypaths(paths);
%
%    See also: GETRAYPATHS, CRUSTLESS_RAYPATHS, MANCOR, PLOTRAYPATHS,
%              INSERT_DEPTHS_IN_RAYPATHS, TAUPPATH, TRIM_DEPTHS_RAYPATHS

%     Version History:
%        May  31, 2010 - initial version
%        June  3, 2010 - updated the example, name changed from GET_*,
%                        reduced number of inputs by making path generation
%                        external
%        Feb. 27, 2012 - update for tauppath changes
%        Jan. 23, 2014 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 02:45 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% input should be a tauppath struct
test=tauppath('ph','P','ev',[0 0],'st',[0 10]);
if(~isstruct(paths) || any(~ismember(fieldnames(paths),fieldnames(test))))
    error('seizmo:extract_upswing_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember(fieldnames(paths(1).path),fieldnames(test(1).path))))
    error('seizmo:extract_upswing_raypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
end

% specific locations or not?
% - this fails if some are and some are not
if(any(~ismember({'latitude' 'longitude'},fieldnames(paths(1).path))))
    specific=false;
else % not specific
    specific=true;
end

% check the depth
if(~isreal(depth) || ~isequalsizeorscalar(paths,depth))
    error('seizmo:extract_upswing_raypaths:badDEPTH',...
        'DEPTH must be a real-valued scalar or array (1 depth per path)!');
end
if(isscalar(depth)); depth=depth(ones(numel(paths),1),1); end

% loop over each path extracting the upswing portion above depth
for i=1:numel(paths)
    % last point in path below upswing depth
    idx=find(paths(i).path.depth>depth(i),1,'last');
    
    % go to next if did not reach upswing depth
    if(isempty(idx)); continue; end
    
    % get fraction of segment crossing upswing depth below
    frac=(depth(i)-paths(i).path.depth(idx))...
        /(diff(paths(i).path.depth(idx:idx+1)));
    
    % remove all before it
    fields=fieldnames(paths(i).path);
    for j=1:numel(fields)
        paths(i).path.(fields{j})(1:(idx-1))=[];
    end
    
    % shift values of last point below to values at upswing depth
    % - assuming that the values are linear between the points
    % - we handle lat/lon specially so no problems there
    % - could do better but honestly not worth it yet
    paths(i).path.depth(1)=depth(i);
    paths(i).path.time(1)=paths(i).path.time(1)...
        +frac*diff(paths(i).path.time(1:2));
    paths(i).path.distance(1)=paths(i).path.distance(1)...
        +frac*diff(paths(i).path.distance(1:2));
    if(specific)
        [dist,az]=sphericalinv(paths(i).path.latitude(1),...
            paths(i).path.longitude(1),...
            paths(i).path.latitude(2),paths(i).path.longitude(2));
        [paths(i).path.latitude(1),...
            paths(i).path.longitude(1)]=sphericalfwd(...
            paths(i).path.latitude(1),...
            paths(i).path.longitude(1),dist*frac,az);
    end
end

end
