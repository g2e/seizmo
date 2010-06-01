function [paths]=get_upswing_raypaths(depth,varargin)
%GET_UPSWING_RAYPATHS    Returns the "upswing" portion of seismic raypaths
%
%    Usage:    paths=get_upswing_raypaths(depth,phase,mod,...
%                                      evla,evlo,evdp,stla,stlo)
%
%    Description:
%     PATHS=GET_UPSWING_RAYPATHS(DEPTH,PHASE,MOD,EVLA,EVLO,EVDP,STLA,STLO)
%     returns the "upswing" portion of a ray path (which is just the last
%     section of the path on the receiver side that is above the cutoff
%     depth DEPTH).  This is mainly for getting core-diffracted mantle
%     corrections.  Depths must be in kilometers, lat/lons must be in
%     degrees.  PHASE & MOD must be recognizable by TAUPPATH!
%
%    Notes:
%
%    Examples:
%     There really is nothing worth putting here yet...
%
%    See also: GETRAYPATHS, CRUST2LESS_RAYPATHS, MANCOR, TAUPPATH

%     Version History:
%        May  31, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  31, 2010 at 20:05 GMT

% todo:

% check nargin
error(nargchk(8,8,nargin));

% pass to get_phase_paths
paths=getraypaths(varargin{:});

% check the depth
if(~isreal(depth) || ~isequalsizeorscalar(paths,depth))
    error('seizmo:get_upswing_raypaths:badDEPTH',...
        ['DEPTH must be a real valued & must be scalar ' ...
        'or equal-sized with the other inputs!']);
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
    paths(i).path.rayparameter(1)=paths(i).path.rayparameter(1)...
        +frac*diff(paths(i).path.rayparameter(1:2));
    paths(i).path.time(1)=paths(i).path.time(1)...
        +frac*diff(paths(i).path.time(1:2));
    paths(i).path.distance(1)=paths(i).path.distance(1)...
        +frac*diff(paths(i).path.distance(1:2));
    [dist,az]=sphericalinv(...
        paths(i).path.latitude(1),paths(i).path.longitude(1),...
        paths(i).path.latitude(2),paths(i).path.longitude(2));
    [paths(i).path.latitude(1),paths(i).path.longitude(1)]=sphericalfwd(...
        paths(i).path.latitude(1),paths(i).path.longitude(1),dist*frac,az);
end

end
