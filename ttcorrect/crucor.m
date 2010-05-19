function [corr]=crucor(lat,lon,rayp,pwave,varargin)

% todo
% options
% - elev (km above sea level - defaults to topo_points)
% - top (km below surface - defaults to 0km)
% - bottom (km below sea level - defaults to getmoho)
% - bounce (true/false - default is false, true will double correction)
% - stretch (true/false - default is true, will use elev)
% - refmod (default is prem)

% check nargin
error(nargchk(4,inf,nargin));

% check required inputs

% go through optional inputs

% check sizes & expand

% get crust for crust2.0

% get 1D crustal model

% loop over points

% get travel times

% get correction

end
