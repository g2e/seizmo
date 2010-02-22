function [topo,lat,lon]=srtm30plus_region(latrange,lonrange)
%SRTM30PLUS_REGION    Returns SRTM30plus data for the region indicated

% check nargin
msg=nargchk(0,2,nargin);
if(~isempty(msg)); error(msg); end

% check lat/lon

% take care of wrap-around

% what regions do we load

% load regions

% combine

% get lat/lon (are we grid/cell referenced?)

% truncate?

end