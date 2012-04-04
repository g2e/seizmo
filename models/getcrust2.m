function [crust]=getcrust2(lat,lon)
%GETCRUST2    Returns Crust2.0 info at specified location(s)
%
%    Usage:    crust=getcrust2(lat,lon)
%
%    Description:
%     CRUST=GETCRUST2(LAT,LON) returns struct CRUST containing info about
%     the crustal structure at the positions given by LAT & LON.  LAT & LON
%     are expected to be scalars or equal sized arrays and have units in
%     degrees.  CRUST has a layout as follows:
%      CRUST.type        -- 2 character crust id
%           .description -- string describing structure
%           .vp          -- isotropic p-velocities (1x8)
%           .vs          -- isotropic s-velocities (1x8)
%           .rho         -- densities (1x8)
%           .thick       -- layer kilometer thicknesses (1x7)
%           .elev        -- elevation (in meters)
%           .moho        -- depth of moho from sealevel (in km)
%     The vp/vs/rho fields contain info about the layers:
%      ice water softseds hardseds uppercrust middlecrust lowercrust moho
%     in that order.  The thick field contains info about all but the Moho
%     layer.
%
%    Notes:
%     - ICE    WATER SOFTSED HARDSED UPPERCRUST MIDDLECRUST LOWERCRUST MOHO
%       vp(1)   vp(2)  vp(3)   vp(4)   vp(5)      vp(6)       vp(7)   vp(8)
%       vs(1)   vs(2)  vs(3)   vs(4)   vs(5)      vs(6)       vs(7)   vs(8)
%       rho(1) rho(2) rho(3)  rho(4)  rho(5)     rho(6)      rho(7)  rho(8)
%       thk(1) thk(2) thk(3)  thk(4)  thk(5)     thk(6)      thk(7)    ---
%     - Please do not use the thickness of the water layer given by
%       CRUST.thick(2).  It is an average for that crust type.  Use the
%       elevation given in CRUST.elev for local bathymetry.
%
%    Examples:
%     % What is the crustal structure below Saint Louis?
%     getcrust2(38.649,-90.305)
%
%    See also: GETC2MOHO, GETC2THICK, GETC2ELEV, CRUCOR, MANCOR

%     Version History:
%        May  16, 2010 - initial version
%        May  17, 2010 - improved saved crust2, added moho
%        May  19, 2010 - added moho to doc
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 18:35 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check lat/lon
if(~isreal(lat) || ~isreal(lon))
    error('seizmo:getcrust2:badInput',...
        'LAT & LON must be real-valued arrays!');
elseif(~isequalsizeorscalar(lat,lon))
    error('seizmo:getcrust2:badInput',...
        'LAT & LON must be equal sized or scalar!');
end

% make the same size
[lat,lon]=expandscalars(lat,lon);

% make sure lat/lon values are in proper ranges
[lat,lon]=fixlatlon(lat,lon);

% get dimensions for output
spts=size(lat);
npts=numel(lat);

% load crust2.0
mod=load('crust2');

% lat/lon to i/j to idx
i=fix((90-lat)/2)+1;
i(i==91)=90;
j=fix((lon+180)/2)+1;
j(j==181)=180;
idx=(j-1)*90+i;

% get crust
[crust(1:npts,1)]=deal(mod.key(mod.type(idx(:))));
elev=num2cell(mod.elev(idx));
moho=num2cell(mod.moho(idx));
[crust.elev]=deal(elev{:});
[crust.moho]=deal(moho{:});
crust=reshape(crust,spts);

end
