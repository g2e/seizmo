function delta=mkxflat2sfer(xf,rp);
% mkxflat2sfer........transform epicentral distance from flat to spherical world
%
% call: delta=mkxflat2sfer(xf,rp);
%
%       xf: epicentral; distance in Flat world [km]
%       rp: planetary radius [km]
%
% result: delta: epicentral distance in spherical world [deg]
%
% when transforming ray geometry from flat back to spherical, 
% you have to transfrom depths prior to distances!
%
% Martin Knapmeyer, 02.05.2002



% transform km into deg - kmfactor contains the rad-deg-conversion
kmfactor=360./(2*pi*rp); % on earth surface: 1deg per 111.19km
delta=xf.*kmfactor;

