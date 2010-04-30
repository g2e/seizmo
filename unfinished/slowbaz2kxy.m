function [kx,ky]=slowbaz2kxy(s,baz,w)
%SLOWBAZ2KXY    Converts slowness and back-azimuth to wavenumbers in x & y

baz=baz*pi/180;
kx=w.*s.*sin(baz);
ky=w.*s.*cos(baz);

end
