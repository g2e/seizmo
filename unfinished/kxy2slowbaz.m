function [s,baz]=kxy2slowbaz(kx,ky,w)
%KXY2SLOWBAZ    Converts wavenumbers in x & y to slowness and back-azimuth

s=sqrt(kx.^2+ky.^2)./w;
baz=atan2(kx,ky)*180/pi;

end
