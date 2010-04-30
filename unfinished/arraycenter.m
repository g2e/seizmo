function [clat,clon]=arraycenter(lat,lon)
%ARRAYCENTER    Returns the center of an array on a sphere

[x,y,z]=geocentric2xyz(lat,lon,1);
n=numel(x);
[clat,clon]=xyz2geocentric(sum(x(:))/n,sum(y(:))/n,sum(z(:))/n);

end
