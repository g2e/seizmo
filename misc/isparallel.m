function [lgc]=isparallel(o1,o2)
%ISPARALLEL    TRUE if orientations are parallel
%
%    Usage:    lgc=isparallel(o1,o2)
%
%    Description: LGC=ISPARALLEL(O1,O2) returns TRUE for each pair of
%     orientations formed by O1 & O2 that are parallel.  O1 & O2 must be
%     Nx2 real arrays where the first column gives the inclination from
%     vertical and the second column gives the azimuth from North.  Must be
%     in degrees!  If either O1 or O2 is a scalar orientation (ie 1x2) then
%     it is expanded to the size of the other orientation array.
%
%    Notes:
%     - O1 & O2 must be in DEGREES!
%
%    Examples:
%     Some simple negative examples:
%      isparallel([  0   0],[ 90   0]) % vertical vs north
%      isparallel([  0   0],[ 90  90]) % vertical vs east
%      isparallel([ 90   0],[ 90  90]) % north vs east
%
%     Some wrap-around examples:
%      isparallel([  0   0],[  0 360])
%      isparallel([ 90  90],[ 90 450])
%      isparallel([ 90 -90],[ 90 270])
%
%    See also: ISORTHOGONAL, DOT

%     Version History:
%        Nov.  2, 2009 - initial version
%        Feb. 22, 2010 - bug fix for size checks, single-prec tolerance
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 22, 2010 at 21:35 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check orientations
sz1=size(o1); sz2=size(o2);
if(~isreal(o1) || sz1(2)~=2)
    error('seizmo:isparallel:badInput',...
        'O1 must be a Nx2 matrix of reals!');
end
if(~isreal(o2) || sz2(2)~=2)
    error('seizmo:isparallel:badInput',...
        'O2 must be a Nx2 matrix of reals!');
end

% tile scalar & check size matches
if(prod(sz1([1 3:end]))==1)
    o1=repmat(o1,[sz2(1) 1 sz2(3:end)]);
    sz1=sz2;
elseif(prod(sz2([1 3:end]))==1)
    o2=repmat(o2,[sz1(1) 1 sz1(3:end)]);
    sz2=sz1;
end
if(~isequal(sz1,sz2))
    error('seizmo:isparallel:badInput',...
        'O1 & O2 must be single orientations or equal sized!');
end

% permute
nd=ndims(o1);
o1=permute(o1,[2 1 3:nd]);
o2=permute(o2,[2 1 3:nd]);

% convert to radians then to unit vector in xyz
o1=o1*pi/180; o2=o2*pi/180;
o1=cat(1,sin(o1(2,:)).*sin(o1(1,:)),...
    cos(o1(2,:)).*sin(o1(1,:)),cos(o1(1,:)));
o2=cat(1,sin(o2(2,:)).*sin(o2(1,:)),...
    cos(o2(2,:)).*sin(o2(1,:)),cos(o2(1,:)));

% check parallel
lgc=ipermute(sum(abs(o1-o2)<(2*pi*eps(single(1))),1)==3,[2 1 3:nd]);

end
