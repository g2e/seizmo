function [lgc]=isorthogonal(o1,o2)
%ISORTHOGONAL    TRUE if orientations are orthogonal
%
%    Usage:    lgc=isorthogonal(o1,o2)
%
%    Description: LGC=ISORTHOGONAL(O1,O2) returns TRUE/FALSE for each pair
%     of orientations formed by O1 & O2.  O1 & O2 must be Nx2 real arrays
%     where the first column gives the inclination from vertical and the
%     second column gives the azimuth from North.  Must be in degrees!  If
%     either O1 or O2 is a scalar orientation (ie 1x2) then it is expanded
%     to the size of the other orientation array.  
%
%    Notes:
%     - O1 & O2 must be in DEGREES!
%
%    Examples:
%     Some simple affirmative examples:
%      isorthogonal([ 0  0],[90  0]) % vertical vs north
%      isorthogonal([ 0  0],[90 90]) % vertical vs east
%      isorthogonal([90  0],[90 90]) % north vs east
%
%     Some wrap-around examples:
%      isorthogonal([  0   0],[270   0])
%      isorthogonal([ 90  90],[ 90 360])
%      isorthogonal([ 90 -90],[ 90 180])
%
%    See also: ISPARALLEL, DOT

%     Version History:
%        Nov.  2, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  2, 2009 at 19:15 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check orientations
sz1=size(o1); sz2=size(o2);
if(~isreal(o1) || sz1(2)~=2)
    error('seizmo:isorthogonal:badInput',...
        'O1 must be a Nx2 matrix of reals!');
end
if(~isreal(o2) || sz2(2)~=2)
    error('seizmo:isorthogonal:badInput',...
        'O2 must be a Nx2 matrix of reals!');
end

% tile scalar & check size matches
if(prod(sz1([1 3:end]))==1)
    o1=repmat(o1,[sz2(1) 1 sz2(3:end)]);
elseif(prod(sz2([1 3:end]))==1)
    o2=repmat(o2,[sz1(1) 1 sz1(3:end)]);
end
if(~isequal(sz1([1 3:end]),sz1([1 3:end])))
    error('seizmo:isorthogonal:badInput',...
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

% check orthogonality
lgc=ipermute(abs(sum(o1.*o2,1))<3*2*pi*eps,[2 1 3:nd]);

end
