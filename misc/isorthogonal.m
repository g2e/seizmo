function [lgc]=isorthogonal(o1,o2)
%ISORTHOGONAL    TRUE if orientations are orthogonal
%
%    Usage:    lgc=isorthogonal(o1,o2)
%
%    Description:
%     LGC=ISORTHOGONAL(O1,O2) returns TRUE/FALSE for each pair of
%     orientations formed by O1 & O2.  O1 & O2 must be Nx2 real arrays
%     where the first column gives the inclination from vertical and the
%     second column gives the azimuth from North.  Must be in degrees!  If
%     either O1 or O2 is a scalar orientation (ie 1x2) then it is expanded
%     to the size of the other orientation array.  
%
%    Notes:
%     - O1 & O2 must be in DEGREES!
%
%    Examples:
%     % Some simple affirmative examples:
%     isorthogonal([ 0  0],[90  0]) % vertical vs north
%     isorthogonal([ 0  0],[90 90]) % vertical vs east
%     isorthogonal([90  0],[90 90]) % north vs east
%
%     % Some wrap-around examples:
%     isorthogonal([  0   0],[270   0])
%     isorthogonal([ 90  90],[ 90 360])
%     isorthogonal([ 90 -90],[ 90 180])
%
%    See also: ISPARALLEL, DOT

%     Version History:
%        Nov.  2, 2009 - initial version
%        Feb. 22, 2010 - bug fix for size checks, single-prec tolerance
%        Feb. 11, 2011 - mass nargchk fix
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

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
    sz1=sz2;
elseif(prod(sz2([1 3:end]))==1)
    o2=repmat(o2,[sz1(1) 1 sz1(3:end)]);
    sz2=sz1;
end
if(~isequal(sz1,sz2))
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
lgc=ipermute(abs(sum(o1.*o2,1))<3*2*pi*eps(single(1)),[2 1 3:nd]);

end
