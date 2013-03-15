function [strike,dip,rake]=auxplane(strike,dip,rake)
%AUXPLANE    Returns strike-dip-slip of 2nd (auxiliary) focal plane
%
%    Usage:  sdr=auxplane(sdr)
%            sdr=auxplane(strike,dip,rake)
%            [strike,dip,rake]=auxplane(...)
%
%    Description:
%     SDR=AUXPLANE(SDR) finds the auxiliary plane(s) and rake(s) that match
%     the strain(s) produced by the input plane(s) and rake(s).  SDR must
%     be given as an Nx3 array arranged as [strike dip rake] where N allows
%     for multiple planes & rakes to be processed simultaneously.  Strike
%     is clockwise from North, dip is positive downward from the horizontal
%     and rake is counter-clockwise in the fault plane from the strike
%     direction.  Note that the strike must be such that when you look
%     along the direction of the strike the fault dips to your right.
%
%     SDR=AUXPLANE(STRIKE,DIP,RAKE) allows strike, dip & rake to be given
%     separately.
%
%     [STRIKE,DIP,RAKE]=AUXPLANE(...) returns strike, dip & rake
%     separately.
%
%    Notes:
%     - The auxiliary plane is perpendicular to the input fault plane.
%     - The normal & slip vector of the auxiliary plane is the slip &
%       normal vector of the fault plane
%     - The plane perpendicular to both the fault plane and the auxiliary
%       plane has its normal along the null axis of the focal mechanism.
%
%    Examples:
%     % Corresponding auxiliary fault plane and rake for
%     % dip-slip on a normal fault that strikes North:
%     auxplane(0,45,-90)
%
%    See also: NEU2STRIKEDIP, STRIKEDIP2NEU, SDR2SLIP, MT2SDR, SDR2MT

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 22, 2010 - added docs
%        June  1, 2011 - updated docs
%        Mar. 14, 2013 - rewrite
%        Mar. 15, 2013 - more cleanup
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% one or both inputs
switch nargin
    case 1
        if(size(strike,2)~=3 || ndims(strike)>2)
            error('seizmo:auxplane:badInput',...
                'SDR must be a Nx3 array as [STRIKE DIP RAKE] !');
        elseif(~isnumeric(strike) || ~isreal(strike))
            error('seizmo:auxplane:badInput',...
                'SDR must be a real-valued Nx3 array!');
        end
        [strike,dip,rake]=deal(strike(:,1),strike(:,2),strike(:,3));
    case 3
        if(~isnumeric(strike) || ~isreal(strike) ...
                || ~isnumeric(dip) || ~isreal(dip) ...
                || ~isnumeric(rake) || ~isreal(rake))
            error('seizmo:auxplane:badInput',...
                'STRIKE/DIP/RAKE must be real-valued arrays!');
        end
        [strike,dip,rake]=expandscalars(strike,dip,rake);
end

% get strike/dip/rake of auxiliary plane
% -> sdr to normal, slip
% -> switch & convert back to sdr
normal=strikedip2norm(strike,dip);
slip=sdr2slip(strike,dip,rake);
[strike,dip,rake]=normslip2sdr(slip,normal);

%{
% make strike relative to west (why?)
strike=strike+90;

% cos & sin
scos=cosd(strike); ssin=sind(strike);
dcos=cosd(dip);    dsin=sind(dip);
rcos=cosd(rake);   rsin=sind(rake);

% vector pointing in direction of slip
% is normal to the auxiliary plane
sl1=-rcos.*scos-rsin.*ssin.*dcos;
sl2=rcos.*ssin-rsin.*scos.*dcos;
sl3=rsin.*dsin;

% strike & dip of auxiliary plane
[strike,dip]=neu2strikedip(sl2,sl1,sl3);

% rake direction is normal to the primary focal plane
% - getting normal to primary focal plane
n1=ssin.*dsin;
n2=scos.*dsin;
%n3=dcos;

% a vector along intersection of aux. plane with horizontal
h1=-sl2;
h2=sl1;
%h3=0*sl1;

% angle between normal vector and horizontal vector gives rake (+/-)
rake=acosd((h1.*n1+h2.*n2)./sqrt(h1.^2+h2.^2));
    
% fix sense of orientation
j=find(sl3<=0);
rake(j)=-rake(j);
%}

% combine if only one output
if(nargout<=1); strike=[strike(:) dip(:) rake(:)]; end

end
