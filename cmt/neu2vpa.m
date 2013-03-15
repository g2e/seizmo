function [v,p,a]=neu2vpa(n,e,u)
%NEU2VPA    Convert vector in North/East/Up to Value/Plunge/Azimuth
%
%    Usage:    vpa=neu2vpa(neu)
%              vpa=neu2vpa(n,e,u)
%              [v,p,a]=neu2vpa(...)
%
%    Description:
%     VPA=NEU2VPA(NEU) converts vectors given as [North East Up] in NEU to
%     [value plunge azimuth] in VPA.  Plunge is the angle down from the
%     horizontal plane and azimuth is the angle clockwise from North.  Both
%     are in degrees.  The input and output are Nx3 arrays where N allows
%     for multiple vectors to be converted simultaneously.
%
%     VPA=NEU2VPA(N,E,U) allows North, East & Up to be given separately.
%
%     [V,P,A]=NEU2VPA(...) returns value, plunge & azimuth separately.
%
%    Notes:
%     - If the vector points up then the plunge will be negative.  So all
%       normal vectors from STRIKEDIP2NEU will have negative plunge!
%
%    Examples:
%     % Vector pointing halfway between North & Down will
%     % give a plunge or 45deg and an azimuth of 0:
%     neu2vpa([1 0 -1])
%
%     % How to get a positive plunge from the normal to a fault plane:
%     neu2vpa(-1*strikedip2neu(0,30))
%
%    See also: VPA2NEU, NEU2STRIKEDIP, STRIKEDIP2NEU, SDR2TPB, TPB2SDR,
%              MT2TPB, TPB2MT

%     Version History:
%        Mar. 13, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 13, 2013 at 23:55 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% single or trip input
switch nargin
    case 1
        if(size(n,2)~=3 || ndims(n)>2)
            error('seizmo:neu2vpa:badInput',...
                'NEU must be a Nx3 array as [N E U] !');
        elseif(~isnumeric(n) || ~isreal(n))
            error('seizmo:neu2vpa:badInput',...
                'NEU must be a real-valued Nx3 array!');
        end
        [n,e,u]=deal(n(:,1),n(:,2),n(:,3));
    case 3
        if(~isnumeric(n) || ~isreal(n) || ~isnumeric(e) || ~isreal(e) ...
                || ~isnumeric(u) || ~isreal(u))
            error('seizmo:neu2vpa:badInput',...
                'N/E/U must be real-valued arrays!');
        end
        [n,e,u]=expandscalars(n,e,u);
    otherwise
        error('seizmo:neu2vpa:badNumInputs',...
            'Incorrect number of inputs (only 1 or 3)!');
end

% calculate value/plunge/azimuth using cart2sph
R2D=180/pi;
[a,p,v]=cart2sph(e,n,u);
a=R2D*(pi/2-a);
p=-R2D*p;

% combine if only one output
if(nargout<=1); v=[v(:) p(:) a(:)]; end

end
