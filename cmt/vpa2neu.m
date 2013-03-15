function [n,e,u]=vpa2neu(v,p,a)
%VPA2NEU    Convert vector in Value/Plunge/Azimuth to North/East/Up
%
%    Usage:    neu=vpa2neu(vpa)
%              neu=vpa2neu(v,p,a)
%              [n,e,u]=vpa2neu(...)
%
%    Description:
%     NEU=VPA2NEU(VPA) converts vectors given as [value plunge azimuth] in
%     in VPA to [North East Up] in NEU.  Plunge is the angle down from the
%     horizontal plane and azimuth is the angle clockwise from North.  Both
%     are in degrees.  The input and output are Nx3 arrays where N allows
%     for multiple vectors to be converted simultaneously.
%
%     NEU=VPA2NEU(V,P,A) allows value, plunge & azimuth to be given
%     separately.
%
%     [N,E,U]=VPA2NEU(...) returns North, East & Up separately.
%
%    Notes:
%
%    Examples:
%     % Vector with an azimuth of 45deg and plunging at 45deg will
%     % not give equal magnitude North, East, & Up components:
%     vpa2neu([1 45 45])
%
%    See also: NEU2VPA, NEU2STRIKEDIP, STRIKEDIP2NEU, SDR2TPB, TPB2SDR,
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
        if(size(v,2)~=3 || ndims(v)>2)
            error('seizmo:vpa2neu:badInput',...
                'VPA must be a Nx3 array as [V P A] !');
        elseif(~isnumeric(v) || ~isreal(v))
            error('seizmo:vpa2neu:badInput',...
                'VPA must be a real-valued Nx3 array!');
        end
        [v,p,a]=deal(v(:,1),v(:,2),v(:,3));
    case 3
        if(~isnumeric(v) || ~isreal(v) || ~isnumeric(p) || ~isreal(p) ...
                || ~isnumeric(a) || ~isreal(a))
            error('seizmo:vpa2neu:badInput',...
                'V/P/A must be real-valued arrays!');
        end
        [v,p,a]=expandscalars(v,p,a);
    otherwise
        error('seizmo:vpa2neu:badNumInputs',...
            'Incorrect number of inputs (only 1 or 3)!');
end

% calculate north/east/up using sph2cart
D2R=pi/180;
[e,n,u]=sph2cart(pi/2-D2R*a,-D2R*p,v);

% combine if only one output
if(nargout<=1); n=[n(:) e(:) u(:)]; end

end
