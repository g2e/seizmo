function [Qp]=qkqu2qp(Qk,Qu,Vp,Vs)
%QKQU2QP    Calculate corresponding Qalpha for given Qk, Qu, Vp, Vs
%
%    Usage:    [Qp]=qkqu2qp(Qk,Qu,Vp,Vs)
%
%    Description:
%     [Qp]=QKQU2QP(Qk,Qu,Vp,Vs) calculates the p-wave quality factor from
%     the bulk & shear moduli quality factors and the P & S velocities.
%     This is from equation (32) on page 192 of Stein & Wysession.
%
%    Notes:
%
%    Examples:
%     % Qp for the lower crust in PREM:
%     Qk=57283; Qu=600;
%     Vp=6.8; Vs=3.9;
%     Qp=qkqu2qp(Qk,Qu,Vp,Vs)
%
%    See also: QPQS2QK

%     Version History:
%        May  19, 2010 - initial version
%        Aug.  9, 2010 - minor doc fix
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 02:55 GMT

% todo:

% check nargin
error(nargchk(4,4,nargin));

% check inputs
if(~isreal(Qk) || ~isreal(Qu) || ~isreal(Vp) || ~isreal(Vs))
    error('seizmo:qkqu2qp:badInput',...
        'Qk, Qu, Vp, Vs all must be real-valued arrays!');
elseif(any(Qk<=0) || any(Qu<=0) || any(Vp<=0) || any(Vs<0))
    error('seizmo:qkqu2qp:badInput',...
        'Qk, Qu, Vp must be > 0, Vs must be >=0!');
elseif(~isequalsizeorscalar(Qk,Qu,Vp,Vs))
    error('seizmo:qkqu2qp:badInput',...
        'Qk, Qu, Vp & Vs must be equal sized or scalar!');
end

% get Qa
L=4/3*(Vs./Vp).^2;
Qp=1./((1-L)./Qk+L./Qu);

end
