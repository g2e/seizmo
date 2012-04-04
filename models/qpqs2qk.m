function [Qk]=qpqs2qk(Qp,Qs,Vp,Vs)
%QPQS2QK    Calculate corresponding Qkappa for given Qp, Qs, Vp, Vs
%
%    Usage:    [Qk]=qpqs2qk(Qp,Qs,Vp,Vs)
%
%    Description:
%     [Qk]=QPQS2QK(Qp,Qs,Vp,Vs) calculates the bulk modulus (eg. 
%     incompressibility) quality factor from the P & S wave quality factors
%     and velocities.  This is from equation (32) on page 192 of Stein &
%     Wysession.
%
%    Notes:
%
%    Examples:
%     % Qk for the lower crust in PREM:
%     Qp=1350; Qs=600;
%     Vp=6.8; Vs=3.9;
%     Qk=qpqs2qk(Qp,Qs,Vp,Vs)
%
%    See also: QKQU2QP

%     Version History:
%        May  19, 2010 - initial version
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 02:55 GMT

% todo:

% check nargin
error(nargchk(4,4,nargin));

% check inputs
if(~isreal(Qp) || ~isreal(Qs) || ~isreal(Vp) || ~isreal(Vs))
    error('seizmo:qpqs2qk:badInput',...
        'Qp, Qs, Vp, Vs all must be real-valued arrays!');
elseif(any(Qp<=0) || any(Qs<=0) || any(Vp<=0) || any(Vs<0))
    error('seizmo:qpqs2qk:badInput',...
        'Qp, Qs, Vp must be > 0, Vs must be >=0!');
elseif(~isequalsizeorscalar(Qp,Qs,Vp,Vs))
    error('seizmo:qpqs2qk:badInput',...
        'Qp, Qs, Vp & Vs must be equal sized or scalar!');
end

% get Qk
L=4/3*(Vs./Vp).^2;
Qk=(1-L)./(1./Qp-L./Qs);

end
