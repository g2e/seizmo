function [u]=radpat(mt,inc,az,type)
%RADPAT    Calculates moment tensor radiation pattern
%
%    Usage:    u=radpat(mt,inc,az)
%              u=radpat(mt,inc,az,type)
%
%    Description:
%     U=RADPAT(MT,INC,AZ) calculates the dimensionless displacement for a
%     moment tensor given by MT at the angles INC & AZ.  MT must be in
%     Harvard convention with one of the following formats:
%        - 1x6  [Mrr Mtt Mpp Mrt Mrp Mtp]
%        - 3x3  [Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp]
%        - scalar struct with fields .mrr .mtt .mpp .mrt .mrp .mtp
%     INC gives the angle from down (ie. the inclination) and AZ gives the
%     clockwise angle from North.  U contains the P, SV, & SH displacements
%     as U(:,1), U(:,2), & U(:,3).
%
%     U=RADPAT(MT,THETA,PHI,TYPE) calculates the displacement for the
%     specified wavetype TYPE.  TYPE must be either 'P', 'SH', 'SV' or
%     'ALL'.  The default is 'ALL'.
%
%    Notes:
%     - Calculations are based on Aki & Richards (2002) eqs 4.88 & 4.96
%     - INC is related to the rayparameter (s/deg), source radius (km)
%       & velocity (km/s) by:  INC=asind(180/pi*P*V0/R)
%
%    Examples:
%     % Get the dimensionless P & S wave displacements at the T/P/B axes
%     % of the lowest magnitude moment tensor from the GlobalCMT project:
%     cmt=mt_change('g',findcmt('magnitude',4));
%     [t,p,b]=mt2tpb(cmt);
%     radpat(cmt,90-[t(2) p(2) b(2)],[t(3) p(3) b(3)])
%
%    See also: FINDCMT, FINDCMTS, PLOTMT, PLOTMT3, MAPCMTS, RAYP2INC

%     Version History:
%        Mar.  8, 2011 - initial version modified from Ken Creager's Coral
%        Mar. 11, 2011 - improved struct checking
%        Feb.  6, 2012 - minor doc update, handle single couple (no norm)
%        Mar. 25, 2013 - theta & phi renamed to inc & az for clarity, minor
%                        doc update
%
%     Written by Ken Creager (kcc+ess/washington/edu)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 23:55 GMT

% todo:

% check nargin
error(nargchk(3,4,nargin));

% default type to all
if(nargin==3 || isempty(type)); type='all'; end

% check inputs
error(mt_check(mt));
mt=hrv2ar(mt_change('g',mt));
if(size(mt,3)~=1)
    error('seizmo:radpat:badInput',...
        'MT must be only contain one moment tensor for RADPAT!');
end
if(~isnumeric(inc) || ~isreal(inc) || any(abs(inc(:))>180))
    error('seizmo:radpat:badInput',...
        ['INC must be an array of real-valued angles ' ...
        'from 0 (down) to 180 (up)!']);
end
if(~isnumeric(az) || ~isreal(az))
    error('seizmo:radpat:badInput',...
        ['AZ must be an array of real-valued azimuth angles ' ...
        '(0 = North, 90 = East)!']);
elseif(~isequalsizeorscalar(inc,az))
    error('seizmo:radpat:badInput',...
        'INC & AZ must be scalar or equally sized!');
end
if(~ischar(type) || size(type,1)~=1 ...
        || ~any(strcmpi(type,{'p' 'sh' 'sv' 'all'})))
    error('seizmo:radpat:badInput',...
        'TYPE must be ''p'' ''sh'' ''sv'' or ''all''!');
end

% expand scalar inc/az & force to be column vectors
[inc,az]=expandscalars(inc,az);
inc=inc(:); az=az(:);

% normalize mt
% - handle moment=0 case (single force?)
mo=sqrt(trace(mt*mt)/2);
if(mo); mt=mt/mo; end

% compute p/sv/sh directions
gamma=[sind(inc).*cosd(az) sind(inc).*sind(az) cosd(inc)];
sv=[cosd(inc).*cosd(az) cosd(inc).*sind(az) -sind(inc)];
sh=[-sind(az) cosd(az) zeros(size(az))];

% get displacements
u=nan(numel(inc),1);
switch lower(type)
    case 'p'
        for i=1:numel(inc); u(i,1)=gamma(i,:)*mt*gamma(i,:).'; end
    case 'sv'
        for i=1:numel(inc); u(i,1)=sv(i,:)*mt*gamma(i,:).'; end
    case 'sh'
        for i=1:numel(inc); u(i,1)=sh(i,:)*mt*gamma(i,:).'; end
    case 'all'
        u=[u u u];
        for i=1:numel(inc)
            u(i,1)=gamma(i,:)*mt*gamma(i,:).';
            u(i,2)=sv(i,:)*mt*gamma(i,:).';
            u(i,3)=sh(i,:)*mt*gamma(i,:).';
        end
end

end

