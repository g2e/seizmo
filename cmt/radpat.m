function [u]=radpat(mt,theta,phi,type)
%RADPAT    Calculates moment tensor radiation pattern
%
%    Usage:    u=radpat(mt,theta,phi)
%              u=radpat(mt,theta,phi,type)
%
%    Description:
%     U=RADPAT(MT,THETA,PHI) calculates the dimensionless displacement
%     for a moment tensor given by MT at the angles given by THETA & PHI.
%     MT must be in Harvard convention and is one of the following formats:
%       1x6  [Mrr Mtt Mpp Mrt Mrp Mtp]
%       3x3  [Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp]
%       scalar struct with scalar fields .mrr .mtt .mpp .mrt .mrp .mtp
%     THETA & PHI give the angle from down (ie the inclination) and the
%     azimuth (from North) respectively.  U contains the P, SV, & SH
%     displacements as U(:,1), U(:,2), & U(:,3).
%
%     U=RADPAT(MT,THETA,PHI,TYPE) calculates the displacement for the
%     specified wavetype TYPE.  TYPE must be either 'P', 'SH', 'SV' or
%     'ALL'.  The default is 'ALL'.
%
%    Notes:
%     - Calculations are based on Aki & Richards (2002) eqs 4.88 & 4.96
%     - THETA is related to the rayparameter, source depth & velocity by:
%        THETA=asind(P*V0/R0)
%
%    Examples:
%     % Grab a cmt and get the P-wave displacement pattern:
%     [theta,phi]=meshgrid(0:180,0:360);
%     cmt=findcmt();
%     u_p=radpat(cmt,theta,phi,'p');
%     figure; s=surface(x,y,z,'facecolor','texturemap','cdata',u_p');
%     axis vis3d
%
%    See also: FINDCMT, BB

%     Version History:
%        Mar.  8, 2011 - initial version
%
%     Written by Ken Creager (kcc+ess/washington/edu)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2011 at 23:55 GMT

% todo:

% check nargin
error(nargchk(3,4,nargin));

% default type to all
if(nargin==3 || isempty(type)); type='all'; end

% check inputs
if(isstruct(mt) && isscalar(mt)) % global cmt struct
    % convert to 3x3 in ray coordinates (aki & richards)
    if(~all(isfield(mt,{'mrr', 'mtt' 'mpp' 'mrt' 'mrp' 'mtp'})))
        error('seizmo:radpat:badInput',...
            'MT struct must be have fields: mrr mtt mpp mrt mrp mtp !');
    end
    mt=[mt.mtt mt.mtp mt.mrt; mt.mtp mt.mpp mt.mrp; mt.mrt mt.mrp mt.mrr];
    mt(2:2:8)=-mt(2:2:8);
elseif(isequal(size(mt),[1 6])) % (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
    % convert to 3x3 in ray coordinates (aki & richards)
    mt=mt([2 6 4; 6 3 5; 4 5 1]);
    mt(2:2:8)=-mt(2:2:8);
elseif(isequal(size(mt),[3 3])) % [Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp]
    % convert to ray coordinates (aki & richards)
    mt=mt([2 3 1],[2 3 1]);
    mt(2:2:8)=-mt(2:2:8);
else
    error('seizmo:radpat:badInput',...
        'MT is not in a valid format!');
end
if(~isreal(theta) || ~isvector(theta) || any(theta<0 | theta>180))
    error('seizmo:radpat:badInput',...
        ['THETA must be a vector of real-valued angles ' ...
        'from 0 (down) to 180 (up)!']);
end
if(~isreal(phi) || ~isvector(phi))
    % 0 = North, 90 = East
    error('seizmo:radpat:badInput',...
        ['PHI must be a vector of real-valued azimuth angles ' ...
        '(0 = North, 90 = East)!']);
elseif(~isequalnumelorscalar(theta,phi))
    error('seizmo:radpat:badInput',...
        'THETA & PHI must have equal number of elements or be scalar!');
end
if(~ischar(type) || size(type,1)~=1 ...
        || ~any(strcmpi(type,{'p' 'sh' 'sv' 'all'})))
    error('seizmo:radpat:badInput',...
        'TYPE must be ''p'' ''sh'' ''sv'' or ''all''!');
end

% expand scalar theta/phi
[theta,phi]=expandscalars(theta,phi);

% force to be column vectors
theta=theta(:); phi=phi(:);

% normalize mt
mt=mt/sqrt(trace(mt*mt)/2);

% compute p/sv/sh directions
gamma=[sind(theta).*cosd(phi) sind(theta).*sind(phi) cosd(theta)];
sv=[cosd(theta).*cosd(phi) cosd(theta).*sind(phi) -sind(theta)];
sh=[-sind(phi) cosd(phi) zeros(size(phi))];

% get displacements
u=nan(numel(theta),1);
switch lower(type)
    case 'p'
        for i=1:numel(theta); u(i,1)=gamma(i,:)*mt*gamma(i,:).'; end
    case 'sv'
        for i=1:numel(theta); u(i,1)=sv(i,:)*mt*gamma(i,:).'; end
    case 'sh'
        for i=1:numel(theta); u(i,1)=sh(i,:)*mt*gamma(i,:).'; end
    case 'all'
        u=[u u u];
        for i=1:numel(theta)
            u(i,1)=gamma(i,:)*mt*gamma(i,:).';
            u(i,2)=sv(i,:)*mt*gamma(i,:).';
            u(i,3)=sh(i,:)*mt*gamma(i,:).';
        end
end

end

