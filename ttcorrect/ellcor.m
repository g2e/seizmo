function [varargout]=ellcor(evla,evdp,gcarc,az,varargin)
%ELLCOR    Returns ellipticity corrections
%
%    Usage:    corr=ellcor(evla,evdp,gcarc,az,phase)
%              [corr1,...,corrN]=ellcor(evla,evdp,gcarc,az,ph1,...,phN)
%
%    Description:
%     CORR=ELLCOR(EVLA,EVDP,GCARC,AZ,PHASE) calculates the effect of an
%     elliptical Earth (WGS-84 ellipsoid & AK135 velocity structure) on the
%     travel times of seismic phases.  The effect is given so that
%     TTsphere+CORR=TTellipsoid where TTsphere is the travel time of the
%     phase for the AK135 1D Earth model.  PHASE may be any of the
%     following phases (given as a string):
%      p,P,Pdiff,PcP,PP,PS,PcS,SP,ScP,s,S,Sdiff,ScS,SS.
%     EVLA, EVDP, GCARC, & AZ should all be scalar or equal sized.  EVLA
%     is the geocentric latitude of the earthquake(s) and should be in
%     degrees.  EVDP is the kilometer depth of the earthquake(s).  GCARC is
%     the distance from the earthquake to the station in degrees and AZ is
%     the earthquake to station azimuth in degrees.  Note that if the
%     earthquake-station pair is outside the phase tabulation the
%     correction will be NaN!
%
%     [CORR1,...,CORRN]=ELLCOR(EVLA,EVDP,GCARC,AZ,PHASE1,...,PHASEN)
%     retrieves ellipticity corrections for multiple phases.  CORR1..N
%     correspond to the corrections for PHASE1..N.
%
%    Notes:
%     - Phases ellipticity tables are in AK135_ELLIP_TABLES, so go there to
%       add more phases to the tabulation.
%
%    Examples:
%     % Plot up an ellipticity correction cloud for P-wave arrivals (note
%     % we use RANDLATLON to get a spacially uniform random distribution of
%     % events and stations):
%     n=100000;
%     [evla,evlo]=randlatlon(n);
%     [stla,stlo]=randlatlon(n);
%     [gcarc,az]=sphericalinv(evla,evlo,stla,stlo);
%     evdp=700*rand(n,1);
%     [m,c]=hist3([evla ellcor(evla,evdp,gcarc,az,'P')],[100 100]);
%     figure; imagesc(c{1},fliplr(c{2}),rot90(m)); axis xy;
%     [m,c]=hist3([evdp ellcor(evla,evdp,gcarc,az,'P')],[100 100]);
%     figure; imagesc(c{1},fliplr(c{2}),rot90(m)); axis xy;
%     [m,c]=hist3([gcarc ellcor(evla,evdp,gcarc,az,'P')],[100 100]);
%     figure; imagesc(c{1},fliplr(c{2}),rot90(m)); axis xy;
%     [m,c]=hist3([az ellcor(evla,evdp,gcarc,az,'P')],[100 100]);
%     figure; imagesc(c{1},fliplr(c{2}),rot90(m)); axis xy;
%
%     % Same, but for Pdiff:
%     n=100000;
%     [evla,evlo]=randlatlon(n);
%     [stla,stlo]=randlatlon(n);
%     [gcarc,az]=sphericalinv(evla,evlo,stla,stlo);
%     evdp=700*rand(n,1);
%     [m,c]=hist3([evla ellcor(evla,evdp,gcarc,az,'Pdiff')],[100 100]);
%     figure; imagesc(c{1},fliplr(c{2}),rot90(m)); axis xy;
%     [m,c]=hist3([evdp ellcor(evla,evdp,gcarc,az,'Pdiff')],[100 100]);
%     figure; imagesc(c{1},fliplr(c{2}),rot90(m)); axis xy;
%     [m,c]=hist3([gcarc ellcor(evla,evdp,gcarc,az,'Pdiff')],[100 100]);
%     figure; imagesc(c{1},fliplr(c{2}),rot90(m)); axis xy;
%     [m,c]=hist3([az ellcor(evla,evdp,gcarc,az,'Pdiff')],[100 100]);
%     figure; imagesc(c{1},fliplr(c{2}),rot90(m)); axis xy;
%
%    See also: AK135_ELLIP_TABLES, CRUCOR, MANCOR

%     Version History:
%        May  16, 2010 - initial version
%        Jan.  7, 2011 - update example to use RANDLATLON
%        Apr.  2, 2012 - minor doc update
%        Aug.  6, 2012 - nargin check added, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  6, 2012 at 12:25 GMT

% todo:

% check number of inputs
error(nargchk(5,inf,nargin));

% check location inputs
if(~isreal(evla) || ~isreal(evdp) || ~isreal(gcarc) || ~isreal(az))
    error('seizmo:ellcor:badInput',...
        'EVLA, EVDP, GCARC & AZ must be real-valued arrays!');
elseif(~isequalsizeorscalar(evla,evdp,gcarc,az))
    error('seizmo:get_scripps_value:badInput',...
        'EVLA, EVDP, GCARC & AZ must be equal sized or scalar!');
end

% expand scalar location info
[evla,evdp,gcarc,az]=expandscalars(evla,evdp,gcarc,az);

% check phases
phases={'p' 'P' 'Pdiff' 'PcP' 'PP' 'PS' 'PcS' ...
    'SP' 'ScP' 's' 'S' 'Sdiff' 'ScS' 'SS'};
if(any(~ismember(varargin,phases)))
    tf=~ismember(varargin,phases);
    error('seizmo:ellcor:badPhase',...
        ['Unknown Phase(s): ' sprintf('%s ',varargin{tf})]);
end

% load ellipticity tables
table=ak135_ellip_tables;

% event lat to colat & deg to rad (for azimuth too)
d2r=pi/180;
evcola=(90-evla)*d2r;
az=az*d2r;

% for efficiency precompute coefficients based on event location
threetwo=sqrt(3)/2;
corr0=1/4*(1+3*cos(2*evcola));
corr1=threetwo*sin(2*evcola).*cos(az);
corr2=threetwo*sin(evcola).^2.*cos(2*az);

% loop over phases
nphases=nargin-4;
for i=1:nphases
    % get tau values for corrections
    tau0=interp2(table.(varargin{i}).depth,table.(varargin{i}).distance,...
        table.(varargin{i}).t0,evdp,gcarc);
    tau1=interp2(table.(varargin{i}).depth,table.(varargin{i}).distance,...
        table.(varargin{i}).t1,evdp,gcarc);
    tau2=interp2(table.(varargin{i}).depth,table.(varargin{i}).distance,...
        table.(varargin{i}).t2,evdp,gcarc);
    
    % debugging
    % [tau0 tau1 tau2 corr0 corr1 corr2]
    
    % get corrections
    varargout{i}=corr0.*tau0+corr1.*tau1+corr2.*tau2;
end

end
