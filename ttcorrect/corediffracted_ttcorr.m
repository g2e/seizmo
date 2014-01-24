function [ecor,ccor,mcorf,mcoru]=corediffracted_ttcorr(data,phase,mod)
%COREDIFFRACTED_TTCORR    Traveltime corrections for core-diffracted waves
%
%    Usage:    [ecor,ccor,mcorfull,mcorup]=corediffracted_ttcorr(data,...
%                                                              phase,model)
%
%    Description:
%     [ECOR,CCOR,MCORFULL,MCORUP]=COREDIFFRACTED_TTCORR(DATA,PHASE,MODEL)
%     gets travel time corrections for core-diffracted arrivals expected
%     at stations in SEIZMO struct DATA.  PHASE must either be 'Pdiff' or
%     'Sdiff'.  MODEL is optional and may be any listed by
%     AVAILABLE_1DMODELS.  ECOR is the ellipsoid travel time correction,
%     CCOR is the crustal travel time correction, MCORFULL is the full
%     raypath mantle travel time correction and MCORUP is the upswing
%     raypath mantle travel time correction.  Corrections follow this rule:
%                         TT3D=TT1D+CORRECTIONS
%     All corrections are in seconds.  If you do not choose a model, PREM
%     is the default.
%
%    Notes:
%     - Ellipticity corrections use AK135, not the selected model.
%
%    Examples:
%     % Get corrections for a core-diffracted dataset and apply them:
%     [ecor,ccor,mcorfull,mcorup]=corediffracted_ttcorr(data,'Pdiff');
%     data=makearrivals(data,'P,Pdiff','prem');
%     Pdiff=findpicks(arrivals2picks(data),'P,Pdiff',1);
%     cddata=timeshift(data,-(Pdiff+ecor+ccor+mcorfull));
%     plot0(cddata);
%
%    See also: ELLCOR, CRUCOR, MANCOR, GETRAYPATHS, PLOTRAYPATHS,
%              EXTRACT_UPSWING_RAYPATHS, CRUSTLESS_RAYPATHS, TAUPPATH

%     Version History:
%        June 11, 2010 - initial version
%        Aug.  8, 2010 - cleaned up docs & code
%        Feb. 27, 2012 - doc update
%        Mar. 15, 2012 - fix example
%        Jan. 23, 2014 - minor doc update, bugfix: call CRUSTLESS_RAYPATHS
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 14:25 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% check struct & header
data=checkheader(data,...
    'UNSET_ST_LATLON','ERROR',...
    'UNSET_EV_LATLON','ERROR',...
    'UNSET_ELEV','FIX',...
    'UNSET_DEPTH','FIX');

% check phase
if(~ischar(phase))
    error('seizmo:corediffracted_ttcorr:badPHASE',...
        'PHASE must be either ''Pdiff'' or ''Sdiff''!');
end

% default model (prem)
if(nargin==2 || isempty(mod)); mod='prem'; end
if(~ischar(mod) || ~ismember(mod,available_1dmodels))
    error('seizmo:corediffracted_ttcorr:badMOD',...
        'MODEL must be one listed by AVAILABLE_1DMODELS!');
end

% cmb depth
switch lower(mod)
    case 'iasp91'
        cmb=2889;
    case 'ak135'
        cmb=2891.5;
    case 'prem'
        cmb=2891;
end

% get info from header
[evla,evlo,evdp,stla,stlo,stel,stdp,gcarc,az]=getheader(data,...
    'evla','evlo','evdp','stla','stlo','stel','stdp','gcarc','az');

% convert meters to kilometers
evdp=evdp/1000;
stel=stel/1000;
stdp=stdp/1000;

switch lower(phase)
    case 'pdiff'
        % get raypaths and rayparameters
        paths=getraypaths('P,Pdiff',mod,evla,evlo,evdp,stla,stlo);
        paths=crustless_raypaths(paths);
        uppaths=extract_upswing_raypaths(paths,cmb-350);
        rp=[paths.rayparameter]';
        
        % get ellipsoid corrections
        ecor=ellcor(evla,evdp,gcarc,az,'Pdiff');
        
        % get crustal corrections
        ccor=crucor(stla,stlo,rp,'p','elev',stel,'hole',stdp,'refmod',mod);
        
        % get mantle corrections
        mcorf=mancor(paths,'hmsl06p');
        mcoru=mancor(uppaths,'hmsl06p');
    case {'sdiff' 'svdiff' 'shdiff'}
        % get raypaths and rayparameters
        paths=getraypaths('S,Sdiff',mod,evla,evlo,evdp,stla,stlo);
        paths=crustless_raypaths(paths);
        uppaths=extract_upswing_raypaths(paths,cmb-350);
        rp=[paths.rayparameter]';
        
        % get ellipsoidal correction
        ecor=ellcor(evla,evdp,gcarc,az,'Sdiff');
        
        % get crustal corrections
        ccor=crucor(stla,stlo,rp,'s','elev',stel,'hole',stdp,'refmod',mod);
        
        % get mantle corrections
        mcorf=mancor(paths,'hmsl06s');
        mcoru=mancor(uppaths,'hmsl06s');
    otherwise
        error('seizmo:corediffracted_ttcorr:badPhase',...
            'Unknown Phase: %s',phase);
end

end
