function [ns,es]=fssphasecorr(phase,mod3d,evla,evlo,evdp,stla,stlo)
%FSSPHASECORR    Calculate array slowness corrections for a seismic phase
%
%    Usage:    [ns,es]=fssphasecorr(phase,mod3d,evla,evlo,evdp,stla,stlo)
%
%    Description:
%     [NS,ES]=FSSPHASECORR(PHASE,MOD3D,EVLA,EVLO,EVDP,STLA,STLO) returns
%     the north and east slowness deviations for a seismic phase PHASE from
%     a source EVLA/EVLO/EVDP passing through 3D mantle model MODEL and
%     recorded by a seismic array defined by STLA/STLO.  Additional
%     corrections for crustal heterogeneity in Crust1.0 and Earth
%     ellipticity are included.  The slownesses are the apparent moveout of
%     the phase from heterogeneity and ellipticity, which can be applied to
%     the moveout for the phase through a 1D Earth model to better match
%     the observed moveout for a phase across an array.
%
%    Notes:
%     - Latitudes are assumed to be geographic.
%     - EVDP is expected in kilometers!
%     - Phase handling is basic - does not handle exotic phases!
%
%    Examples:
%     % Get the expected slowness of a P wave passing through the
%     % Yellowknife array from a source near Japan, then get the effect of
%     % 3D mantle heterogeneity on travel times and thus slowness across
%     % the array.  Finally calculate the apparent source location:
%     evla=40; evlo=140; evdp=0;
%     st=yellowknife;
%     [clat,clon]=arraycenter(st(:,1),st(:,2));
%     evlas=geographic2geocentriclat(evla);
%     clat=geographic2geocentriclat(clat);
%     [gcarc,az,baz]=sphericalinv(evlas,evlo,clat,clon);
%     slow=deg2slowness('P',gcarc);
%     [es,ns]=slowbaz2kxy(slow,baz);
%     [dns,des]=fssphasecorr('P','HMSL06P',evla,evlo,evdp,st(:,1),st(:,2));
%     [newslow,baz]=kxy2slowbaz(es+des,ns+dns);
%     deg=slowness2deg('P',newslow);
%     [evla2,evlo2]=sphericalfwd(clat,clon,deg,baz);
%     [geocentric2geographiclat(evla2) evlo2]
%
%    See also: FSS, FSSXC, SLOWNESS2DEG, DEG2SLOWNESS, GETRAYPATHS,
%              TAUPPATH, MANCOR, MANTLEDV, SELECT_FK_PEAKS

%     Version History:
%        Aug.  6, 2012 - initial version
%        Aug. 30, 2012 - fixed forced phase/wavetype inputs
%        Jan. 31, 2014 - update for CRUST2.0 to CRUST1.0 switch, switched
%                        bad phase removal to match that of TauP 2.1.1
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 31, 2014 at 23:18 GMT

% todo:

% check number of inputs
error(nargchk(7,7,nargin));

% add core-diff for P or S
if(strcmp(phase,'P'))
    phase='P,Pdiff';
    ephase='P';
    wavetype='p';
elseif(strcmp(phase,'S'))
    phase='S,Sdiff';
    ephase='S';
    wavetype='s';
else
    ephase=phase;
    wavetype=lower(phase(1));
end

% get paths through prem
[paths,idx]=getraypaths(phase,'prem',evla,evlo,evdp,stla,stlo);

% deal with missing paths or triplications
%[idx1,idx2]=unique(idx,'last'); % last b/c earlier paths are crap
[idx1,idx2]=unique(idx,'first'); % first b/c later paths are crap
stla=stla(idx1); stlo=stlo(idx1);
paths=paths(idx2);

% get mantle paths
cmb=2891; % prem based cmb depth
gpaths=crustless_raypaths(trim_depths_raypaths(paths,[0 cmb]));

% mantle corrections
mc=mancor(gpaths,mod3d);

% crustal corrections
rayp=[gpaths.rayparameter]; rayp=rayp(:);
cc=crucor(stla,stlo,rayp,wavetype,'s',false);

% ellipticity corrections
[gcarc,az]=sphericalinv(geographic2geocentriclat(evla),evlo,...
    geographic2geocentriclat(stla),stlo);
ec=ellcor(geographic2geocentriclat(evla),evdp,gcarc,az,ephase);

% total corrections
tc=mc+cc+ec;

% east & north positioning
[clat,clon]=arraycenter(stla,stlo);
[e,n]=geographic2enu(stla,stlo,0,clat,clon,0);

% convert km to deg
deg2km=6371*pi/180;
e=e/deg2km;
n=n/deg2km;

% linear fits to get moveout
me=wlinem(e,tc); es=me(2);
mn=wlinem(n,tc); ns=mn(2);

% debugging
%d=[n e mc cc ec];

end
