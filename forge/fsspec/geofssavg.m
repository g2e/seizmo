function [s]=geofssavg(s,frng,srng)
%GEOFSSAVG    Average geofss output over slowness and frequency
%
%    Usage:    s=geofssavg(s,frng,srng)
%
%    Description:
%     S=GEOFSSAVG(S,FRNG,SRNG) averages the frequency-slowness-position
%     spectra in S over the frequency and slowness ranges given by FRNG and
%     SRNG.  S must be a struct as output by GEOFSS.  FRNG is the frequency
%     range to average as [FREQLOW FREQHIGH] in Hz.  SRNG is the horizontal
%     slowness range to average as [SLOWLOW SLOWHIGH] in sec/deg.
%
%    Notes:
%
%    Examples:
%     % Global migration of Rayleigh waves at 25s periods:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     s=geofss(d,[lat(:) lon(:)],30,1./[25.1 24.9],'center');
%     plotgeofss(geofssavg(s));
%
%    See also: GEOFSS, GEOFSSXC, GEOFSSSUB, PLOTGEOFSS

%     Version History:
%        June 24, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%        June  7, 2012 - altered from geofkvol2map
%        Sep. 29, 2012 - handle no vector field and slow func
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check struct
error(chkgeofss(s));

% check frequency range
if(nargin<2); frng=[]; end
if(nargin<3); srng=[]; end
sf=size(frng);
ss=size(srng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<0)))
    error('seizmo:geofssavg:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofssavg:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
end

% average spectrum over frequencies/slownesses
for i=1:numel(s)
    % handle empty frng
    if(isempty(frng))
        fmin=min(s(i).freq);
        fmax=max(s(i).freq);
    else
        fmin=frng(1);
        fmax=frng(2);
    end
    
    % handle empty srng
    if(isempty(srng) && ~isa(s(i).slow,'function_handle'))
        smin=min(s(i).slow);
        smax=max(s(i).slow);
    else
        smin=srng(1);
        smax=srng(2);
    end
    
    % average (if possible)
    if(numel(s(i).freq)==size(s(i).spectra,3))
        fidx=s(i).freq>=fmin & s(i).freq<=fmax;
        if(~any(fidx))
            error('seizmo:fssavg:badFRNG',...
                'FRNG outside S(%d) frequency range!',i);
        end
        s(i).spectra=nanmean(s(i).spectra(:,:,fidx),3);
        s(i).freq=s(i).freq(fidx);
    end
    if(~isa(s(i).slow,'function_handle') ...
            && numel(s(i).slow)==size(s(i).spectra,2))
        sidx=s(i).slow>=smin & s(i).slow<=smax;
        if(~any(fidx))
            error('seizmo:fssavg:badSRNG',...
                'SRNG outside S(%d) slowness range!',i);
        end
        s(i).spectra=nanmean(s(i).spectra(:,sidx,:),2);
        s(i).slow=s(i).slow(sidx);
    end
end

end
