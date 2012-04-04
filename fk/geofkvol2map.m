function [vol]=geofkvol2map(vol,frng,srng)
%GEOFKVOL2MAP    Converts a geofk volume to a geofk map
%
%    Usage:    map=geofkvol2map(vol,frng,srng)
%
%    Description:
%     MAP=GEOFKVOL2MAP(VOL,FRNG,SRNG) averages the geofk beam data in VOL
%     across the frequency range given by FRNG and the horizontal slowness
%     range SRNG, returning the result in MAP.  MAP may be then plotted
%     using PLOTGEOFKMAP.  VOL is a geofk struct produced by a geofk
%     function (like GEOFKXCVOLUME).  FRNG gives the frequency range to
%     average over as [FREQLOW FREQHIGH] in Hz.  SRNG gives the horizontal
%     slowness range to average over as [SLOWLOW SLOWHIGH] in sec/deg.
%
%    Notes:
%
%    Examples:
%
%    See also: GEOFKFREQSLIDE, GEOFKSLOWSLIDE, PLOTGEOFKMAP, GEOFKXCVOLUME,
%              GEOFKXCHORZVOLUME, CHKGEOFKSTRUCT, GEOFKSUBVOL

%     Version History:
%        June 24, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check fk struct
error(chkgeofkstruct(vol));

% don't allow map
if(any(~sum(reshape([vol.volume],[2 numel(vol)]))))
    error('seizmo:geofkvol2map:badInput',...
        'VOL must be a geofk beam volume!');
end

% check frequency range
if(nargin<2); frng=[]; end
if(nargin<3); srng=[]; end
sf=size(frng);
ss=size(srng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:geofkvol2map:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofkvol2map:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
end

% average frequencies/slownesses
for i=1:numel(vol)
    % handle empty frng
    if(isempty(frng))
        fmin=min(vol(i).freq);
        fmax=max(vol(i).freq);
    else
        if(~vol(i).volume(2))
            error('seizmo:geofkvol2map:badVol',...
                'Can not average a freq range from averaged volume!');
        end
        fmin=frng(1);
        fmax=frng(2);
    end
    
    % handle empty srng
    if(isempty(srng))
        smin=min(vol(i).horzslow);
        smax=max(vol(i).horzslow);
    else
        if(~vol(i).volume(1))
            error('seizmo:geofkvol2map:badVol',...
                'Can not average a slowness range from averaged volume!');
        end
        smin=srng(1);
        smax=srng(2);
    end
    
    % average (if possible)
    if(vol(i).volume(2))
        fidx=vol(i).freq>=fmin & vol(i).freq<=fmax;
        vol(i).beam=sum(vol(i).beam(:,:,fidx),3)/sum(fidx);
        vol(i).freq=vol(i).freq(fidx);
    end
    if(vol(i).volume(1))
        sidx=vol(i).horzslow>=smin & vol(i).horzslow<=smax;
        vol(i).beam=sum(vol(i).beam(:,sidx),2)/sum(sidx);
        vol(i).horzslow=vol(i).horzslow(sidx);
    end
    vol(i).volume=[false false];
    
    % rescale
    maxdb=max(vol(i).beam(:));
    vol(i).beam=vol(i).beam-maxdb;
    vol(i).normdb=vol(i).normdb+maxdb;
end

end
