function [arf]=geofkarf2map(arf,srng)
%GEOFKARF2MAP    Converts a geofk ARF volume to a geofk ARF map
%
%    Usage:    arf=geofkarf2map(arf,srng)
%
%    Description:
%     ARF=GEOFKARF2MAP(ARF,SRNG) averages the geofk array beam response
%     data in ARF across the horizontal slowness range SRNG, returning the
%     result in ARF.  ARF may then be plotted using PLOTGEOFKARF.  ARF is a
%     GEOFKARF struct (see that function for details on the struct layout).
%     SRNG gives the horizontal slowness range to average over as
%     [SLOWLOW SLOWHIGH] in sec/deg.
%
%    Notes:
%
%    Examples:
%     % Normal usage:
%     plotgeofkarf(geofkarf2map(geofkarf(...)))
%
%    See also: GEOFKARF, PLOTGEOFKARF, GEOFKSUBARF, UPDATEGEOFKARF,
%              GEOFKARFSLOWSLIDE, CHKGEOFKARFSTRUCT

%     Version History:
%        July  7, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check fk struct
error(chkgeofkarfstruct(arf));

% don't allow map
if(any(~sum(reshape([arf.volume],[2 numel(arf)]))))
    error('seizmo:geofkarf2map:badInput',...
        'ARF must be a GEOFKARF volume!');
end

% check frequency range
if(nargin<2); srng=[]; end
ss=size(srng);
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofkarf2map:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
end

% average frequencies/slownesses
for i=1:numel(arf)
    % handle empty srng
    if(isempty(srng))
        smin=min(arf(i).horzslow);
        smax=max(arf(i).horzslow);
    else
        if(~arf(i).volume(1))
            error('seizmo:geofkvol2map:badVol',...
                'Can not average a slowness range from averaged volume!');
        end
        smin=srng(1);
        smax=srng(2);
    end
    
    % average (if possible)
    if(arf(i).volume(1))
        sidx=arf(i).horzslow>=smin & arf(i).horzslow<=smax;
        arf(i).beam=sum(arf(i).beam(:,sidx),2)/sum(sidx);
        arf(i).horzslow=arf(i).horzslow(sidx);
    end
    arf(i).volume=[false false];
    
    % rescale
    maxdb=max(arf(i).beam(:));
    arf(i).beam=arf(i).beam-maxdb;
    arf(i).normdb=arf(i).normdb+maxdb;
end

end
