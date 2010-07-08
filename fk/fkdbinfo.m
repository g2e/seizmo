function [mindb,meddb,maxdb]=fkdbinfo(fk,frng,bazrng,srng,esrng,nsrng)
%FKDBINFO    Returns the min/median/max dB for a FK struct
%
%    Usage:    [mindb,meddb,maxdb]=fkdbinfo(fk)
%             [mindb,meddb,maxdb]=fkdbinfo(fk,frng,bazrng,srng,esrng,nsrng)
%
%    Description: [MINDB,MEDDB,MAXDB]=FKDBINFO(FK) returns the limits and
%     median of the beam(s) in FK in decibels.  This is useful for
%     quick identification of elements with strong plane wave coherency.
%
%     [MINDB,MEDDB,MAXDB]=FKDBINFO(FK,FRNG,BAZRNG,SRNG,ESRNG,NSRNG) returns
%     dB info about a particular subsection of the slowness volume(s) or
%     map(s).  FRNG is the frequency range in Hz.  BAZRNG is the back-
%     azimuth range in degrees.  SRNG is the horizontal slowness magnitude
%     range in sec/deg.  ESRNG & NSRNG are the East & North slowness ranges
%     in sec/deg.  The defaults for all ranges will include all data.
%
%    Notes:
%
%    Examples:
%     Analyze a 4D fk dataset:
%      s4d=fk4d(data,[],[],50,201,[1/100 1/4]);
%      [mindb,meddb,maxdb]=fkdbinfo(s4d);
%      figure; plot(mindb);
%      hold on;
%      plot(maxdb);
%      plot(meddb);
%
%    See also: FKMAP, FKVOLUME, FK4D

%     Version History:
%        May  26, 2010 - initial version
%        June 30, 2010 - added range arguments
%        July  1, 2010 - range bugfix
%        July  6, 2010 - update for new struct
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July  6, 2010 at 16:45 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check fk struct
error(chkfkstruct(fk));

% check ranges
if(nargin<2); frng=[]; end
if(nargin<3); bazrng=[]; end
if(nargin<4); srng=[]; end
if(nargin<5); esrng=[]; end
if(nargin<6); nsrng=[]; end
sf=size(frng);
sbaz=size(bazrng);
ss=size(srng);
ses=size(esrng);
sns=size(nsrng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:fkdbinfo:badInput',...
        'FRNG must be a 1x2 array of [FREQLO FREQHI] in Hz!');
end
if(~isempty(bazrng) && (~isreal(bazrng) || numel(sbaz)~=2 || sbaz(1)~=1 ...
        || sbaz(2)~=2))
    error('seizmo:fkdbinfo:badInput',...
        'BAZRNG must be a 1x2 array of [BAZLO BAZHI] in deg!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:fkdbinfo:badInput',...
        'SRNG must be a 1x2 array of [SLO SHI] in sec/deg!');
end
if(~isempty(esrng) && (~isreal(esrng) || numel(ses)~=2 || ses(1)~=1 ...
        || ses(2)~=2))
    error('seizmo:fkdbinfo:badInput',...
        'ESRNG must be a 1x2 array of [EAST_SLO EAST_SHI] in sec/deg!');
end
if(~isempty(nsrng) && (~isreal(nsrng) || numel(sns)~=2 || sns(1)~=1 ...
        || sns(2)~=2))
    error('seizmo:fkdbinfo:badInput',...
        'NSRNG must be a 1x2 array of [NORTH_SLO NORTH_SHI] in sec/deg!');
end

% loop over every element
mindb=nan(size(fk));
meddb=mindb;
maxdb=mindb;
for i=1:numel(fk)
    % get backazimuths and slownesses
    if(fk(i).polar)
        nx=numel(fk(i).x);
        ny=numel(fk(i).y);
        baz=fk(i).x(ones(ny,1),:);
        s=fk(i).y(:,ones(nx,1));
        [es,ns]=slowbaz2kxy(s,baz);
    else % cartesian
        nx=numel(fk(i).x);
        ny=numel(fk(i).y);
        es=fk(i).x(ones(ny,1),:);
        ns=fk(i).y(:,ones(nx,1));
        [s,baz]=kxy2slowbaz(es,ns);
    end
    
    % get extraction indices
    if(fk(i).volume)
        if(isempty(frng))
            fidx=true(1,numel(fk(i).freq));
        else
            fidx=fk(i).freq>=frng(1) & fk(i).freq<=frng(2);
        end
    else
        if(~isempty(frng))
            warning('seizmo:fkdbinfo:badInput',...
                'FRNG does not work for maps!');
        end
        fidx=true;
    end
    if(isempty(bazrng))
        bazidx=true(size(baz));
    else
        bazidx=baz>=bazrng(1) & baz<=bazrng(2);
    end
    if(isempty(srng))
        sidx=true(size(s));
    else
        sidx=s>=srng(1) & s<=srng(2);
    end
    if(isempty(esrng))
        esidx=true(size(es));
    else
        esidx=es>=esrng(1) & es<=esrng(2);
    end
    if(isempty(nsrng))
        nsidx=true(size(ns));
    else
        nsidx=ns>=nsrng(1) & ns<=nsrng(2);
    end
    
    % subsection
    idx=repmat(bazidx & sidx & esidx & nsidx,[1 1 numel(fidx)]);
    idx(:,:,~fidx)=false;
    sub=fk(i).beam(idx);
    
    % stats
    mindb(i)=min(sub(:))+fk(i).normdb;
    meddb(i)=median(sub(:))+fk(i).normdb;
    maxdb(i)=max(sub(:))+fk(i).normdb;
end

end
