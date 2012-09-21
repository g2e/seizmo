function [idx]=fssidx(s,frng,bazrng,srng,esrng,nsrng)
%FSSIDX    Returns logical index matrix for fss spectra
%
%    Usage:    idx=fssidx(s,frng,bazrng,srng,esrng,nsrng)
%
%    Description:
%     IDX=FSSIDX(S,FRNG,BAZRNG,SRNG,ESRNG,NSRNG) returns the logical
%     indexing matrix IDX which is the same size as the frequency-slowness
%     spectra in fss struct S for the portion of the spectra that is within
%     the specified ranges.  S must be a struct as output from function
%     FSS.  FRNG is the frequency range as [FREQLO FREQHI] in Hz.  BAZRNG
%     is the back-azimuth range as [BAZLO BAZHI] in degrees.  SRNG is the
%     horizontal slowness magnitude range as [SLO SHI] in sec/deg.  ESRNG
%     is the east horizontal slowness range as [SLO SHI] in sec/deg.  NSRNG
%     is the north horizontal slowness range as [SLO SHI] in sec/deg.
%
%    Notes:
%     - If S has multiple elements IDX is a cell array of equal size to S
%       where each cell gives the corresponding logical matrix.
%
%    Examples:
%     % Set spectra outside of 0 to 90 degrees to 0:
%     idx=fssidx(s,[],[0 90]);
%     plotfss(fssavg(s))
%     s.spectra(~idx)=0;
%     plotfss(fssavg(s))
%
%    See also: FSSGRIDS, FSSDBINFO, FSSSUB, FSSAVG,
%              FSS, FSSXC, FSSHORZ, FSSHORZXC

%     Version History:
%        Sep. 12, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 12, 2012 at 14:25 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check fk struct
error(chkfss(s));

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
    error('seizmo:fssidx:badInput',...
        'FRNG must be a 1x2 array of [FREQLO FREQHI] in Hz!');
end
if(~isempty(bazrng) && (~isreal(bazrng) || numel(sbaz)~=2 || sbaz(1)~=1 ...
        || sbaz(2)~=2))
    error('seizmo:fssidx:badInput',...
        'BAZRNG must be a 1x2 array of [BAZLO BAZHI] in deg!');
end
if(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:fssidx:badInput',...
        'SRNG must be a 1x2 array of [SLO SHI] in sec/deg!');
end
if(~isempty(esrng) && (~isreal(esrng) || numel(ses)~=2 || ses(1)~=1 ...
        || ses(2)~=2))
    error('seizmo:fssidx:badInput',...
        'ESRNG must be a 1x2 array of [EAST_SLO EAST_SHI] in sec/deg!');
end
if(~isempty(nsrng) && (~isreal(nsrng) || numel(sns)~=2 || sns(1)~=1 ...
        || sns(2)~=2))
    error('seizmo:fssidx:badInput',...
        'NSRNG must be a 1x2 array of [NORTH_SLO NORTH_SHI] in sec/deg!');
end

% loop over every spectra
idx=cell(size(s));
for i=1:numel(s)
    % get grids of backazimuth and slowness
    if(s(i).polar)
        nx=numel(s(i).x);
        ny=numel(s(i).y);
        baz=s(i).x(ones(ny,1),:);
        smag=s(i).y(:,ones(nx,1));
        [es,ns]=slowbaz2kxy(smag,baz);
    else % cartesian
        nx=numel(s(i).x);
        ny=numel(s(i).y);
        es=s(i).x(ones(ny,1),:);
        ns=s(i).y(:,ones(nx,1));
        [smag,baz]=kxy2slowbaz(es,ns);
    end
    
    % get individual logical indices
    if(numel(s(i).freq)==size(s(i).spectra,3))
        if(isempty(frng))
            fidx=true(1,numel(s(i).freq));
        else
            fidx=s(i).freq>=frng(1) & s(i).freq<=frng(2);
        end
    else
        if(~isempty(frng))
            warning('seizmo:fssidx:badInput',...
                'FRNG does not work for averaged spectra!');
        end
        fidx=true;
    end
    if(isempty(bazrng))
        bazidx=true(size(baz));
    else
        bazidx=(baz>=bazrng(1) & baz<=bazrng(2)) ...
            | (baz>=(bazrng(1)-360) & baz<=(bazrng(2)-360)) ...
            | (baz>=(bazrng(1)+360) & baz<=(bazrng(2)+360));
    end
    if(isempty(srng))
        sidx=true(size(smag));
    else
        sidx=smag>=srng(1) & smag<=srng(2);
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
    
    % combine
    idx{i}=repmat(bazidx & sidx & esidx & nsidx,[1 1 numel(fidx)]);
    idx{i}(:,:,~fidx)=false;
end

% uncell scalar spectra
if(isscalar(s)); idx=idx{1}; end

end
