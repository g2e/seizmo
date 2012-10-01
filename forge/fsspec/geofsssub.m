function [s]=geofsssub(s,frng,srng,latrng,lonrng,fidx,sidx)
%GEOFSSSUB    Extracts a subspectra of a geofss struct
%
%    Usage:    s=geofsssub(s,frng,srng,latrng,lonrng)
%              s=geofsssub(s,frng,srng,latrng,lonrng,fidx,sidx)
%
%    Description:
%     S=GEOFSSSUB(S,FRNG,SRNG,LATRNG,LONRNG) reduces the geofss spectra in
%     S to the ranges given by FRNG, SRNG, LATRNG, & LONRNG and returns S
%     with all pertinent info updated.  S is a struct (see GEOFSS for
%     details).  FRNG specifies the frequency range as [FREQLOW FREQHIGH]
%     in Hz.  SRNG is the horizontal slowness range as [SLOWLOW SLOWHIGH]
%     in sec/deg.  LATRNG & LONRNG provide the parallels and meridians that
%     define a window in the geographic domain and they are expected to be
%     in degrees as [LOW HIGH].  Note that longitudes positions will be
%     wrapped so all possible points are included in the LONRNG.
%
%     S=GEOFSSSUB(S,FRNG,SRNG,LATRNG,LONRNG,FIDX,SIDX) also reduces the
%     geofss spectra in S to the frequencies & horizontal slownesses with
%     indices corresponding to FIDX and SIDX.  Note that only those with
%     within the range defined by FRNG & SRNG are returned.
%
%     By default all ranges are set to [] (full range).
%
%    Notes:
%     - Although this will wrap the longitudes into the longitude range,
%       they are NOT rearranged to correspond to a monotonic grid as
%       required by PLOTGEOFSS (ie the dateline will be an issue if your
%       grid did not extend across it in the first place).
%
%    Examples:
%     % Split a spectra into 10s period windows:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     s=geofss(d,[lat(:) lon(:)],27:33,[1/50 1/20],'center');
%     s1=geofsssub(s,[1/50 1/40]);
%     s2=geofsssub(s,[1/40 1/30]);
%     s3=geofsssub(s,[1/30 1/20]);
%
%    See also: GEOFSS, GEOFSSXC, GEOFSSAVG, PLOTGEOFSS, GEOFSSFREQSLIDE,
%              GEOFSSSLOWSLIDE, GEOFSSFRAMESLIDE

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Apr.  4, 2012 - minor doc update
%        June  8, 2012 - adapted from geofksubvol
%        Sep. 29, 2012 - handle no vector field and slow function handle,
%                        added fidx & sidx options
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,7,nargin));

% check struct
error(chkgeofss(s));

% don't allow averaged spectrums
% - allowed spectra averaged across only one domain previously
for a=1:numel(s)
    if(size(s(a).spectra,2)~=numel(s(a).slow) ...
            || size(s(a).spectra,3)~=numel(s(a).freq))
        error('seizmo:geofsssub:badInput',...
            'S(%d) is already averaged: cannot extract a subspectrum!',a);
    end
end

% check frequency range
if(nargin<2); frng=[]; end
if(nargin<3); srng=[]; end
if(nargin<4); latrng=[]; end
if(nargin<5); lonrng=[]; end
if(nargin<6); fidx=[]; end
if(nargin<7); sidx=[]; end
sf=size(frng);
ss=size(srng);
sla=size(latrng);
slo=size(lonrng);
if(~isempty(frng) && (~isreal(frng) || numel(sf)~=2 || sf(1)~=1 ...
        || sf(2)~=2 || any(frng(:)<=0)))
    error('seizmo:geofsssub:badInput',...
        'FRNG must be a 1x2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(~isempty(srng) && (~isreal(srng) || numel(ss)~=2 || ss(1)~=1 ...
        || ss(2)~=2 || any(srng(:)<=0)))
    error('seizmo:geofsssub:badInput',...
        'SRNG must be a 1x2 array of [SLOWLOW SLOWHIGH] in sec/deg!');
elseif(~isempty(latrng) && (~isreal(latrng) || numel(sla)~=2 ...
        || sla(1)~=1 || sla(2)~=2 || any(abs(latrng(:))>90)))
    error('seizmo:geofsssub:badInput',...
        'LATRNG must be a 1x2 array of [LATLOW LATHIGH] in degrees!');
elseif(~isempty(lonrng) && (~isreal(lonrng) || numel(slo)~=2 ...
        || slo(1)~=1 || slo(2)~=2))
    error('seizmo:geofsssub:badInput',...
        'LONRNG must be a 1x2 array of [LONLOW LONHIGH] in degrees!');
elseif(~islogical(fidx) && ~isnumeric(fidx))
    error('seizmo:geofsssub:badInput',...
        'FIDX must be indices!');
elseif(~islogical(sidx) && ~isnumeric(sidx))
    error('seizmo:geofsssub:badInput',...
        'SIDX must be indices!');
end

% extract subspectra
for i=1:numel(s)
    % freq range
    if(isempty(frng))
        fmin=min(s(i).freq);
        fmax=max(s(i).freq);
    else
        fmin=frng(1);
        fmax=frng(2);
    end
    fidx0=s(i).freq>=fmin & s(i).freq<=fmax;
    
    % direct frequency indexing
    if(isempty(fidx))
        % do nothing
    elseif(islogical(fidx))
        if(~isequal(size(fidx0),size(fidx)))
            error('seizmo:geofsssub:badInput',...
                'FIDX logical indexing does not match frequency size!');
        end
        fidx0=fidx0 & fidx;
    else
        if(any(fidx>numel(fidx0)) || any(fidx<0))
            error('seizmo:geofsssub:badInput',...
                'FIDX indices outside indexing range!');
        end
        fidx0=intersect(find(fidx0),fidx);
    end
    
    % slow range
    if(isa(s(i).slow,'function_handle'))
        if(~isempty(srng) || ~isempty(sidx))
            error('seizmo:geofsssub:badInput',...
                'Cannot index into slowness function handle!');
        end
        sidx0=true;
    else
        if(isempty(srng))
            smin=min(s(i).slow);
            smax=max(s(i).slow);
        else
            smin=srng(1);
            smax=srng(2);
        end
        sidx0=s(i).slow>=smin & s(i).slow<=smax;
        
        % direct slowness indexing
        if(isempty(sidx))
            % do nothing
        elseif(islogical(sidx))
            if(~isequal(size(sidx0),size(sidx)))
                error('seizmo:geofsssub:badInput',...
                    'SIDX logical indexing does not match slowness size!');
            end
            sidx0=sidx0 & sidx;
        else
            if(any(sidx>numel(sidx0)) || any(sidx<0))
                error('seizmo:geofsssub:badInput',...
                    'SIDX indices outside indexing range!');
            end
            sidx0=intersect(find(sidx0),sidx);
        end
    end
    
    % lat range
    if(isempty(latrng))
        latmin=min(s(i).latlon(:,1));
        latmax=max(s(i).latlon(:,1));
    else
        latmin=latrng(1);
        latmax=latrng(2);
    end
    latidx=s(i).latlon(:,1)>=latmin & s(i).latlon(:,1)<=latmax;
    
    % lon range
    if(isempty(lonrng))
        lonmin=min(s(i).latlon(:,2));
        lonmax=max(s(i).latlon(:,2));
    else
        lonmin=lonrng(1);
        lonmax=lonrng(2);
    end
    while(any(lonmin-s(i).latlon(:,2)>180))
        moveme=lonmin-s(i).latlon(:,2)>180;
        s(i).latlon(moveme,2)=s(i).latlon(moveme,2)+360;
    end
    while(any(s(i).latlon(:,2)-lonmax>180))
        moveme=s(i).latlon(:,2)-lonmax>180;
        s(i).latlon(moveme,2)=s(i).latlon(moveme,2)-360;
    end
    lonidx=s(i).latlon(:,2)>=lonmin & s(i).latlon(:,2)<=lonmax;
    
    % check for output
    if(isempty(fidx0) || ~any(fidx0))
        error('seizmo:geofsssub:badInput',...
            'FRNG+FIDX selects nothing in S(%d) frequency range!',i);
    elseif(isempty(sidx0) || ~any(sidx0))
        error('seizmo:geofsssub:badInput',...
            'SRNG+SIDX selects nothing in S(%d) slowness range!',i);
    elseif(~any(latidx) || ~any(lonidx))
        error('seizmo:geofsssub:badInput',...
            'LATRNG/LONRNG outside S(%d) lat/lon range!',i);
    end
    
    % extract subspectra
    s(i).spectra=s(i).spectra(latidx & lonidx,sidx0,fidx0);
    
    % update info
    s(i).freq=s(i).freq(fidx0);
    s(i).slow=s(i).slow(sidx0);
    s(i).latlon=s(i).latlon(latidx & lonidx,:);
end

end
