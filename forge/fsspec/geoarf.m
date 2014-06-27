function [s]=geoarf(stlalo,lalo,slow,lalo0,slow0,f0,varargin)
%GEOARF    Calculate the geographic array response function for an array
%
%    Usage:    s=geoarf(stlalo,latlon,slow,latlon0,slow0,f0)
%              s=geoarf(...,'method',string,...)
%              s=geoarf(...,'weights',w,...)
%              s=geoarf(...,'avg',true|false,...)
%
%    Description:
%     S=GEOARF(STLALO,LATLON,SLOW,lalo0,slow0,f0) computes the array
%     response function (ARF) for an array at the locations in STLALO with
%     waves passing through the array defined by their source positions
%     latlon0, horizontal slownesses slow0, and frequencies f0.  GEOARF
%     shows how well the array is able to resolve the source locations with
%     beamforming methods.  The response is computed at the positions
%     LATLON (to see the aliasing to other positions) & horizontal
%     slownesses SLOW (to see how slowness error affects the positional
%     aliasing) at the given f0.  Latitudes & longitudes are expected in
%     degrees and are expected to be input as [LAT LON] arrays.  Positions
%     may also be given as [LAT LON ELEV DEPTH] where elevation and depth
%     are given in meters.  Generally you will set ELEV to 0 for the
%     beamforming positions in LATLON & lalo0.  Slownesses SLOW & slow0 are
%     the magnitude of the horizontal slowness in sec/deg.  SLOW & slow0
%     may also specify functions and this is described further in the Notes
%     section below.  The ARF is created using the 'center' method by
%     default as this is faster (the method can be altered -- see below).
%     The output struct S has the same fields as GEOFSS output but also
%     includes a S.source field which has the following layout:
%      .source.nsrc   - number of sources
%      .source.latlon - source positions as [LAT LON] (deg)
%      .source.slow   - source horizontal slowness (sec/deg)
%      .source.freq   - source frequency (Hz)
%
%     S=GEOARF(...,'METHOD',STRING,...) defines the beamforming method.
%     STRING may be 'center', 'coarray', 'full', or [LAT LON].  Note that
%     the 'capon' method is not available as the array response for that
%     method is dependent on the data.  The default is 'center'.
%
%     S=GEOARF(...,'WEIGHTS',W,...) specifies the relative weights for each
%     station (must have the same number of rows as STLALO) or pairing (if
%     METHOD is 'coarray' or 'full').  You may also specify weighting
%     based on pair geometry using 'ddiff' or 'azdiff' (but only if METHOD
%     is 'coarray' or 'full').  The default is even weighting.
%
%     S=GEOARF(...,'AVG',TRUE|FALSE,...) indicates if the spectra is
%     averaged across slowness & source during computation.  This can save
%     a significant amount of memory.  The default is false.
%
%    Notes:
%     - Latitudes are assumed to be geographic!
%     - Slowness inputs may also be specified as functions that output
%       travel time and phase shift given the following style of input:
%        [tt,shift]=func([evla evlo ...],[stla stlo ...],freq)
%       where tt is the travel time in seconds between the locations and
%       shift is the phase delay in radians (a Hilbert transform is -pi/2).
%       The lat/lon inputs are Nx2+ arrays as [lat lon], [lat lon elev] or
%       [lat lon elev depth].  The frequency input must be a scalar.  The
%       outputs are therefore Nx1 arrays.  This allows GEOARF to show the
%       source resolution from beamforming of just about any wave.
%
%    Examples:
%     % Global geoARF for a random array:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     st=randlatlon(100);
%     s=geoarf(st,[lat(:) lon(:)],25,[0 0],25,1/500);
%     plotgeoarf(s);
%
%     % Multi-source geoARF for a global array:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     st=randlatlon(100);
%     s=geoarf(st,[lat(:) lon(:)],25,randlatlon(5),25,1/500);
%     plotgeoarf(geofssavg(s));
%
%     % The array response of LASA to a 3s source near Lake Fort Peck:
%     src=[47.9 -106.5];
%     [lat,lon]=meshgrid(45:.025:48.5,-110:.025:-104.5);
%     plotgeofss(geoarf(lasa,[lat(:) lon(:)],33,src,33,1/3),...
%         [-3 0],[],'g','i');
%     % Now compare that to an adaptively weighted version:
%     plotgeofss(geoarf(lasa,[lat(:) lon(:)],33,src,33,1/3,...
%         'm','coarray','w','azdiff'),[-3 0],[],'g','i');
%
%    See also: PLOTGEOARF, GEOFSS, GEOFSSXC, GEOFSSHORZ, GEOFSSHORZXC,
%              GEOFSSAVG, GEOFSSSUB, PLOTGEOFSS, GEOFSSFRAMESLIDE,
%              GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE, GEOFSSCORR, GEOFSSDBINFO

%     Version History:
%        July  7, 2010 - initial version
%        Nov. 18, 2010 - add weighting
%        June 12, 2012 - adapt from geofkarf, alter output format, default
%                        to center method, always output struct, don't
%                        average the sources
%        Sep. 28, 2012 - pv pair inputs, doc update, tt functions
%        Sep. 29, 2012 - special weighting 'azdiff' & 'ddiff', avg option
%        Oct.  2, 2012 - minor doc update, minor bugfix for func handle,
%                        variable collision bugfix
%        Jan. 15, 2013 - history update
%        Mar. 21, 2014 - minor doc update, phase shift now also expected to
%                        be returned by function handles
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 21, 2014 at 17:00 GMT

% todo:
% - make functions to test
%   - teleseismic & regional P (travel time curve based)
%   - surface waves
% - better warnings from checking...

% check nargin
error(nargchk(6,inf,nargin));

% check inputs
if(~isreal(stlalo) || ndims(stlalo)~=2 || size(stlalo,2)<2)
    error('seizmo:geoarf:badInput',...
        'STLALO must be a Nx2 real array of [STLA STLO] in deg!');
elseif(~isreal(lalo) || ndims(lalo)~=2 ...
        || size(lalo,2)<2 || size(lalo,1)<1)
    error('seizmo:geoarf:badInput',...
        'LATLON must be a Mx2 real array of [LAT LON] in deg!');
elseif(~isa(slow,'function_handle') ...
        && (~isreal(slow) || any(slow<=0) || numel(slow)<1))
    error('seizmo:geoarf:badInput',...
        'SLOW must be positive real vector in s/deg!');
elseif(~isa(slow0,'function_handle') ...
        && (~isreal(slow0) || any(slow0<=0) || numel(slow0)<1))
    error('seizmo:geoarf:badInput',...
        'slow0 must be positive slowness in s/deg!');
elseif(~isreal(lalo0) || ndims(lalo0)~=2 ...
        || size(lalo0,2)<2 || size(lalo0,1)<1)
    error('seizmo:geoarf:badInput',...
        'latlon0 must be a NSRCx2 real array of [LAT LON] in deg!');
elseif(~isreal(f0) || any(f0<=0) || numel(f0)<1)
    error('seizmo:geoarf:badInput',...
        'f0 must be positive frequency in Hz!');
elseif(~isequalsizeorscalar(slow0(:),lalo0(:,1),f0(:)))
    error('seizmo:geoarf:badInput',...
        'slow0(:), latlon0(:,1), f0(:) must be equal sized or scalar!');
end

% slowness function logicals
sisafunc=isa(slow,'function_handle');
s0isafunc=isa(slow0,'function_handle');

% number of stations
nrecs=size(stlalo,1);

% need 2+ stations
if(nrecs<2)
    error('seizmo:geoarf:arrayTooSmall',...
        'GEOARF requires 2+ stations!');
end

% parse options
pv=parse_geoarf_pv_pairs(varargin{:});

% defaults for optionals
if(isempty(pv.w)); pv.w=ones(nrecs,1); end

% fix lat/lon
[stlalo(:,1),stlalo(:,2)]=fixlatlon(stlalo(:,1),stlalo(:,2));
[lalo0(:,1),lalo0(:,2)]=fixlatlon(lalo0(:,1),lalo0(:,2));
[lalo(:,1),lalo(:,2)]=fixlatlon(lalo(:,1),lalo(:,2));
nll=size(lalo,1);

% expand spherical wave details
nsrc=max([numel(slow0) size(lalo0,1) numel(f0)]);
if(isscalar(slow0) && ~s0isafunc); slow0(1:nsrc,1)=slow0; end
if(isscalar(f0)); f0(1:nsrc,1)=f0; end
if(size(lalo0,1)==1); lalo0=lalo0(ones(nsrc,1),:); end
if(~s0isafunc); slow0=slow0(:); end % force column vector
f0=f0(:); % force column vector

% fix method/center/npairs
if(ischar(pv.method))
    method=lower(pv.method);
    [clalo(1),clalo(2)]=arraycenter(stlalo(:,1),stlalo(:,2));
    clalo=[clalo mean(stlalo(:,3:end))];
    switch pv.method
        case 'coarray'
            [master,slave]=find(triu(true(nrecs),1));
            npairs=nrecs*(nrecs-1)/2;
        case 'full'
            [master,slave]=find(true(nrecs));
            npairs=nrecs*nrecs;
        case 'center'
            npairs=nrecs;
    end
else
    clalo=pv.method;
    pv.method='user';
    npairs=nrecs;
end

% check weights again
if(isnumeric(pv.w) && ~any(numel(pv.w)==[nrecs npairs]))
    error('seizmo:geoarf:badInput',...
        'Number of WEIGHTS must match the number of stations or pairs!');
end

% column vector slownesses
if(sisafunc)
    nslow=1;
else % explicit
    slow=slow(:);
    nslow=numel(slow);
end

% create geoARF struct
s.nsta=nrecs;
s.st=stlalo;
s.butc=[0 0 0 0 0];
s.eutc=[0 0 0 0 0];
s.delta=nan;
s.npts=0;
s.latlon=lalo;
s.slow=slow;
s.freq=f0;
s.method=pv.method;
s.npairs=npairs;
s.center=clalo;
s.whiten=true;
s.weights=pv.w;
s.source.nsrc=nsrc;
s.source.latlon=lalo0;
s.source.slow=slow0;
s.source.freq=f0;
if(pv.avg); s.spectra=zeros(nll,1,'single');
else s.spectra=zeros(nll,nslow,nsrc,'single');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of defaults, checks, & data setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% expand & normalize weights
if(isnumeric(pv.w))
    pv.w=pv.w(:);
    switch pv.method
        case {'coarray' 'full'}
            if(numel(pv.w)==nrecs)
                pv.w=pv.w(slave).*conj(pv.w(master));
            end
    end
    pv.w=pv.w./sum(abs(pv.w));
    z=1;
else
    switch pv.method
        case {'user' 'center'}
            error('seizmo:geoarf:badInput',...
                'Cannot use a pair weighting scheme for solo methods!');
    end
    switch pv.w
        case 'azdiff'
            [az,az]=sphericalinv(...
                lalo(:,ones(nrecs,1)),lalo(:,2*ones(nrecs,1)),...
                stlalo(:,ones(nll,1))',stlalo(:,2*ones(nll,1))');
            pv.w=ones(npairs,1);
            z=abs(azdiff(az(:,master),az(:,slave)));
            z=z./repmat(sum(z,2),1,npairs);
        case 'ddiff'
            pv.w=sphericalinv(...
                stlalo(master,1),stlalo(master,2),...
                stlalo(slave,1),stlalo(slave,2));
            pv.w=pv.w./sum(abs(pv.w));
            z=1;
    end
end

% detail message
verbose=seizmoverbose;
if(verbose)
    fprintf('Getting geoARF for %d spherical wave sources\n',nsrc);
    print_time_left(0,nsrc);
end

% loop over spherical wave sources
for a=1:nsrc
    % Get distance or traveltime from each
    % position of interest to every station
    % tt is NLLxNPAIRS
    if(a==1 || ttredo)
        [tt,ctt,shift,cshift,ttredo]=gettt(...
            lalo,stlalo,clalo,slow,f0,a,nll,nrecs,sisafunc,pv);
        % get degree distance difference or traveltime difference
        switch pv.method
            case {'coarray' 'full'}
                tt=tt(:,slave)-tt(:,master);
                shift=shift(:,slave)-shift(:,master);
            otherwise % {'center' 'user'}
                tt=tt-ctt(:,ones(1,nrecs));
                shift=shift-cshift(:,ones(1,nrecs));
        end
    end
    
    % Get traveltime from source to every station
    % tt0 is 1xNPAIRS
    [tt0,ctt0,shift0,cshift0]=gettt0(...
        lalo0(a,:),stlalo,clalo,slow0,f0,a,nrecs,s0isafunc,pv);
    % get traveltime difference
    switch pv.method
        case {'coarray' 'full'}
            tt0=tt0(:,slave)-tt0(:,master);
            shift0=shift0(:,slave)-shift0(:,master);
        otherwise % {'center' 'user'}
            tt0=tt0-ctt0;
            shift0=shift0-cshift0;
    end
    
    % loop over slownesses
    for b=1:nslow
        % modify distance to traveltime?
        if(sisafunc); hs=1; else hs=slow(b); end
        
        % beamforming method
        switch method
            case {'full' 'coarray'}
                if(pv.avg)
                    s.spectra=s.spectra+real((z.*exp(-2*pi*1i*f0(a)...
                        *(hs*tt-tt0(ones(nll,1),:))...
                        -1i*(shift-shift0(ones(nll,1),:))))*pv.w);
                else
                    s.spectra(:,b,a)=real((z.*exp(-2*pi*1i*f0(a)...
                        *(hs*tt-tt0(ones(nll,1),:))...
                        -1i*(shift-shift0(ones(nll,1),:))))*pv.w);
                end
            otherwise % {'center' 'user'}
                if(pv.avg)
                    s.spectra=s.spectra+abs((z.*exp(-2*pi*1i*f0(a)...
                        *(hs*tt-tt0(ones(nll,1),:))...
                        -1i*(shift-shift0(ones(nll,1),:))))*pv.w).^2;
                else
                    s.spectra(:,b,a)=abs((z.*exp(-2*pi*1i*f0(a)...
                        *(hs*tt-tt0(ones(nll,1),:))...
                        -1i*(shift-shift0(ones(nll,1),:))))*pv.w).^2;
                end
        end
    end
    if(verbose); print_time_left(a,nsrc); end
end
if(pv.avg); s.spectra=s.spectra/(nsrc*nslow); end

end


function [tt,ctt,shift,cshift,ttredo]=gettt(...
    lalo,stlalo,clalo,slow,f0,a,nll,nrecs,sisafunc,pv)
% super overcomplex distance/traveltime retriever
ctt=[];
cshift=[];
ttredo=true;
if(a==1)
    % function gives traveltime but numeric input gives slowness
    if(sisafunc)
        % travel time from function
        if(isscalar(unique(f0)))
            % only one frequency so only need to get tt once
            tt=nan(nll,nrecs);
            shift=nan(nll,nrecs);
            for i=1:nrecs
                [tt(:,i),shift(:,i)]=slow(...
                    lalo,stlalo(i*ones(nll,1),:),f0(1));
            end
            switch pv.method
                case {'user' 'center'}
                    [ctt,cshift]=slow(lalo,clalo,f0(1));
            end
            ttredo=false;
        else % multiple frequencies
            try
                % see if function depends on f
                [tt,shift]=slow(lalo(1,:),stlalo(1,:)); %#ok<NASGU>
                
                % still here?
                % guess that means frequency is ignored...
                tt=nan(nll,nrecs);
                shift=nan(nll,nrecs);
                for i=1:nrecs
                    [tt(:,i),shift(:,i)]=slow(...
                        lalo,stlalo(i*ones(nll,1),:));
                end
                switch pv.method
                    case {'user' 'center'}
                        [ctt,cshift]=slow(lalo,clalo);
                end
                ttredo=false;
            catch
                % slow is a function of frequency
                tt=nan(nll,nrecs);
                shift=nan(nll,nrecs);
                for i=1:nrecs
                    [tt(:,i),shift(:,i)]=slow(...
                        lalo,stlalo(i*ones(nll,1),:),f0(a));
                end
                switch pv.method
                    case {'user' 'center'}
                        [ctt,cshift]=slow(lalo,clalo,f0(a));
                end
            end
        end
    else
        % degree distance (scaled by the slowness gives the traveltime)
        tt=sphericalinv(...
            lalo(:,ones(nrecs,1)),lalo(:,2*ones(nrecs,1)),...
            stlalo(:,ones(nll,1))',stlalo(:,2*ones(nll,1))');
        shift=zeros(size(tt));
        ctt=sphericalinv(lalo(:,1),lalo(:,2),clalo(1),clalo(2));
        cshift=zeros(size(ctt));
        ttredo=false;
    end
else
    % slow is a function of frequency
    tt=nan(nll,nrecs);
    shift=nan(nll,nrecs);
    for i=1:nrecs
        [tt(:,i),shift(:,i)]=slow(lalo,stlalo(i*ones(nll,1),:),f0(a));
    end
    switch pv.method
        case {'user' 'center'}
            [ctt,cshift]=slow(lalo,clalo,f0(a));
    end
end

end


function [tt0,ctt0,shift0,cshift0]=gettt0(...
    lalo0,stlalo,clalo,slow0,f0,a,nrecs,s0isafunc,pv)
% super overcomplex distance/traveltime retriever
% function gives traveltime but numeric input gives slowness
ctt0=[];
cshift0=[];
if(s0isafunc)
    [tt0,shift0]=slow0(lalo0(ones(nrecs,1),:),stlalo,f0(a));
    tt0=tt0'; shift0=shift0';
    switch pv.method
        case {'user' 'center'}
            [ctt0,cshift0]=slow0(lalo0,clalo,f0(a));
    end
else
    % degree distance (scaled by the slowness gives the traveltime)
    tt0=slow0(a)*sphericalinv(...
        lalo0(:,1),lalo0(:,2),stlalo(:,1),stlalo(:,2))';
    shift0=zeros(size(tt0));
    ctt0=slow0(a)*sphericalinv(lalo0(:,1),lalo0(:,2),clalo(1),clalo(2));
    cshift0=zeros(size(ctt0));
end
end


function [pv]=parse_geoarf_pv_pairs(varargin)

% defaults
pv.avg=false;
pv.method='center';
pv.w=[];

% require pv pairs
if(mod(nargin,2))
    error('seizmo:geoarf:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:geoarf:badInput',...
        'Parameters must be specified using a string!');
end
for i=1:2:nargin
    switch varargin{i}
        case {'method' 'meth' 'm'}
            pv.method=varargin{i+1};
        case {'weights' 'w' 'wgt' 'wgts' 'weight'}
            pv.w=varargin{i+1};
        case {'average' 'av' 'avg' 'a'}
            pv.avg=varargin{i+1};
        otherwise
            error('seizmo:geoarf:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% valid method strings
valid.METHOD={'center' 'coarray' 'full'};
valid.W={'azdiff' 'ddiff'};

% check values
if((isnumeric(pv.method) ...
        && (~isreal(pv.method) || ~numel(pv.method)==2)) ...
        || (ischar(pv.method) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:geoarf:badInput',...
        ['METHOD must be one of the following:\n' ...
        '''CENTER'', ''COARRAY'', ''FULL'', or [LAT LON]!']);
elseif((~isnumeric(pv.w) && ~ischar(pv.w)) ...
        || (isnumeric(pv.w) && any(pv.w(:)<0)) ...
        || (ischar(pv.w) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:geoarf:badInput',...
        'WEIGHTS must be positive!');
elseif(~isscalar(pv.avg) || ~islogical(pv.avg))
    error('seizmo:geoarf:badInput',...
        'AVG must be TRUE or FALSE!');
end

end
