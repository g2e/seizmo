function [varargout]=fkmap(data,smax,spts,frng,polar,center)
%FKMAP    Returns a map of energy in frequency-wavenumber space
%
%    Usage:    smap=fkmap(data,smax,spts,frng)
%              smap=fkmap(data,smax,spts,frng,polar)
%              smap=fkmap(data,smax,spts,frng,polar,center)
%
%    Description: SMAP=FKMAP(DATA,SMAX,SPTS,FRNG) calculates the energy
%     moving through an array in frequency-wavenumber space.  Actually, to
%     allow for easier interpretation and averaging across frequencies, the
%     energy is mapped into slowness space.  The array info and data are
%     derived from the SEIZMO struct DATA.  Make sure station location and
%     timing fields are set!  The range of the slowness space is given by
%     SMAX (in s/deg) and extends from -SMAX to SMAX for both East/West and
%     North/South directions.  SPTS controls the number of slowness points
%     for both directions (SPTSxSPTS grid).  FRNG gives the frequency range
%     to average over as [FREQLOW FREQHIGH] in Hz.  SMAP is a struct
%     containing relevant info and the slowness map itself.  The struct
%     layout is:
%          .response - slowness map
%          .nsta     - number of stations utilized in making map
%          .stla     - station latitudes
%          .stlo     - station longitudes
%          .stel     - station elevations (surface)
%          .stdp     - station depths (from surface)
%          .butc     - UTC start time of data
%          .eutc     - UTC end time of data
%          .npts     - number of time points
%          .delta    - number of seconds between each time point
%          .x        - east/west slowness or azimuth values
%          .y        - north/south or radial slowness values
%          .z        - frequency values
%          .polar    - true if slowness is sampled in polar coordinates
%          .center   - array center or method
%          .normdb   - what 0dB actually corresponds to
%          .volume   - true if frequency-slowness volume (false for FKMAP)
%
%     Calling FKMAP with no outputs will automatically plot the slowness
%     map using PLOTFKMAP.
%
%     SMAP=FKMAP(DATA,SMAX,SPTS,FRNG,POLAR) specifies if the slowness map
%     is sampled regularly in cartesian or polar coordinates.  Polar coords
%     are useful for slicing the volume by azimuth (pie slice) or slowness
%     (rings).  Cartesian coords (the default) samples the slowness space
%     regularly in the East/West & North/South directions and so exhibits
%     less distortion of the slowness space.
%
%     SMAP=FKMAP(DATA,SMAX,SPTS,FRNG,POLAR,CENTER) defines the array
%     center.  CENTER may be [LAT LON], 'center', 'coarray', or 'full'.
%     The default is 'coarray'.  The 'center' option finds the center
%     position of the array by averaging the station positions (using
%     ARRAYCENTER).  Both 'coarray' and 'full' are essentially centerless
%     methods using the relative positioning between every possible pairing
%     of stations in the array.  The 'full' method includes redundant and
%     same station pairings (and will always give poorer results compared
%     to 'coarray').
%
%    Notes:
%     - Records in DATA must have equal number of points, equal sample
%       spacing, the same start time (in absolute time), and be evenly
%       spaced time series records.  Use functions SYNCHRONIZE, SYNCRATES,
%       & INTERPOLATE to get the timing/sampling the same.
%
%    Examples:
%     Show slowness map for a dataset at about 50s periods:
%      fkmap(data,50,201,[1/51 1/49])
%
%    See also: FKARF, SNYQUIST, PLOTFKMAP, KXY2SLOWBAZ, SLOWBAZ2KXY

%     Version History:
%        May   3, 2010 - initial version
%        May   7, 2010 - only doing one triangle gives better response and
%                        takes less than half the time
%        May   8, 2010 - array math version (another big speed jump), added
%                        a couple options for doing rose vs grid slowness
%                        plots, also allow changing between coarray & a
%                        specified array center (specified is hugely faster
%                        but suffers in resolution)
%        May   9, 2010 - struct changes
%        May  10, 2010 - use checkheader more effectively
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  10, 2010 at 15:50 GMT

% todo:

% check nargin
msg=nargchk(4,6,nargin);
if(~isempty(msg)); error(msg); end

% define some constants
d2r=pi/180;
d2km=6371*d2r;

% check struct
versioninfo(data,'dep');

% defaults for optionals
if(nargin<5 || isempty(polar)); polar=false; end
if(nargin<6 || isempty(center)); center='coarray'; end

% valid center strings
valid.CENTER={'center' 'coarray' 'full'};

% check inputs
sf=size(frng);
if(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:fkmap:badInput',...
        'SMAX must be a positive real scalar in s/deg!');
elseif(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:fkmap:badInput',...
        'SPTS must be a positive scalar integer >2!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:fkmap:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
elseif(~isscalar(polar) || (~islogical(polar) && ~isnumeric(polar)))
    error('seizmo:fkmap:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif((isnumeric(center) && (~isreal(center) || ~numel(center)==2)) ...
        || (ischar(center) && ~any(strcmpi(center,valid.CENTER))))
    error('seizmo:fkmap:badInput',...
        'CENTER must be [LAT LON], ''CENTER'', ''COARRAY'' or ''FULL''');
end
nrng=sf(1);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'MULTIPLE_NPTS','ERROR',...
        'NONINTEGER_REFTIME','ERROR',...
        'UNSET_REFTIME','ERROR',...
        'OUTOFRANGE_REFTIME','ERROR',...
        'UNSET_ST_LATLON','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
end

% do fk analysis
try
    % number of records
    nrecs=numel(data);
    
    % need 2+ records
    if(nrecs<2)
        error('seizmo:fkmap:arrayTooSmall',...
            'DATA must have 2+ records!');
    end
    
    % verbosity
    verbose=seizmoverbose;
    
    % require all records to have equal npts, delta, b utc, and 1 cmp
    % - we could drop the b UTC requirement but that would require having a
    %   shift term for each record (might be useful for surface waves and
    %   large aperture arrays)
    [npts,delta,butc,eutc,stla,stlo,stel,stdp]=getheader(data,...
        'npts','delta','b utc','e utc','stla','stlo','stel','stdp');
    butc=cell2mat(butc); eutc=cell2mat(eutc);
    if(size(unique(butc,'rows'),1)~=1)
        error('seizmo:fkmap:badData',...
            'Records in DATA must have equal B (UTC)!');
    end
    
    % check nyquist
    fnyq=1/(2*delta(1));
    if(any(frng>=fnyq))
        error('seizmo:fkmap:badFRNG',...
            ['FRNG frequencies must be under the nyquist frequency (' ...
            num2str(fnyq) ')!']);
    end
    
    % setup output
    [smap(1:nrng,1).nsta]=deal(nrecs);
    [smap(1:nrng,1).stla]=deal(stla);
    [smap(1:nrng,1).stlo]=deal(stlo);
    [smap(1:nrng,1).stel]=deal(stel);
    [smap(1:nrng,1).stdp]=deal(stdp);
    [smap(1:nrng,1).butc]=deal(butc(1,:));
    [smap(1:nrng,1).eutc]=deal(eutc(1,:));
    [smap(1:nrng,1).delta]=deal(delta(1));
    [smap(1:nrng,1).npts]=deal(npts(1));
    [smap(1:nrng,1).polar]=deal(polar);
    [smap(1:nrng,1).center]=deal(center);
    [smap(1:nrng,1).volume]=deal(false);
    
    % fix center
    if(ischar(center))
        center=lower(center);
    else
        clat=center(1);
        clon=center(2);
        center='user';
    end
    
    % get frequencies
    nspts=2^(nextpow2(npts(1))+1);
    f=(0:nspts/2)/(delta(1)*nspts);  % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    data=records2mat(data);
    seizmoverbose(verbose);
    
    % get fft
    data=delta(1)*fft(data,nspts,1); % delta for Parseval's
    
    % get relative positions from center
    % r=(x  ,y  )
    %     ij  ij
    %
    % x is km east
    % y is km north
    %
    % r is 2xNR
    switch center
        case 'coarray'
            % centerless (make coarray)
            % [ r   r   ... r
            %    11  12      1N
            %   r   r   ... r
            %    21  22      2N
            %    .   .  .    .
            %    .   .   .   .
            %    .   .    .  .
            %   r   r   ... r   ]
            %    N1  N2      NN
            [dist,az]=vincentyinv(...
                stla(:,ones(nrecs,1))',stlo(:,ones(nrecs,1))',...
                stla(:,ones(nrecs,1)),stlo(:,ones(nrecs,1)));
            idx=triu(true(nrecs),1);
            dist=dist(idx);
            az=az(idx);
        case 'full'
            % centerless too 
            [dist,az]=vincentyinv(...
                stla(:,ones(nrecs,1))',stlo(:,ones(nrecs,1))',...
                stla(:,ones(nrecs,1)),stlo(:,ones(nrecs,1)));
        case 'center'
            % get array center
            [clat,clon]=arraycenter(stla,stlo);
            [dist,az]=vincentyinv(clat,clon,stla,stlo);
        otherwise % user
            % array center was specified
            [dist,az]=vincentyinv(clat,clon,stla,stlo);
    end
    az=az*d2r;
    r=[dist(:).*sin(az(:)) dist(:).*cos(az(:))]';
    nidx=size(r,2);
    clear dist az
    
    % make projection arrays
    % p=2*pi*i*s*r
    %
    % where s is the slowness vector s=(s ,s ) and is NSx2
    %                                    x  y
    % p is NSxNR
    smax=smax/d2km;
    if(polar)
        if(numel(spts)==2)
            bazpts=spts(2);
            spts=spts(1);
        else
            bazpts=181;
        end
        smag=(0:spts-1)/(spts-1)*smax;
        [smap(1:nrng,1).y]=deal(smag'*d2km);
        smag=smag(ones(bazpts,1),:)';
        baz=(0:bazpts-1)/(bazpts-1)*360*d2r;
        [smap(1:nrng,1).x]=deal(baz/d2r);
        baz=baz(ones(spts,1),:);
        p=2*pi*1i*[smag(:).*sin(baz(:)) smag(:).*cos(baz(:))]*r;
        clear smag baz
        [smap(1:nrng,1).response]=deal(zeros(spts,bazpts));
    else % cartesian
        spts=spts(1);
        sx=-smax:2*smax/(spts-1):smax;
        [smap(1:nrng,1).x]=deal(sx*d2km);
        [smap(1:nrng,1).y]=deal(fliplr(sx*d2km)');
        sx=sx(ones(spts,1),:);
        sy=fliplr(sx)';
        p=2*pi*1i*[sx(:) sy(:)]*r;
        clear sx sy
        [smap(1:nrng,1).response]=deal(zeros(spts));
    end
    
    % loop over frequency ranges
    for a=1:nrng
        % get frequencies
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
        smap(a).z=f(fidx);
        nfreq=numel(fidx);
        
        % warning if no frequencies
        if(~nfreq)
            warning('seizmo:fkmap:noFreqs',...
                'No frequencies within the range %g to %g Hz!',...
                frng(a,1),frng(a,2));
            continue;
        end
        
        % build complex spectra array
        % cs is NRxNF
        cs=data(fidx,:)';
        switch center
            case 'coarray'
                cs=cs./sqrt(cs.*conj(cs));
                [row,col]=find(triu(true(nrecs),1));
                cs=cs(row,:).*conj(cs(col,:));
            case 'full'
                cs=cs./sqrt(cs.*conj(cs));
                [row,col]=find(true(nrecs));
                cs=cs(row,:).*conj(cs(col,:));
        end
        
        % detail message
        if(verbose)
            fprintf('Getting fk Map for %g to %g Hz\n',...
                frng(a,1),frng(a,2));
            print_time_left(0,nfreq);
        end
        
        % loop over frequencies
        for b=1:nfreq
            % get response
            switch center
                case {'full' 'coarray'}
                    smap(a).response(:)=smap(a).response(:)...
                        +exp(f(fidx(b))*p)*cs(:,b);
                otherwise
                    smap(a).response(:)=smap(a).response(:)...
                        +10*log10(abs(exp(f(fidx(b))*p)*cs(:,b)).^2/nidx);
            end
            
            % detail message
            if(verbose); print_time_left(b,nfreq); end
        end
        
        % convert to dB
        switch center
            case {'full' 'coarray'}
                % using full and real here gives the exact plots of
                % Koper, Seats, and Benz 2010 in BSSA
                smap(a).response=...
                    10*log10(abs(real(smap(a).response))/(nfreq*nidx));
            otherwise
                smap(a).response=smap(a).response/nfreq;
        end
        
        % normalize so max peak is at 0dB
        smap(a).normdb=max(smap(a).response(:));
        smap(a).response=smap(a).response-smap(a).normdb;
        
        % plot if no output
        if(~nargout); plotfkmap(smap(a)); drawnow; end
    end
    
    % return struct
    if(nargout); varargout{1}=smap; end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end