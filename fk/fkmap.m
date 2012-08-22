function [varargout]=fkmap(data,smax,spts,frng,polar,method,w)
%FKMAP    Returns beamformer map in frequency-wavenumber space
%
%    Usage:    smap=fkmap(data,smax,spts,frng)
%              smap=fkmap(data,smax,spts,frng,polar)
%              smap=fkmap(data,smax,spts,frng,polar,method)
%              smap=fkmap(data,smax,spts,frng,polar,method,weights)
%
%    Description:
%     SMAP=FKMAP(DATA,SMAX,SPTS,FRNG) beamforms wave energy moving through
%     an array in frequency-wavenumber space.  Actually, to allow for
%     easier interpretation and averaging across frequencies, the energy is
%     mapped into slowness space.  The array info and data are derived from
%     the SEIZMO struct DATA.  Make sure station location and timing fields
%     are set!  The range of the slowness space is given by SMAX (in s/deg)
%     and extends from -SMAX to SMAX for both East/West and North/South
%     directions.  SPTS controls the number of slowness points for both
%     directions (SPTSxSPTS grid).  FRNG gives the frequency range to
%     average over as [FREQLOW FREQHIGH] in Hz.  SMAP is a struct
%     containing relevant info and the slowness map itself.  The struct
%     layout is:
%          .beam     - slowness beamforming map
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
%          .freq     - frequency values
%          .polar    - true if slowness is sampled in polar coordinates
%          .npairs   - number of pairs
%          .method   - beamforming method (center, coarray, full, user)
%          .center   - array center as [lat lon]
%          .normdb   - what 0dB actually corresponds to
%          .volume   - true if frequency-slowness volume (false for FKMAP)
%          .weights  - weights used in beamforming
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
%     SMAP=FKMAP(DATA,SMAX,SPTS,FRNG,POLAR,METHOD) defines the beamforming
%     method.  METHOD may be 'center', 'coarray', 'full', or [LAT LON].
%     The default is 'coarray' which utilizes information from all unique
%     record pairings in the beamforming and is the default.  The 'full'
%     method will utilize all possible pairings including pairing records
%     with themselves and pairing records as (1st, 2nd) & (2nd, 1st) making
%     this method quite redundant and slow.  The 'center' option only pairs
%     each record against the array center (found using ARRAYCENTER) and is
%     extremely fast for large arrays compared to the 'coarray' & 'full'
%     methods.  Both 'center' and 'full' methods give slightly degraded
%     results compared to 'coarray'.  Using [LAT LON] for method is
%     essentially the same as the 'center' method but uses the defined
%     coordinates as the center for the array.
%
%     SMAP=FKMAP(DATA,SMAX,SPTS,FRNG,POLAR,METHOD,WEIGHTS) specifies
%     weights for each station pair in the array (this depends on the
%     method) for use in beamforming.  For example, the 'center' method 
%     requires N weights whereas 'coarray' requires (N*N-N)/2 weights.  The
%     weights are normalized internally to sum to 1.  See the Examples
%     section of FKARF for coarray weight indexing.
%
%    Notes:
%     - Records in DATA must have equal number of points, equal sample
%       spacing, the same start time (in absolute time), and be evenly
%       spaced time series records.
%
%    Examples:
%     % Show slowness map for a dataset at about 50s periods:
%     fkmap(data,50,201,[1/51 1/49])
%
%    See also: FKARF, SNYQUIST, PLOTFKMAP, KXY2SLOWBAZ, SLOWBAZ2KXY,
%              FKVOLUME, FK4D, FKVOL2MAP

%     Version History:
%        May   3, 2010 - initial version
%        May   7, 2010 - only doing one triangle gives better beam and
%                        takes less than half the time
%        May   8, 2010 - array math version (another big speed jump), added
%                        a couple options for doing rose vs grid slowness
%                        plots, also allow changing between coarray & a
%                        specified array center (specified is hugely faster
%                        but suffers in resolution slightly)
%        May   9, 2010 - struct changes
%        May  10, 2010 - use checkheader more effectively
%        May  12, 2010 - fixed an annoying bug (2 wrongs giving the right
%                        answer), minor code touches while debugging
%                        fkxcvolume
%        May  24, 2010 - beam is now single precision
%        June 16, 2010 - fixed nargchk, better verbose msg, improved see
%                        also section, create beam as s.p. array
%        July  1, 2010 - high latitude fix, allocation fix
%        July  6, 2010 - major update to struct, doc update
%        Nov. 18, 2010 - added weighting
%        Apr.  3, 2012 - use seizmocheck
%        May  30, 2012 - pow2pad=0 by default
%        June 13, 2012 - add capon method (not working!)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 13, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(4,7,nargin));

% define some constants
d2r=pi/180;
d2km=6371*d2r;

% check struct
error(seizmocheck(data,'dep'));

% defaults for optionals
if(nargin<5 || isempty(polar)); polar=false; end
if(nargin<6 || isempty(method)); method='coarray'; end
if(nargin<7); w=[]; end

% valid center strings
valid.METHOD={'center' 'coarray' 'full' 'capon'};

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
elseif((isnumeric(method) && (~isreal(method) || ~numel(method)==2)) ...
        || (ischar(method) && ~any(strcmpi(method,valid.METHOD))))
    error('seizmo:fkmap:badInput',...
        'METHOD must be ''CENTER'', ''COARRAY'', ''FULL'', or [LAT LON]!');
elseif(~isempty(w) && (any(w(:)<0) || ~isreal(w) || sum(w(:))==0))
    error('seizmo:fkmap:badInput',...
        'WEIGHTS must be positive real values!');
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
    
    % rethrow error
    error(lasterror);
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
    
    % convert weights to column vector and normalize
    w=w(:);
    w=w./sum(w);
    
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
    
    % fix method/center/npairs
    if(ischar(method))
        method=lower(method);
        [clat,clon]=arraycenter(stla,stlo);
        switch method
            case 'coarray'
                npairs=nrecs*(nrecs-1)/2;
            case {'full' 'capon'}
                npairs=nrecs*nrecs;
            case 'center'
                npairs=nrecs;
        end
    else
        clat=method(1);
        clon=method(2);
        method='user';
        npairs=nrecs;
    end

    % check that number of weights is equal to number of pairs
    if(isempty(w)); w=ones(npairs,1)/npairs; end
    if(numel(w)~=npairs)
        error('seizmo:fkmap:badInput',...
            ['WEIGHTS must have ' num2str(npairs) ...
            ' elements for method: ' method '!']);
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
    [smap(1:nrng,1).npairs]=deal(npairs);
    [smap(1:nrng,1).method]=deal(method);
    [smap(1:nrng,1).center]=deal([clat clon]);
    [smap(1:nrng,1).volume]=deal(false);
    [smap(1:nrng,1).weights]=deal(w);
    
    % get frequencies
    nspts=2^nextpow2(npts(1));
    f=(0:nspts/2)/(delta(1)*nspts);  % only +freq
    
    % extract data (silently)
    seizmoverbose(false);
    data=records2mat(data);
    seizmoverbose(verbose);
    
    % get fft
    data=fft(data,nspts,1);
    
    % get relative positions for each pair
    % r=(x  ,y  )
    %     ij  ij
    %
    % position of j as seen from i
    % x is km east
    % y is km north
    %
    % r is a 2xNPAIRS matrix
    switch method
        case 'coarray'
            % [ r   r   ... r
            %    11  12      1N
            %   r   r   ... r
            %    21  22      2N
            %    .   .  .    .
            %    .   .   .   .
            %    .   .    .  .
            %   r   r   ... r   ]
            %    N1  N2      NN
            %
            % then we just use the
            % upper triangle of that
            [e,n]=geographic2enu(stla,stlo,0,clat,clon,0);
            e=e(:,ones(nrecs,1))'-e(:,ones(nrecs,1));
            n=n(:,ones(nrecs,1))'-n(:,ones(nrecs,1));
            idx=triu(true(nrecs),1);
            e=e(idx);
            n=n(idx);
        case {'full' 'capon'}
            % retain full coarray
            [e,n]=geographic2enu(stla,stlo,0,clat,clon,0);
            e=e(:,ones(nrecs,1))'-e(:,ones(nrecs,1));
            n=n(:,ones(nrecs,1))'-n(:,ones(nrecs,1));
        case 'center'
            % each record relative to array center
            [e,n]=geographic2enu(clat,clon,0,stla,stlo,0);
        otherwise % user
            % each record relative to define array center
            [e,n]=geographic2enu(clat,clon,0,stla,stlo,0);
    end
    r=[e(:) n(:)]';
    clear e n
    
    % get phasors corresponding to each wave speed+direction (slowness)
    % to slant stack in the frequency domain
    % p=2*pi*i*s*r
    %
    % where s is the slowness vector s=(s ,s ) and is NSLOWx2
    %                                    x  y
    %
    % x is sec/km east
    % y is sec/km north
    %
    % p is the projection of all slownesses onto all of the position
    % vectors (multiplied by 2*pi*i for efficiency reasons)
    %
    % p is NSLOWxNPAIRS
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
        clear r smag baz
        [smap(1:nrng,1).beam]=deal(zeros(spts,bazpts,'single'));
    else % cartesian
        spts=spts(1);
        sx=-smax:2*smax/(spts-1):smax;
        [smap(1:nrng,1).x]=deal(sx*d2km);
        [smap(1:nrng,1).y]=deal(fliplr(sx*d2km)');
        sx=sx(ones(spts,1),:);
        sy=fliplr(sx)';
        p=2*pi*1i*[sx(:) sy(:)]*r;
        clear r sx sy
        [smap(1:nrng,1).beam]=deal(zeros(spts,'single'));
    end
    
    % loop over frequency ranges
    for a=1:nrng
        % get frequencies
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
        smap(a).freq=f(fidx);
        nfreq=numel(fidx);
        
        % warning if no frequencies
        if(~nfreq)
            warning('seizmo:fkmap:noFreqs',...
                'No frequencies within the range %g to %g Hz!',...
                frng(a,1),frng(a,2));
            continue;
        end
        
        % build complex cross spectra array, cs
        % cs is normalized by the auto spectra.
        %
        % note that center/user do not require
        % an array but do need normalization
        %
        % cs is NPAIRSxNFREQ
        cs=data(fidx,:).';
        cs=cs./abs(cs);
        switch method
            case 'coarray'
                [row,col]=find(triu(true(nrecs),1));
                cs=cs(row,:).*conj(cs(col,:));
            case {'full' 'capon'}
                [row,col]=find(true(nrecs));
                cs=cs(row,:).*conj(cs(col,:));
        end
        
        % capon 1969 minimum length method for weighting
        % - invert the cross power spectra complex matrix
        switch method
            case 'capon'
                % detail message
                if(verbose)
                    fprintf('Inverting Cross Power Spectra Matrix\n');
                    print_time_left(0,nfreq);
                end
                for b=1:nfreq
                    cs(:,b)=reshape(...
                        pinv(reshape(cs(:,b),[nrecs nrecs])),[],1);
                    if(verbose); print_time_left(b,nfreq); end
                end
            otherwise
                % apply weighting
                cs=cs.*w(:,ones(1,nfreq));
        end
        
        % detail message
        if(verbose)
            fprintf('Getting fk Map %d for %g to %g Hz\n',...
                a,frng(a,1),frng(a,2));
            print_time_left(0,nfreq);
        end
        
        % loop over frequencies, adding them together
        for b=1:nfreq
            % form beam
            switch method
                case {'full' 'coarray' 'capon'}
                    smap(a).beam(:)=smap(a).beam(:)...
                        +exp(f(fidx(b))*p)*cs(:,b);
                otherwise
                    smap(a).beam(:)=smap(a).beam(:)...
                        +abs(exp(f(fidx(b))*p)*cs(:,b)).^2;
            end
            
            % detail message
            if(verbose); print_time_left(b,nfreq); end
        end
        
        % convert sum to mean and convert to dB
        switch method
            case {'full' 'coarray' 'capon'}
                % using full and real here gives the exact plots of
                % Koper, Seats, and Benz 2010 in BSSA
                smap(a).beam=...
                    10*log10(abs(real(smap(a).beam))/nfreq);
            otherwise
                smap(a).beam=10*log10(smap(a).beam/nfreq);
        end
        
        % normalize so max peak is at 0dB
        smap(a).normdb=max(smap(a).beam(:));
        smap(a).beam=smap(a).beam-smap(a).normdb;
        
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
    error(lasterror);
end

end
