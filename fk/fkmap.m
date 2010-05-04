function [varargout]=fkmap(data,smax,spts,frng)
%FKMAP    Returns a map of energy in frequency-wavenumber space
%
%    Usage:    map=fkmap(data,smax,spts,frng)
%
%    Description: MAP=FKMAP(DATA,SMAX,SPTS,FRNG) calculates the energy
%     moving through an array in frequency-wavenumber space.  Actually, to
%     allow for easier interpretation and averaging across frequencies, the
%     energy is mapped into slowness space.  The array info and data are
%     derived from the SEIZMO struct DATA.  Make sure station location and
%     timing fields are set!  The range of the slowness space is given by
%     SMAX (in s/deg) and extends from -SMAX to SMAX for both East/West and
%     North/South directions.  SPTS controls the number of slowness points
%     for both directions (SPTSxSPTS grid).  FRNG gives the frequency range
%     to average over as [FREQLOW FREQHIGH] in Hz.  MAP is a struct
%     containing all relevant info and the slowness map itself.  The struct
%     layout is:
%          .map   - slowness map
%          .nsta  - number of stations utilized in making map
%          .stla  - station latitudes
%          .stlo  - station longitudes
%          .stel  - station elevations (surface)
%          .stdp  - station depths (from surface)
%          .butc  - UTC start time of data
%          .eutc  - UTC end time of data
%          .npts  - number of time points
%          .delta - number of seconds between each time point
%          .smax  - maximum slowness in map
%          .spts  - resolution in horizontal and vertical slowness
%          .freqs - frequencies averaged in making map
%     Calling FKMAP with no outputs will automatically plot the slowness
%     map using PLOTFKMAP.  
%
%    Notes:
%     - Records in DATA must have equal number of points, equal sample
%       spacing, the same start time (in absolute time), and be evenly
%       spaced time series records.  Use functions SYNCHRONIZE, SYNCRATES,
%       & INTERPOLATE to get the timing/sampling the same.
%
%     - The circles of the bull's eye in the plot correspond to several
%       surface/diffracted/head seismic phases (which appear depends on the
%       plot's maximum slowness):
%        Phase:       Rg    Lg    Sn    Pn    Sdiff  Pdiff  PKPcdiff
%        Vel (km/s):  3.0   4.0   4.5   8.15  13.3   25.1   54.0
%        S (s/deg):   37.0  27.8  24.7  13.6  8.36   4.43   2.06
%       The radial lines correspond to 30deg steps in backazimuth.
%
%    Examples:
%     Show slowness map for a dataset at about 50s periods:
%      fkmap(data,50,201,[1/51 1/49])
%
%    See also: FKARF, SNYQUIST, PLOTFKMAP, KXY2SLOWBAZ, SLOWBAZ2KXY

%     Version History:
%        May   3, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   3, 2010 at 21:50 GMT

% todo:

% check nargin
msg=nargchk(4,4,nargin);
if(~isempty(msg)); error(msg); end

% check struct
versioninfo(data,'dep');

% check inputs
sf=size(frng);
if(~isreal(smax) || ~isscalar(smax) || smax<=0)
    error('seizmo:fk:badInput',...
        'SMAX must be a positive real scalar in s/deg!');
elseif(~isscalar(spts) || fix(spts)~=spts || spts<=0)
    error('seizmo:fk:badInput',...
        'SPTS must be a positive scalar integer!');
elseif(~isreal(frng) || numel(sf)~=2 || sf(2)~=2 || any(frng(:)<=0))
    error('seizmo:fk:badInput',...
        'FRNG must be a Nx2 array of [FREQLOW FREQHIGH] in Hz!');
end
nrng=sf(1);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data);
    
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
    
    % verbosity
    verbose=seizmoverbose;
    
    % require evenly spaced time series
    iftype=getenumid(data,'iftype');
    leven=getlgc(data,'leven');
    if(any(strcmpi(leven,'false')))
        error('seizmo:fk:badLEVEN',...
            ['Record(s):\n' sprintf('%d ',find(strcmpi(leven,'false'))) ...
            '\nInvalid operation on unevenly sampled record(s)!']);
    elseif(any(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy')))
        error('seizmo:fk:badIFTYPE',...
            ['Record(s):\n' sprintf('%d ',...
            find(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy'))) ...
            '\nDatatype of record(s) in DATA must be Timeseries or XY!']);
    end
    
    % require all records to have equal npts, delta, b utc, and 1 cmp
    % - we could drop the b UTC requirement but that would require having a
    %   shift term for each record (might be useful for surface waves and
    %   large aperture arrays)
    [npts,delta,butc,eutc,ncmp,stla,stlo,stel,stdp]=getheader(data,...
        'npts','delta','b utc','e utc','ncmp','stla','stlo','stel','stdp');
    butc=cell2mat(butc); eutc=cell2mat(eutc);
    if(numel(unique(npts))~=1 || numel(unique(delta))~=1 ...
            || size(unique(butc,'rows'),1)~=1)
        error('seizmo:fk:badData',...
            'Records in DATA must have equal NPTS, DELTA, & B (UTC)!');
    elseif(any(ncmp)~=1)
        error('seizmo:fk:badNCMP',...
            'Records in DATA must be single component only!');
    end
    
    % check nyquist
    fnyq=1/(2*delta(1));
    if(any(frng>=fnyq))
        error('seizmo:fk:badFRNG',...
            ['FRNG frequencies must be under the nyquist frequency (' ...
            num2str(fnyq) ')!']);
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
    
    % calculate coarray (relative positions)
    % r=(x  ,y  )
    %     ij  ij
    %
    % x is km east
    % y is km north
    %
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
    d2r=pi/180;
    [dist,az]=vincentyinv(stla(:,ones(nrecs,1)),stlo(:,ones(nrecs,1)),...
        stla(:,ones(nrecs,1))',stlo(:,ones(nrecs,1))');
    az=az*d2r;
    x=dist.*sin(az);
    y=dist.*cos(az);
    
    % make wavenumber arrays
    d2km=6371*d2r;
    smax=smax/d2km;
    sx=-smax:2*smax/(spts-1):smax;
    sy=fliplr(sx)';
    kx=2*pi*1i*sx(ones(spts,1),:);
    ky=2*pi*1i*sy(:,ones(spts,1));

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
    [smap(1:nrng,1).smax]=deal(smax*d2km);
    [smap(1:nrng,1).spts]=deal(spts);
    
    % loop over frequency ranges
    for a=1:nrng
        % get frequencies
        fidx=find(f>=frng(a,1) & f<=frng(a,2));
        smap(a).freqs=f(fidx);
        nfreq=numel(fidx);
        
        % trim data to necessary freqs, normalize
        fd=data(fidx,:);
        fd=fd./sqrt(fd.*conj(fd));
        cd=conj(fd);
        
        % detail message
        if(verbose)
            fprintf('Getting fk Map for %g to %g Hz\n',...
                f(fidx(1)),f(fidx(end)));
            print_time_left(0,nrecs^2*nfreq);
        end
        
        % loop over frequencies
        smap(a).map=zeros(spts,spts);
        for b=1:nfreq
            % loop over station pairs
            for c=1:nrecs^2
                d=mod(c,nrecs);
                if(d==0); d=nrecs; end
                e=ceil(c/nrecs);
                smap(a).map=smap(a).map+fd(b,d)*cd(b,e)...
                    *exp(f(fidx(b))*(kx*x(c)+ky*y(c)));
                if(verbose)
                    print_time_left(c+nrecs^2*(b-1),nrecs^2*nfreq);
                end
            end
        end
        
        % convert to dB
        % - note the use of abs to handle slightly negative terms
        smap(a).map=10*log10(abs(real(smap(a).map))/(nfreq*nrecs^2));
        
        % normalize so max peak is at 0dB
        smap(a).map=smap(a).map-max(smap(a).map(:));
        
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
