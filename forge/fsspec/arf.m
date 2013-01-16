function [s]=arf(stlalo,smax,spts,baz0,slow0,f0,varargin)
%ARF    Returns the array response function for a seismic array
%
%    Usage:    s=arf(stlalo,smax,spts,baz0,slow0,f0)
%              s=arf(...,'polar',true|false,...)
%              s=arf(...,'method',string,...)
%              s=arf(...,'weights',w,...)
%              s=arf(...,'avg',true|false,...)
%
%    Description:
%     S=ARF(STLALO,SMAX,SPTS,baz0,slow0,f0) computes the array response
%     function (ARF) for an array at the locations in STLALO with plane
%     waves passing through the array defined by backazimuths baz0,
%     horizontal slownesses slow0, and frequencies f0.  The ARF is computed
%     in cartesian slowness space where the range of the sampling is given
%     by SMAX (sec/deg) and extends from -SMAX to SMAX for both East/West
%     and North/South directions.  SPTS controls the number of slowness
%     points for both directions (SPTSxSPTS grid).  Latitudes & longitudes
%     are expected in degrees with the position arrays organized as
%     [LAT LON].  The ARF is created using the 'center' method as this is
%     faster (see below to adjust the method).  The output struct S has the
%     same fields as FSS output but also includes a S.source field which
%     has the following layout:
%      .source.nsrc - number of plane wave sources
%      .source.baz  - backazimuth (deg)
%      .source.slow - horizontal slowness (sec/deg)
%      .source.freq - frequency (Hz)
%     By default SPTS is 101, baz0 is 0, slow0 is 0, & f0 is 1.
%
%     S=ARF(...,'POLAR',TRUE|FALSE,...) decides if the spectra of the array
%     response is sampled regularly with cartesian or polar coordinates.
%     Polar coords are useful for slicing the spectra by azimuth (eg. a pie
%     slice) or slowness (eg. rings).  Cartesian coords (the default)
%     samples the slowness space regularly in the East/West & North/South
%     directions and so exhibits less distortion in plots of the slowness
%     space. If POLAR=TRUE, SPTS may be given as [spnts bazpnts] to control
%     the azimuthal resolution (default is 181 points).
%
%     S=ARF(...,'METHOD',STRING,...) defines the beamforming method.
%     STRING may be 'center', 'coarray', 'full', or [LAT LON].  Note that
%     the 'capon' method is not available as the array response for that
%     method is dependent on the data.
%
%     S=ARF(...,'WEIGHTS',W,...) specifies the relative weights for each
%     station (must have the same number of rows as STLALO) or pairing (if
%     METHOD is 'coarray' or 'full').
%
%     S=ARF(...,'AVG',TRUE|FALSE,...) indicates if the spectra is averaged
%     across frequency during computation.  This can save a significant
%     amount of memory.  The default is false.
%
%    Notes:
%
%    Examples:
%     % Show the 1Hz array response function out to 5s/deg
%     % for the Large Aperture Seismic Array (LASA):
%     plotarf(arf(lasa,5));
%
%     % Get a multi-plane wave response at 0.03Hz for a random array:
%     plotarf(arf(randlatlon(20)/30,50,201,[0 0 45],[20 10 20],0.03,...
%         'a',true,'m','coarray'));
%
%    See also: PLOTARF, ARFHORZ, FSS, FSSXC, FSSHORZ, FSSHORZXC, SNYQUIST

%     Version History:
%        May   1, 2010 - initial version
%        May   3, 2010 - show slowness nyquist ring, doc update
%        May   4, 2010 - circle is its own function now
%        May   7, 2010 - only doing one triangle gives better beam and
%                        takes less than half the time
%        May  10, 2010 - added in options available to FKMAP
%        May  18, 2010 - minor doc touch
%        June 16, 2010 - fixed nargchk, improved see also section and notes
%        July  6, 2010 - major update to struct, doc update, high latitude
%                        fix, arf data is now single precision
%        July  7, 2010 - center/user method for multi-plane wave ARF now is
%                        fixed
%        Nov. 16, 2010 - added weighting (franklin witnessed)
%        Nov. 17, 2010 - forgot .weights field
%        Nov. 18, 2010 - coarray weight indexing example, fixed scaling bug
%        Apr. 27, 2012 - spts now allows input as [spts bazpts]
%        Sep. 14, 2012 - adapted from fkarf
%        Sep. 22, 2012 - allow station or pair weights
%        Sep. 27, 2012 - pv pair inputs, doc update
%        Sep. 30, 2012 - avg option
%        Jan.  9, 2013 - allow options to be any case
%        Jan. 14, 2013 - update history
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 14, 2013 at 14:40 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% check required inputs
if(~isreal(stlalo) || ndims(stlalo)>2 || size(stlalo,2)<2)
    error('seizmo:arf:badInput',...
        'STLALO must be a Nx2 array of [STLA STLO] in deg!');
elseif(~isscalar(smax) || ~isreal(smax) || smax<=0)
    error('seizmo:arf:badInput',...
        'SMAX must be a positive real scalar in sec/deg!');
end

% number of stations
nrecs=size(stlalo,1);

% need 2+ stations
if(nrecs<2)
    error('seizmo:arf:arrayTooSmall',...
        'ARF requires 2+ stations!');
end

% defaults
if(nargin<3 || isempty(spts)); spts=101; end
if(nargin<4 || isempty(baz0)); baz0=0; end
if(nargin<5 || isempty(slow0)); slow0=0; end
if(nargin<6 || isempty(f0)); f0=1; end

% check inputs
if(~any(numel(spts)==[1 2]) || any(fix(spts)~=spts) || any(spts<=2))
    error('seizmo:arf:badInput',...
        'SPTS must be a positive scalar integer >2!');
elseif(~isreal(baz0))
    error('seizmo:arf:badInput',...
        'baz0 must be in degrees!');
elseif(~isreal(slow0))
    error('seizmo:arf:badInput',...
        'slow0 must be horizontal slowness in s/deg!');
elseif(~isreal(f0) || any(f0<=0))
    error('seizmo:arf:badInput',...
        'f0 must be a positive frequency in Hz!');
elseif(~isequalsizeorscalar(slow0,baz0,f0))
    error('seizmo:arf:badInput',...
        'baz0, slow0, f0 must be equal sized or scalar!');
end

% parse options
pv=parse_arf_pv_pairs(varargin{:});

% defaults for optionals
if(isempty(pv.w)); pv.w=ones(nrecs,1); end

% fix backazimuths
baz0=lonmod(baz0);

% expand plane wave details
[baz0,slow0,f0]=expandscalars(baz0,slow0,f0);
npw=numel(f0);
slow0=slow0(:);
baz0=baz0(:);
f0=f0(:);

% fix negative slowness
baz0(slow0<0)=lonmod(baz0(slow0<0)+180);
slow0=abs(slow0);

% fix method/center/npairs
if(ischar(pv.method))
    pv.method=lower(pv.method);
    [clat,clon]=arraycenter(stlalo(:,1),stlalo(:,2));
    switch pv.method
        case 'coarray'
            npairs=nrecs*(nrecs-1)/2;
        case {'full' 'capon'}
            npairs=nrecs*nrecs;
        case 'center'
            npairs=nrecs;
    end
else
    clat=pv.method(1);
    clon=pv.method(2);
    pv.method='user';
    npairs=nrecs;
end

% check weights again
if(~any(numel(pv.w)==[nrecs npairs]))
    error('seizmo:arf:badInput',...
        'Number of WEIGHTS must match the number of stations or pairs!');
end

% setup slowness grid
if(pv.polar)
    if(numel(spts)==1); spts(2)=181; end % default # azimuthal points
    sx=(0:spts(2)-1)/(spts(2)-1)*360; % baz (wedge decided x/y)
    sy=(0:spts(1)-1).'/(spts(1)-1)*smax; % smag
else
    spts(2)=spts(1);
    sx=-smax:2*smax/(spts(1)-1):smax; % east
    sy=fliplr(sx).'; % north
end

% setup output
s.nsta=nrecs;
s.st=stlalo;
s.butc=[0 0 0 0 0];
s.eutc=[0 0 0 0 0];
s.delta=nan;
s.npts=0;
s.polar=pv.polar;
s.x=sx;
s.y=sy;
s.freq=f0;
s.method=pv.method;
s.npairs=npairs;
s.center=[clat clon];
s.whiten=true;
s.weights=pv.w;
s.source.nsrc=npw;
s.source.baz=baz0;
s.source.slow=slow0;
s.source.freq=f0;
if(pv.avg); s.spectra=zeros(spts,'single');
else s.spectra=zeros([spts npw],'single');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of defaults, checks, & data setup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% radius constants
d2r=pi/180;
d2km=6371*d2r;

% get number of slowness points
nslow=prod(spts);

% expand & normalize weights
pv.w=pv.w(:);
switch pv.method
    case 'coarray'
        % row is master index, column is slave index
        [master,slave]=find(triu(true(nrecs),1));
        if(numel(pv.w)==nrecs); pv.w=pv.w(slave).*conj(pv.w(master)); end
    case 'full'
        [master,slave]=find(true(nrecs));
        if(numel(pv.w)==nrecs); pv.w=pv.w(slave).*conj(pv.w(master)); end
end
pv.w=pv.w./sum(abs(pv.w));

% get relative positions for each pair (the coarray)
% r=(x  ,y  )
%     ij  ij
%
% x   is km east, y   is km north
%  ij              ij
%
% r is a 2xNPAIRS matrix
switch pv.method
    case {'coarray' 'full'}
        % [ r   r   ... r
        %    11  12      1N
        %   r   r   ... r
        %    21  22      2N
        %    .   .  .    .
        %    .   .   .   .
        %    .   .    .  .
        %   r   r   ... r   ]
        %    N1  N2      NN
        [e,n]=geographic2enu(stlalo(:,1),stlalo(:,2),0,clat,clon,0);
        e=e(slave)-e(master);
        n=n(slave)-n(master);
    otherwise
        % each station relative to array center
        [e,n]=geographic2enu(stlalo(:,1),stlalo(:,2),0,clat,clon,0);
end
r=[e(:) n(:)]';
clear e n

% Get phasors for each wave speed+direction (slowness) and pair
% which allows us to delay and sum in the frequency domain.
% p=s*r
%
% where r was defined above
%       s is the slowness vector s=(s ,s ) and is NSLOWx2
%                                    x  y
%       and s  is in sec/km east, s  is in sec/km north
%            x                     y
%
% p is the projection of the slowness vectors s onto the
% spatial difference vectors r (called the coarray)
%
% p is NSLOWxNPAIRS
if(pv.polar)
    sx=sx(ones(spts(1),1),:); % baz in degrees
    sy=sy(:,ones(spts(2),1))/d2km; % smag in sec/km
    [sx,sy]=deal(sy.*sind(sx),sy.*cosd(sx));
    p=[sx(:) sy(:)]*r;
else % cartesian
    sx=sx/d2km;
    sx=sx(ones(spts(1),1),:);
    sy=fliplr(sx)';
    p=[sx(:) sy(:)]*r;
end
clear sx sy

% slowness offsets based on plane wave locations
slow0=[slow0(:).*sin(d2r*baz0(:)) slow0(:).*cos(d2r*baz0(:))]/d2km;

% detail message
verbose=seizmoverbose;
if(verbose)
    fprintf('Getting ARF for %d plane wave sources\n',npw);
    print_time_left(0,npw);
end

% loop over plane waves
for a=1:npw
    % projection of plane wave
    p0=slow0(a,:)*r;
    p0=p0(ones(nslow,1),:);
    
    % beamform
    switch pv.method
        case {'center' 'user'}
            if(pv.avg)
                s.spectra=s.spectra+reshape(abs(...
                    exp(-2*pi*1i*f0(a)*(p-p0))*pv.w).^2,spts);
            else
                s.spectra(:,:,a)=reshape(abs(...
                    exp(-2*pi*1i*f0(a)*(p-p0))*pv.w).^2,spts);
            end
        otherwise % {'coarray' 'full'}
            if(pv.avg)
                s.spectra=s.spectra+reshape(real(...
                    exp(-2*pi*1i*f0(a)*(p-p0))*pv.w),spts);
            else
                s.spectra(:,:,a)=reshape(real(...
                    exp(-2*pi*1i*f0(a)*(p-p0))*pv.w),spts);
            end
    end
    
    % progress bar update
    if(verbose); print_time_left(a,npw); end
end
if(pv.avg); s.spectra=s.spectra/npw; end

end


function [pv]=parse_arf_pv_pairs(varargin)

% defaults
pv.avg=false;
pv.polar=false;
pv.method='center';
pv.w=[];

% require pv pairs
if(mod(nargin,2))
    error('seizmo:arf:badInput',...
        'Unpaired parameter/value!');
end

% parse pv pairs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:arf:badInput',...
        'Parameters must be specified using a string!');
end
for i=1:2:nargin
    switch lower(varargin{i})
        case {'polar' 'plr' 'pol' 'p'}
            pv.polar=varargin{i+1};
        case {'method' 'meth' 'm'}
            pv.method=varargin{i+1};
        case {'weights' 'w' 'wgt' 'wgts' 'weight'}
            pv.w=varargin{i+1};
        case {'average' 'av' 'avg' 'a'}
            pv.avg=varargin{i+1};
        otherwise
            error('seizmo:arf:badInput',...
                'Unknown parameter: %s !',varargin{i});
    end
end

% valid method strings
valid.METHOD={'center' 'coarray' 'full'};

% check values
if(~isscalar(pv.polar) || (~islogical(pv.polar) && ~isnumeric(pv.polar)))
    error('seizmo:arf:badInput',...
        'POLAR must be TRUE or FALSE!');
elseif((isnumeric(pv.method) ...
        && (~isreal(pv.method) || ~numel(pv.method)==2)) ...
        || (ischar(pv.method) && ~any(strcmpi(pv.method,valid.METHOD))))
    error('seizmo:arf:badInput',...
        ['METHOD must be one of the following:\n' ...
        '''CENTER'', ''COARRAY'', ''FULL'', or [LAT LON]!']);
elseif(~isnumeric(pv.w) || any(pv.w(:)<0))
    error('seizmo:arf:badInput',...
        'WEIGHTS must be positive!');
elseif(~isscalar(pv.avg) || ~islogical(pv.avg))
    error('seizmo:arf:badInput',...
        'AVG must be TRUE or FALSE!');
end

end