function []=makekernels(f0,v,fs,swin,tprfrac,twin,w,d,lambda,path,post)
%MAKEKERNELS    Makes sensitivity kernels for Yang & Forsyth codes
%
%    Usage:    makekernels(f0,velo,fs,swin,tprfrac,...
%                          twin,width,spacing,lambda)
%              makekernels(f0,velo,fs,swin,tprfrac,...
%                          twin,width,spacing,lambda,path,post)
%
%    Description:
%     MAKEKERNELS(F0,VELO,FS,SWIN,TPRFRAC,TWIN,WIDTH,SPACING,LAMBDA) builds
%     2D phase and amplitude kernels for Rayleigh wave phase velocities
%     using the specifications given and writes them to text files in the
%     current directory with the naming scheme detailed in the Notes
%     section below.  F0 indicates the particular frequency to which the
%     kernels correspond.  F0 may be an array of frequencies, in which case
%     kernels are built at each frequency and all subsequent arguments
%     should be either scalar (replicated to the same value at each
%     frequency) or have the same number of elements as the number of
%     frequencies NF (different value for each frequency).  VELO is the
%     corresponding phase velocity.  FS is the sampling frequency for the
%     windowed sinusoid record.  SWIN is a 1x2 or NFx2 array of [START END]
%     of the sinusoid window limits.  TPRFRAC is the fraction of the window
%     on each end to be tapered.  TWIN is the zero-padded window limits
%     (must encompass SWIN) and is a 1x2 or NFx2 array as well.  WIDTH is
%     the width of the windowed kernel grid in kilometers.  SPACING is the
%     node spacing of the grid in kilometers.  LAMBDA is the characteristic
%     smoothing length in kilometers to apply to the kernels.  LAMBDA may
%     be scalar, NFx1, or NFxNL where NL is the number of smoothing lengths
%     to be used on the kernels at each frequency.
%
%     MAKEKERNELS(F0,VELO,FS,SWIN,TPRFRAC,...
%                          TWIN,WIDTH,SPACING,LAMBDA,PATH,POST)
%     enables changing the path for writing the kernel files to PATH.  PATH
%     must already exist (no directories are created).  POST specifies a
%     string to append to all the kernel filenames.  Both PATH and POST may
%     be strings or cellstr arrays (of size NFxNL where NF is the number of
%     frequencies and NL is the number of smoothing lengths).
%
%    Notes:
%     - Filename format:  kernel.XXXsYYYkmZZZZv
%        where XXX is the period (1/f0)
%              YYY is the smoothing length
%              ZZZZ is the velocity in m/s
%     - More details can be found by reading the help info for related
%       functions in the See Also section below.
%
%    Examples:
%     % First, get some frequency bands:
%     bank=filter_bank([0.0055 0.055],'variable',0.2,0.1);
%
%     % Second, get some phase velocities at those frequencies:
%     phvel=prem_dispersion(bank(:,1));
%
%     % Third, use the beat length as a window width:
%     L_beat=1./(bank(:,3)-bank(:,2));
%
%     % Fourth, modify that window width by an empirically derived
%     % power law to get a more typical window width:
%     win=L_beat.*(2.5+1000*bank(:,1).^2);
%
%     % Now get the kernels for a windowed sinusoid sampled at 1Hz, tapered
%     % on the outer 20% of the window, and zero-padded to a length of
%     % 9000s starting at -2000s.  The kernels are sampled in a 5000km grid
%     % centered on the receiver, grid points are every 10km and a Gaussian
%     % smoother with a characteristic distance of 100km is applied.  The
%     % kernel names are appended with '.example'.
%     makekernels(band(:,1),phvel,1,win*[0 1],0.2,[-2000 7000],5000,10,...
%         100,[],'.example');
%
%    See also: WRITEKERNELS, READKERNELS, RAYLEIGH_2D_PLANE_WAVE_KERNELS,
%              GETMAINLOBE, SMOOTH2D, PLOTKERNELS

%     Version History:
%        Feb.  5, 2010 - rewrite and added documentation
%        July  9, 2010 - fixed nargchk, improved example, output filename
%                        includes velocity in m/s now
%        Apr.  2, 2012 - minor doc update
%        July 26, 2013 - fix bug in example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 26, 2013 at 22:10 GMT

% todo:

% check nargin
error(nargchk(9,11,nargin));

% number of frequencies
nf=numel(f0); nl=size(lambda,2);

% handle file options
if(nargin==9); path={'.'}; post={''}; end
if(nargin==10); post={''}; end
if(isempty(path)); path={'.'}; end
if(isempty(post)); post={''}; end
if(~iscell(path) && ischar(path)); path=cellstr(path); end
if(~iscell(post) && ischar(post)); post=cellstr(post); end

% check inputs
if(~isreal(f0) || any(f0<=0))
    error('seizmo:makekernels:badInput',...
        'F0 must be positive real (in Hz)!');
    
elseif(~any(numel(v)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'VELO must be scalar or equal sized with F0!');
elseif(~isreal(v) || any(v<=0))
    error('seizmo:makekernels:badInput',...
        'VELO must be positive real (in km/s)!');
    
elseif(~any(numel(fs)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'FS must be scalar or equal sized with F0!');
elseif(~isreal(fs) || any(fs<=0) || any(fs<=2*f0))
    error('seizmo:makekernels:badInput',...
        'FS must be positive real (in Hz) & FS>(2*F0) !');
    
elseif(size(swin,2)~=2 || ~isreal(swin) || ~any(size(swin,1)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'SWIN must be an 1x2 or Nx2 of [START END] (in sec)!');
    
elseif(~any(numel(tprfrac)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'TPRFRAC must be scalar or equal sized with F0!');
elseif(~isreal(tprfrac) || any(tprfrac<0) || any(tprfrac>1))
    error('seizmo:makekernels:badInput',...
        'TPRFRAC must be positive real from 0 to 1 !');
    
elseif(size(twin,2)~=2 || ~isreal(twin) || ~any(size(twin,1)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'TWIN must be an 1x2 or Nx2 of [START END] (in sec)!');
    
elseif(~any(numel(w)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'WIDTH must be scalar or equal sized with F0!');
elseif(~isreal(w) || any(w<=0))
    error('seizmo:makekernels:badInput',...
        'WIDTH must be positive real (in kilometers)!');
    
elseif(~any(numel(d)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'SPACING must be scalar or equal sized with F0!');
elseif(~isreal(d) || any(d<=0))
    error('seizmo:makekernels:badInput',...
        'SPACING must be positive real (in kilometers)!');
    
elseif(~any(size(lambda,1)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'LAMBDA must have 1 row or as many as F0 has elements!');
elseif(~isreal(lambda) || any(lambda(:)<0))
    error('seizmo:makekernels:badInput',...
        'LAMBDA must be positive real (in kilometers)!');
    
elseif(~iscellstr(path) || ~any(numel(path)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'PATH must be a string or a cellstr with equal size to F0!');
    
elseif(~iscellstr(post) || ~any(numel(post)==[1 nf]))
    error('seizmo:makekernels:badInput',...
        'POST must be a string or a cellstr with equal size to F0!');
end

% expand scalars
if(isscalar(v)); v=v(ones(nf,1),1); end
if(isscalar(v)); v=v(ones(nf,1),1); end
if(isscalar(fs)); fs=v(ones(nf,1),1); end
if(size(swin,1)==1); swin=swin(ones(nf,1),:); end
if(isscalar(tprfrac)); tprfrac=tprfrac(ones(nf,1),1); end
if(size(twin,1)==1); twin=twin(ones(nf,1),:); end
if(isscalar(w)); w=w(ones(nf,1),1); end
if(isscalar(d)); d=d(ones(nf,1),1); end
if(size(lambda,1)==1); lambda=lambda(ones(nf,1),:); end
if(isscalar(path))
    path=path(ones(nf,nl));
elseif(size(path,1)==1)
    path=path(ones(nf,1),:);
elseif(size(path,2)==1)
    path=path(:,ones(1,nl));
end
if(isscalar(post))
    post=post(ones(nf,nl));
elseif(size(path,1)==1)
    post=post(ones(nf,1),:);
elseif(size(path,2)==1)
    post=post(:,ones(1,nl));
end

% zeropadding
zpad=[swin(:,1)-twin(:,1)  twin(:,2)-swin(:,2)];
if(any(zpad<0))
    error('seizmo:makekernels:badInput',...
        'TWIN must encompass SWIN!');
end

% file separator
sep=filesep;

% verbosity
verbose=seizmoverbose;

% detail message
if(verbose)
    disp('Making 2D Plane Wave Kernels');
    print_time_left(0,nf);
end

% loop over frequencies
for i=1:nf
    [f,a]=getmainlobe(f0(i),fs(i),swin(i,:),tprfrac(i),zpad(i,:));
    [Kph,Kam,x,y]=rayleigh_2d_plane_wave_kernels(w(i),d(i),f,a,v(i));
    for j=1:nl
        if(lambda(i,j))
            Kph2=smooth2d(Kph,lambda(i,j)/d(i),[],'zeropad');
            Kam2=smooth2d(Kam,lambda(i,j)/d(i),[],'zeropad');
        end
        file=[path{i,j} sep 'kernel.' sprintf('%03d',round(1/f0(i))) ...
            's' sprintf('%03d',round(lambda(i,j))) 'km' ...
            sprintf('%04d',round(v(i)*1000)) 'v' post{i,j}];
        writekernels(file,Kph2,Kam2,x,y);
    end
    % detail message
    if(verbose); print_time_left(i,nf); end
end

end
