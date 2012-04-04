function [Kph,Kam,x,y]=rayleigh_2d_plane_wave_kernels(w,d,f,a,v,show)
%RAYLEIGH_2D_PLANE_WAVE_KERNELS    Makes Rayleigh wave sensitivity kernels
%
%    Usage:    [Kph,Kam,x,y]=rayleigh_2d_plane_wave_kernels(...
%                              width,spacing,freq,amp,velo)
%              [Kph,Kam,x,y]=rayleigh_2d_plane_wave_kernels(...
%                              width,spacing,freq,amp,velo,show)
%
%    Description:
%     [Kph,Kam,X,Y]=RAYLEIGH_2D_PLANE_WAVE_KERNELS(...
%                                              WIDTH,SPACING,FREQ,AMP,VELO)
%     returns phase and amplitude sensitivity kernels for 2D Rayleigh wave
%     phase velocities.  The kernels are based on a plane wave assumption.
%     That is, the sensitivities assume that the wave front has propogated
%     through the local region as a plane wave.  WIDTH gives the overall
%     width of the grid (ie the local region) in kilometers, SPACING gives
%     the node spacing of the grid in kilometers, FREQ & AMP contain the
%     frequencies and their corresponding weight in the kernel, and VELO is
%     the average Rayleigh wave phase velocity for the frequencies in FREQ.
%     The phase and amplitude kernels are returned in Kph and Kam,
%     respectively. The grid coordinates are in X and Y, which are in a
%     local, receiver-centered system where the positive x-axis points
%     directly away from the source and the positive y-axis points 90deg
%     clockwise from the source (along the plane wave front).
%
%     [Kph,Kam,X,Y]=RAYLEIGH_2D_PLANE_WAVE_KERNELS(...
%                                         WIDTH,SPACING,FREQ,AMP,VELO,SHOW)
%     toggles the plotting of the kernels using SHOW.  SHOW must be TRUE
%     (plots the kernels) or FALSE (no plotting, the default).
%
%    Notes:
%     - References:
%       Yang and Forsyth, 2006, Regional Tomographic Inversion Of The
%       Amplitude And Phase Of Rayleigh Waves With 2-D Sensitivity Kernels,
%       Geophys. J. Int., Vol. 166, pp. 1148-1160
%       
%       Zhou, Dahlen, and Nolet, 2004, Three-Dimensional Sensitivity
%       Kernels For Surface Wave Observables, Geophys. J. Int., Vol. 158,
%       pp. 142-168
%
%    Examples:
%     % Get frequency-amplitude values for a windowed 100s Rayleigh wave:
%     [f,a]=getmainlobe(1/100,1,[1000 2000],200/1000,[1000 1000]);
%
%     % Show the corresponding kernels assuming a Vph=4km/s:
%     rayleigh_2d_plane_wave_kernels(3000,10,f,a,4,true);
%
%    See also: GETMAINLOBE, SMOOTH2D, READKERNELS, WRITEKERNELS,
%              MAKEKERNELS, PLOTKERNELS

%     Version History:
%        Feb.  4, 2010 - rewrite and added documentation
%        July  9, 2010 - fixed see also section, fixed major bug (thanks to
%                        davidh), fixed nargchk
%        Mar. 24, 2012 - minor doc update, plot calls use handles now
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 16:30 GMT

% todo:

% check nargin
error(nargchk(5,6,nargin));

% check inputs
if(nargin==5 || isempty(show)); show=false; end
if(~isscalar(w) || ~isreal(w) || w<=0)
    error('seizmo:rayleigh_2d_plane_wave_kernels:badInput',...
        'WIDTH must be a positive real scalar (in kilometers)!');
elseif(~isscalar(d) || ~isreal(d) || d<=0)
    error('seizmo:rayleigh_2d_plane_wave_kernels:badInput',...
        'SPACING must be a positive real scalar (in kilometers)!');
elseif(isempty(f) || ~isreal(f) || any(f<=0))
    error('seizmo:rayleigh_2d_plane_wave_kernels:badInput',...
        'FREQ must be an array of positive reals (in Hz)!');
elseif(isempty(a) || ~isreal(a) || any(a<=0))
    error('seizmo:rayleigh_2d_plane_wave_kernels:badInput',...
        'AMP must be an array of positive reals (in Hz)!');
elseif(numel(f)~=numel(a))
    error('seizmo:rayleigh_2d_plane_wave_kernels:badInput',...
        'FREQ & AMP must have the same number of elements!');
elseif(~isscalar(v) || ~isreal(v) || v<=0)
    error('seizmo:rayleigh_2d_plane_wave_kernels:badInput',...
        'VELO must be a positive real scalar (in km/s)!');
elseif(~isscalar(show) || ~islogical(show))
    error('seizmo:rayleigh_2d_plane_wave_kernels:badInput',...
        'SHOW must be a scalar logical!');
end

% setup the grid
x=d:d:w/2;
x=[-x(end:-1:1) 0 x];
x0=x;
n=numel(x);
x=x(ones(n,1),:);
y=x';

% normalize amps
a=a/sum(a);

% make sensitivity kernel
k=2*pi*6371*f/v;
r=sqrt(x.^2+y.^2);
r(floor(n*n/2)+1)=sqrt(2)*d; % avoid divide by zero at receiver
c=d^2./(4*6371.^2*sqrt(2*pi*abs(sin(r./6371))));
Kph=zeros(n); Kam=zeros(n);
for i=1:numel(f)
    Kph=Kph-c.*a(i).*k(i).^1.5.*sin(k(i).*(x+r)./6371+pi/4);
    Kam=Kam-c.*a(i).*k(i).^1.5.*cos(k(i).*(x+r)./6371+pi/4);
end

% plotting
if(show)
    % use peak frequency
    [i,i]=max(a); p=1/f(i);
    
    % first plot the phase kernel
    fh(1)=figure;
    ax(1)=axes('parent',fh(1));
    imagesc(x0,x0,Kph,'parent',ax(1));
    colormap(ax(1),jet);
    xlabel(ax(1),'RADIAL POSITION (KM)');
    ylabel(ax(1),'TRANSVERSE POSITION (KM)');
    zlabel(ax(1),'SENSITIVITY (S/KM^2)');
    title(ax(1),{'2D PLANE-WAVE PHASE DELAY KERNEL' ...
        'FOR RAYLEIGH WAVE PHASE VELOCITIES' ...
        [sprintf('PERIOD: %gS ',p) sprintf('VELOCITY: %gKM/S',v)]});
    box(ax(1),'on');
    grid(ax(1),'on');
    colorbar('peer',ax(1));
    
    % now plot the amplitude kernel
    fh(2)=figure;
    ax(2)=axes('parent',fh(2));
    imagesc(x0,x0,Kam,'parent',ax(2));
    colormap(ax(2),jet);
    xlabel(ax(2),'RADIAL POSITION (KM)');
    ylabel(ax(2),'TRANSVERSE POSITION (KM)');
    zlabel(ax(2),'SENSITIVITY (S/KM^2)');
    title(ax(2),{'2D PLANE-WAVE AMPLITUDE DEVIATION KERNEL' ...
        'FOR RAYLEIGH WAVE PHASE VELOCITIES' ...
        [sprintf('PERIOD: %gS ',p) sprintf('VELOCITY: %gKM/S',v)]});
    box(ax(2),'on');
    grid(ax(2),'on');
    colorbar('peer',ax(2));
end

end
