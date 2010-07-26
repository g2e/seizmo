function [mp] = ColorSpiral(nca,npa,pfa);
%ColorSpiral: Generates a monotonic colormap with maximum color depth
%
%   [m] = ColorSpiral(n,np,pf);
%
%   nc   Number of colors (length of the colormap). Default = 64.
%   np   Number of sinusoidal periods. Default = 2.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   m    Color map.
%
%   This function returns an n x 3 matrix containing the RGB entries
%   used for colormaps in MATLAB figures. The colormap is designed
%   to have a monotonically increasing intensity, while maximizing
%   the color depth. This is achieved by generating a spiral through
%   the RGB cube that ranges from RBG = [0 0 0] to RGB = [1 1 1].
%
%   Example: Create a world map of the GEOID datat with the new color
%   scale. 
%      load geoid;
%      figure; 
%      worldmap(geoid,geoidlegend);
%      contourcmap([0:2.5:50],'ColorSpiral');
%      h = colorbar('SouthOutside');
%
%   J. McNames, "An effective color scale for simultaneous color and 
%   gray-scale publications," IEEE Signal Processing Magazine, in
%   press (January 2006).
%
%   Version 1.01 JM
%
%   See also colormap, jet, and caxis.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help ColorSpiral;
    return;
    end;
    
%====================================================================
% Process Function Arguments
%====================================================================    
nc = 64;                                                   % Default number of colors in the colormap
if exist('nca') & ~isempty(nca),
    nc = nca;
    end;

np = 2;                                                    % Default number of sinusoidal periods
if exist('npa') & ~isempty(npa),
    np = npa;
    end;    
    
pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Preprocessing
%====================================================================    
%wn = sqrt(3/8)*[0;triang(nc-2);0];                        % Triangular window function
wn = sqrt(3/8)*Hyperbola(nc);                              % Hyperbolic window function
a12 = asin(1/sqrt(3));                                     % First  rotation angle (radians)
a23 = pi/4;                                                % Second rotation angle (radians)

%====================================================================
% Main Routine
%==================================================================== 
t = linspace(sqrt(3),0,nc).';                              % Independent variable
r0 = t;                                                    % Initial red values = independent variable (t)
g0 = wn.*cos(((t-sqrt(3)/2)*np*2*pi/sqrt(3)));             % Initial green values = real part of complex sinusoid
b0 = wn.*sin(((t-sqrt(3)/2)*np*2*pi/sqrt(3)));             % Initial blue values = imaginary part of complex sinusoid

[ag,rd] = cart2pol(r0,g0);                                 % Convert to RG polar coordinates
[r1,g1] = pol2cart(ag+a12,rd);                             % First rotation & conversion back to cartesian coordiantes            
b1      = b0;

[ag,rd] = cart2pol(r1,b1);                                 % Convert RB to polar coordinates
[r2,b2] = pol2cart(ag+a23,rd);                             % Second rotation & conversion back to cartesian coordinates
g2      = g1;

%====================================================================
% Postprocessing
%====================================================================    
r  = max(min(r2,1),0);                                     % Make sure rotated color cubes don't excede the unit
g  = max(min(g2,1),0);                                     % color cube boundaries due to finite-precision effects
b  = max(min(b2,1),0);

mp = [r g b];                                              % The final colormap matrix

%====================================================================
% Plot Default Figure
%====================================================================
if pf,
    figure;
    h = plot([mp,sum(mp,2)/3]);
    set(h(1),'Color','r');
    set(h(2),'Color','g');
    set(h(3),'Color','b');
    set(h(4),'Color','k');
    set(h,'LineWidth',1.5);
    xlim([1 nc]);
    ylim([0 1]);
    box off;
    xlabel('Map Index');
    legend('Red','Green','Blue','Intensity');
    if nc<=256,                                            % MATLAB doesn't display colormaps with more than 256 colors correctly
        colormap(mp);
        colorbar;
        end;
    ylim([0 1.03]);
    end;
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,                                             % If no output arguments, don't return anything.
    clear('mp');
    end;    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Hyperbola Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [y] = Hyperbola(x,ymaxa,pfa);
%Hyperbola: Generates a hyperbolic window function
%
%   [y] = Hyperbola(x);
%
%   x    If scalar, window length. If vector, indices of window.
%   ymax Maximum value of the window amplitude. Default = 0.95.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   y    Window.
%
%   This function returns a vector that represents a hyperbolic 
%   window function. Visually, this is very similar to a triangular
%   or tent window function. However, the hyperbola is analytic
%   (all it's derivatives exist at all points) and has a rounded
%   peak. The parameter ymax controls how rounded the peak is.
%   This function is used in the ColorSpiral colormap to prevent
%   a discontinuity at the midpoint of the colormap.
%
%   Example: Generate the spectrogram of an intracranial pressure
%   signal using a Hyperbola window that is 45 s in duration.
%
%      load ICP.mat; 
%      icpd = decimate(icp,15);
%      wl   = round(45*fs/15);
%      Spectrogram(icpd,fs/15,Hyperbola(wl));
%
%   C. H. Edwards, D. E. Penney, "Calculus and Analytic Geometry,"
%   2nd edition, Prentice-Hall, 1986.
%
%   Version 1.00 JM
%
%   See also triang, window, and ColorSpiral.

%   See http://mathworld.wolfram.com/Hyperbola.html for details. 
   
%====================================================================
% Process Function Arguments
%====================================================================     
if length(x)==1,                                           % If is an integer, make it into an array
    x = 1:x;
    end;

ymax = 0.95;
if exist('ymaxa') & ~isempty(ymaxa),
    ymax = ymaxa;
    end;    

pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Preprocessing
%====================================================================  
a    = sqrt((1-ymax).^2/(1-(1-ymax).^2));                  % Pick a to obtain desired maximum
xmin = min(x);
xmax = max(x);
xs   = 2*(x-xmin)/(xmax-xmin) - 1;                         % Scale so it ranges from -1 to 1
nx   = length(x);

%====================================================================
% Main Routine
%====================================================================   
y = 1-sqrt(xs.^2+a^2)/sqrt(1+a^2);                         % Constrained to range from 0 to ymax

%====================================================================
% Postprocessing
%====================================================================   
y(y<0) = 0;                                                % Make sure y is not negatative due to finite precision effects
y      = y(:);                                             % Convert into a column vector
