function [ax]=plotkernels(Kph,Kam,x,y,ax)
%PLOTKERNELS  Plots Rayleigh wave sensitivity kernels
%
%    Usage:    plotkernels(Kph,Kam,x,y)
%              plotkernels(Kph,Kam,x,y,ax)
%              ax=plotkernel(...)
%
%    Description:
%     PLOTKERNELS(Kph,Kam,X,Y) plots sensitivity kernels made with
%     RAYLEIGH_2D_PLANE_WAVE_KERNELS or read in with READKERNELS.  The
%     kernels are plotted as images (previously they were plotted as
%     surfaces but this was really slow).
%
%     PLOTKERNELS(Kph,Kam,X,Y,AX) sets the axes to draw in.  This is useful
%     for subplots, guis, etc.  Note that AX must be a 2 element vector!
%
%     AX=PLOTKERNELS(...) returns the handles to the axes that the kernels
%     are plotted in.  AX is a 2 element vector.
%
%    Notes:
%     - Unfortunately this does not include the velocity & period in the
%       title of the figures as this info is not saved in the Yang &
%       Forsyth kernel format.
%
%    Examples:
%     % Read in kernels and plot them up:
%     [Kph,Kam,x,y]=readkernels();
%     plotkernels(Kph,Kam,x,y);
%
%    See also: READKERNELS, RAYLEIGH_2D_PLANE_WAVE_KERNELS, GETMAINLOBE,
%              SMOOTH2D, MAKEKERNELS, WRITEKERNELS

%     Version History:
%        July  9, 2010 - initial version
%        Feb. 11, 2011 - improve axes calls in plotting
%        Apr.  2, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 16:30 GMT

% todo:

% check nargin
error(nargchk(4,5,nargin));

% check handle
if(nargin<5 || isempty(ax) || numel(ax)~=2 || ~isreal(ax) ...
        || ~all(ishandle(ax)) || ~all(strcmp('axes',get(ax,'type'))))
    badax=true;
else
    badax=false;
end

% check inputs
if(ndims(Kph)~=2 || ~isnumeric(Kph))
    error('seizmo:plotkernels:badInput',...
        'Kph must be a 2D array of numeric values!');
elseif(ndims(Kam)~=2 || ~isnumeric(Kam) || ~isequal(size(Kam),size(Kph)))
    error('seizmo:plotkernels:badInput',...
        'Kam must be a numeric array equal in size to Kph!');
elseif(ndims(x)~=2 || ~isnumeric(x) || ~isequal(size(x),size(Kph)))
    error('seizmo:plotkernels:badInput',...
        'X must be a numeric array equal in size to Kph!');
elseif(ndims(y)~=2 || ~isnumeric(y) || ~isequal(size(y),size(Kph)))
    error('seizmo:plotkernels:badInput',...
        'Y must be a numeric array equal in size to Kph!');
end
dx=unique(diff(x,1,2));
dy=unique(diff(y,1,1));
if(~isscalar(dx) || dx<=0)
    error('seizmo:plotkernels:badInput',...
        'X step size is not uniform or is <=0!');
elseif(~isscalar(dy) || dy<=0)
    error('seizmo:plotkernels:badInput',...
        'Y step size is not uniform or is <=0!');
end

% get x0
x0=unique(x(:));

% first plot the phase kernel
if(badax)
    fh=figure;
    ax(1)=axes('parent',fh);
end
imagesc(x0,x0,Kph,'parent',ax(1));
colormap(ax(1),jet);
xlabel(ax(1),'RADIAL POSITION (KM)');
ylabel(ax(1),'TRANSVERSE POSITION (KM)');
zlabel(ax(1),'SENSITIVITY (S/KM^2)');
title(ax(1),{'2D PLANE-WAVE PHASE DELAY KERNEL' ...
    'FOR RAYLEIGH WAVE PHASE VELOCITIES'});
box(ax(1),'on');
grid(ax(1),'on');
colorbar('peer',ax(1));

% now plot the amplitude kernel
if(badax)
    fh=figure;
    ax(2)=axes('parent',fh);
end
imagesc(x0,x0,Kam,'parent',ax(2));
colormap(ax(2),jet);
xlabel(ax(2),'RADIAL POSITION (KM)');
ylabel(ax(2),'TRANSVERSE POSITION (KM)');
zlabel(ax(2),'SENSITIVITY (S/KM^2)');
title(ax(2),{'2D PLANE-WAVE AMPLITUDE DEVIATION KERNEL' ...
    'FOR RAYLEIGH WAVE PHASE VELOCITIES'});
box(ax(2),'on');
grid(ax(2),'on');
colorbar('peer',ax(2));

end
