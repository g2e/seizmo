function [Kph,Kam,x,y]=readkernels(file)
%READKERNELS    Reads sensitivity kernels for Yang & Forsyth codes
%
%    Usage:    [Kph,Kam,x,y]=readkernels(file)
%
%    Description:
%     [Kph,Kam,X,Y]=READKERNELS(FILE) reads in the kernels stored in FILE.
%     FILE should be a file on the system.  FILE is expected to be
%     formatted to work with Yang & Forsyth fortran routines.  FILE may be
%     empty to bring up a graphical file selection menu.  Kph, Kam, X, Y
%     are equal-sized numeric arrays giving the phase sensitivity kernel,
%     the amplitude sensitivity kernel, the x-direction (radial) position
%     and the y-direction (azimuthal or transverse) position respectively.
%     Checks are performed to assure the kernel file is properly formatted.
%
%    Notes:
%
%    Examples:
%     % Check that a kernel is formatted correctly:
%     readkernels('my.kernel');
%
%    See also: WRITEKERNELS, MAKEKERNELS, RAYLEIGH_2D_PLANE_WAVE_KERNELS,
%              GETMAINLOBE, SMOOTH2D, PLOTKERNELS

%     Version History:
%        Feb.  5, 2010 - initial version
%        July  9, 2010 - fixed nargchk
%        Feb. 11, 2011 - mass nargchk fix
%        Mar. 24, 2012 - minor doc update
%        Jan. 26, 2014 - abs path exist fix
%        Feb.  8, 2014 - use readtxt, fix warning ids
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  8, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% get text
if(nargin<1); file=[]; end
txt=readtxt(file,{'*.kernel;*.KERNEL;kernel.*;KERNEL.*;' ...
    'Kernel Files (*.kernel,KERNEL.*)';
    '*.kern;*.KERN;kern.*;KERN.*;' ...
    'Kern Files (*.kern,KERN.*)';
    '*.dat;*.DAT' 'DAT Files (*.dat,*.DAT)';
    '*.*' 'All Files (*.*)'});

% parse and convert to double
v=str2double(getwords(txt));

% check number of elements
nv=numel(v);
if(mod(nv-6,4) || nv<10)
    error('seizmo:readkernels:malformedKernel',...
        'File: %s\nKernel is malformed!',file);
end

% get header info
nx=v(1); bx=v(2); dx=v(3);
ny=v(4); by=v(5); dy=v(6);

% cut off header
v(1:6)=[];

% push values into properly oriented arrays
x=reshape(v(1:4:end),[ny nx]);
y=reshape(v(2:4:end),[ny nx]);
Kph=reshape(v(3:4:end),[ny nx]);
Kam=reshape(v(4:4:end),[ny nx]);

% check position consistency
dx2=unique(diff(x,1,2));
dy2=unique(diff(y,1,1));
if(~isscalar(dx2) || dx2<=0)
    error('seizmo:readkernels:badInput',...
        'X step size is not uniform or is <=0!');
elseif(~isscalar(dy2) || dy2<=0)
    error('seizmo:readkernels:badInput',...
        'Y step size is not uniform or is <=0!');
elseif(dx2~=dx || x(1)~=bx)
    error('seizmo:readkernels:badInput',...
        'X header info does not match data!');
elseif(dy2~=dy || y(1)~=by)
    error('seizmo:readkernels:badInput',...
        'Y header info does not match data!');
end

end
