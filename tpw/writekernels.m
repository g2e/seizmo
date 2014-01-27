function []=writekernels(file,Kph,Kam,x,y,o)
%WRITEKERNELS    Writes out sensitivity kernels for Yang & Forsyth codes
%
%    Usage:    writekernels(file,Kph,Kam,x,y)
%              writekernels(file,Kph,Kam,x,y,overwrite)
%
%    Description:
%     WRITEKERNELS(FILE,KPH,KAM,X,Y) writes out kernel info created with
%     RAYLEIGH_2D_PLANE_WAVE_KERNELS to a file for use with Yang & Forsyth
%     fortran routines.  FILE should be a string giving the filename and
%     path.  FILE may be empty to bring up a graphical file creation menu.
%     Kph, Kam, X, Y should all be equal-sized numeric arrays (see
%     RAYLEIGH_2D_PLANE_WAVE_KERNELS for details).
%
%     WRITEKERNELS(FILENAME,KPH,KAM,X,Y,OVERWRITE) quietly overwrites the
%     pre-existing file without confirmation when OVERWRITE is set to TRUE.
%     By default OVERWRITE is FALSE.  OVERWRITE is ignored in the graphical
%     file creation menu.
%
%    Notes:
%
%    Examples:
%     % Create kernels and write them out:
%     [f,a]=getmainlobe(1/100,1,[1000 2000],200/1000,[1000 1000]);
%     [Kph,Kam,x,y]=rayleigh_2d_plane_wave_kernels(3000,10,f,a,4);
%     Kph=smooth2d(Kph,100/10,[],'zeropad');     % 100km characteristic
%     Kam=smooth2d(Kam,100/10,[],'zeropad');     %   falloff distance
%     writekernels([],Kph,Kam,x,y); % choose your own filename
%
%    See also: READKERNELS, RAYLEIGH_2D_PLANE_WAVE_KERNELS, GETMAINLOBE,
%              SMOOTH2D, MAKEKERNELS, PLOTKERNELS

%     Version History:
%        Feb.  5, 2010 - rewrite and added documentation
%        July  9, 2010 - fixed nargchk
%        Feb. 11, 2011 - use fprintf
%        Apr.  2, 2012 - minor doc update
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 17:10 GMT

% todo:

% check nargin
error(nargchk(5,6,nargin));

% check inputs
if(nargin==5 || isempty(o)); o=false; end
if(ndims(Kph)~=2 || ~isnumeric(Kph))
    error('seizmo:writekernels:badInput',...
        'Kph must be a 2D array of numeric values!');
elseif(ndims(Kam)~=2 || ~isnumeric(Kam) || ~isequal(size(Kam),size(Kph)))
    error('seizmo:writekernels:badInput',...
        'Kam must be a numeric array equal in size to Kph!');
elseif(ndims(x)~=2 || ~isnumeric(x) || ~isequal(size(x),size(Kph)))
    error('seizmo:writekernels:badInput',...
        'X must be a numeric array equal in size to Kph!');
elseif(ndims(y)~=2 || ~isnumeric(y) || ~isequal(size(y),size(Kph)))
    error('seizmo:writekernels:badInput',...
        'Y must be a numeric array equal in size to Kph!');
elseif(~isscalar(o) || ~islogical(o))
    error('seizmo:writekernels:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end
dx=unique(diff(x,1,2));
dy=unique(diff(y,1,1));
if(~isscalar(dx) || dx<=0)
    error('seizmo:writekernels:badInput',...
        'X step size is not uniform or is <=0!');
elseif(~isscalar(dy) || dy<=0)
    error('seizmo:writekernels:badInput',...
        'Y step size is not uniform or is <=0!');
end

% directory separator
fs=filesep;

% handle file
if(isempty(file))
    % graphical selection
    [file,path]=uiputfile(...
        {'*.kernel;*.KERNEL;kernel.*;KERNEL.*;' ...
        'Kernel Files (*.kernel,KERNEL.*)';
        '*.kern;*.KERN;kern.*;KERN.*;' ...
        'Kern Files (*.kern,KERN.*)';
        '*.dat;*.DAT' 'DAT Files (*.dat,*.DAT)';
        '*.*' 'All Files (*.*)'},...
        'Select Kernel File for Writing');
    if(isequal(0,file))
        error('seizmo:writekernels:noFileSelected',...
            'No input file selected!');
    end
    file=[path fs file];
else
    % check file
    if(~isstring(file))
        error('seizmo:writekernels:badInput',...
            'FILENAME must be a string!');
    end
    if(~isabspath(file)); file=[pwd fs file]; end
    if(exist(file,'dir'))
        error('seizmo:writekernels:dirConflict',...
            'File: %s\nIs A Directory!',file);
    end
    if(exist(file,'file'))
        if(~o)
            fprintf('KERNEL File: %s\nFile Exists!\n',file);
            reply=input('Overwrite? Y/N [N]: ','s');
            if(isempty(reply) || ~strncmpi(reply,'y',1))
                disp('Not overwriting!');
                return;
            end
            disp('Overwriting!');
        end
    end
end

% combine (note the switch)
m=[x(:) y(:) Kph(:) Kam(:)];

% write out header portion
n=size(x,1); d=y(2)-y(1);
dlmwrite(file,[n x(1) d; n x(1) d],' ')

% write kernels
dlmwrite(file,m,'-append','delimiter',' ');

end