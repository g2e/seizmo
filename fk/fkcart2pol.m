function [fk]=fkcart2pol(fk,ass,sss,method)
%FKCART2POL    Converts a cartesian space based fk grid to polar space
%
%    Usage:    fk=fkcart2pol(fk)
%              fk=fkcart2pol(fk,azistepsize,slowstepsize)
%              fk=fkcart2pol(fk,azistepsize,slowstepsize,method)
%
%    Description:
%     FK=FKCART2POL(FK) converts an fk struct regularly sampled in
%     cartesian space (ie. Sx & Sy) to one regularly sampled in polar space
%     (ie. Baz & |S|).  Please note that this requires 2D-based
%     interpolation!
%
%     FK=FKCART2POL(FK,AZISTEPSIZE,SLOWSTEPSIZE) specifies the azimuthal
%     and horizontal slowness step sizes.  The azimuthal step size should
%     be in degrees and the slowness step size in sec/deg.  The defaults
%     are 1 for AZISTEPSIZE while SLOWSTEPSIZE is equal to the slowness
%     sampling interval in the cartesian grid.
%
%     FK=FKCART2POL(FK,AZISTEPSIZE,SLOWSTEPSIZE,METHOD) specifies the 2D
%     interpolation method.  See INTERP2D for possible options.
%
%    Notes:
%
%    Examples:
%     % Plot an fk struct already in cartesian, convert to polar and plot
%     % again to compare:
%     plotfkmap(fk)
%     plotfkmap(fkcart2pol(fk))
%
%    See also: FKMAP, FKVOLUME, FK4D, FKSUBVOL, FKVOL2MAP, PLOTFKMAP

%     Version History:
%        July 14, 2010 - initial version
%        Apr.  4, 2012 - minor doc update
%        Mar.  5, 2014 - fix bugs for array fk input
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  5, 2014 at 16:05 GMT

% todo:
% - significant information is lost
% - no sister function fs_pol2cart
% - step size should be replaced with nsteps

% check nargin
error(nargchk(1,4,nargin));

% check fk struct
error(chkfkstruct(fk));
nfk=numel(fk);

% error if any polar
if(any([fk.polar]))
    error('seizmo:fkcart2pol:badInput',...
        'FK struct is already polar!');
end

% check optional inputs
if(nargin<2 || isempty(ass)); ass=1; end
if(nargin<3 || isempty(sss))
    sss=nan(nfk,1);
    for i=1:nfk
        sss(i)=abs(min(diff(fk(i).x)));
    end
end
if(nargin<4 || isempty(method)); method='linear'; end
if(~isreal(ass) || ~isequalnumelorscalar(fk,ass) || any(ass(:)<=0))
    error('seizmo:fkcart2pol:badInput',...
        'AZISTEPSIZE must be a positive and in degrees!');
end
if(~isreal(sss) || ~isequalnumelorscalar(fk,sss) || any(sss(:)<=0))
    error('seizmo:fkcart2pol:badInput',...
        'SLOWSTEPSIZE must be a positive and in sec/deg!');
end
if(~isstring(method))
    error('seizmo:fkcart2pol:badInput',...
        'METHOD must be a string!');
end

% expand step sizes
if(isscalar(ass)); ass=ass(ones(nfk,1),1); end
if(isscalar(sss)); sss=sss(ones(nfk,1),1); end

% loop over each fk element
for i=1:nfk
    % true x/y grid
    [x,y]=meshgrid(fk(i).x,fk(i).y);
    
    % azi/slow grid
    maxs=max(fk(i).x);
    fk(i).x=0:ass(i):360-ass(i);
    nbaz=numel(fk(i).x);
    fk(i).y=(0:sss(i):maxs)';
    nslow=numel(fk(i).y);
    [baz,slow]=meshgrid(fk(i).x,fk(i).y);
    
    % convert to x/y (east/north)
    [xi,yi]=slowbaz2kxy(slow,baz);
    
    % interpolate
    nfreq=size(fk(i).beam,3);
    zi=nan(nslow,nbaz,nfreq);
    fk(i).beam=fk(i).beam+fk(i).normdb;
    for j=1:nfreq
        zi(:,:,j)=interp2(x,y,fk(i).beam(:,:,j),xi,yi,method);
    end
    
    % assign
    fk(i).beam=zi;
    fk(i).polar=true;
    
    % recompute normdb
    fk(i).normdb=max(fk(i).beam(:));
    fk(i).beam=fk(i).beam-fk(i).normdb;
end

end
