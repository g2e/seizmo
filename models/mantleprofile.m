function [varargout]=mantleprofile(varargin)
%MANTLEPROFILE    3D mantle model profile (aka radial slice)
%
%    Usage:    mantleprofile('spnt',[lat lon],'epnt',[lat lon])
%              mantleprofile('spnt',[lat lon],'delaz',[arclen az])
%              mantleprofile('cpnt',[lat lon],'delaz',[arclen az])
%              mantleprofile(...,'model',model,...)
%              mantleprofile(...,'deprng',deplimits,...)
%              mantleprofile(...,'startpoint',[lat lon],...)
%              mantleprofile(...,'endpoint',[lat lon],...)
%              mantleprofile(...,'centerpoint',[lat lon],...)
%              mantleprofile(...,'delaz',[gcarc az],...)
%              mantleprofile(...,'dtick',ticks,...)
%              mantleprofile(...,'dlabel',labels,...)
%              mantleprofile(...,'distinc',increment,...)
%              mantleprofile(...,'gcinc',increment,...)
%              mantleprofile(...,'depinc',increment,...)
%              mantleprofile(...,'distpad',kmpad,...)
%              mantleprofile(...,'gcpad',degpad,...)
%              mantleprofile(...,'dvrng',dvlimits,...)
%              mantleprofile(...,'colormap',cmap,...)
%              mantleprofile(...,'fgcolor',color,...)
%              mantleprofile(...,'bgcolor',color,...)
%              mantleprofile(...,'showcolorbar',logical,...)
%              mantleprofile(...,'wedgeopt',{'opt',val,...},...)
%              mantleprofile(...,'axis',ax,...)
%              ax=mantleprofile(...)
%
%    Description:
%     MANTLEPROFILE requires exactly 2 points or 1 point & distance/azimuth
%     info to create a profile.  More or less inputs will lead to an error
%     message.  Some suitable usage forms:
%      MANTLEPROFILE('STARTPNT',[LAT LON],'ENDPNT',[LAT LON])
%      MANTLEPROFILE('STARTPNT',[LAT LON],'DELAZ',[ARCLEN AZ])
%      MANTLEPROFILE('CENTERPNT',[LAT LON],'DELAZ',[ARCLEN AZ])
%
%     MANTLEPROFILE(...,'MODEL',MODEL,...) indicates the 3D mantle model
%     for this map.  The default is 'S20RTS'.  Call AVAILABLE_3DMODELS for
%     a list.
%
%     MANTLEPROFILE(...,'DEPRNG',DEPLIMITS,...) specifies the range of
%     depth values to sample the mantle model at.  The default is [0 2981].
%     The model is sampled every 25km unless altered by the DEPINC option.
%     Units are in kilometers.
%
%     MANTLEPROFILE(...,'STARTPOINT',[LAT LON],...) gives the starting
%     location of the mantle profile.  LAT & LON must be in degrees.  Must
%     be coupled with either the ENDPOINT, CENTERPOINT, or DELAZ options.
%
%     MANTLEPROFILE(...,'ENDPOINT',[LAT LON],...) gives the ending location
%     of the mantle profile.  LAT & LON must be in degrees.  Must be
%     coupled with either the STARTPOINT, CENTERPOINT, or DELAZ options.
%
%     MANTLEPROFILE(...,'CENTERPOINT',[LAT LON],...) gives the center
%     location of the mantle profile.  LAT & LON must be in degrees.  Must
%     be coupled with either the STARTPOINT, ENDPOINT, or DELAZ options.
%
%     MANTLEPROFILE(...,'DELAZ',[GCARC AZ],...) gives the distance and
%     azimuth to the other end of the mantle profile.  If used with
%     CENTERPOINT, then GCARC & AZ indicate the distance and azimuth to the
%     ending position of the profile.  Both GCARC & AZ are in degrees!
%
%     MANTLEPROFILE(...,'DTICK',TICKS,...) sets the depth tick marks and
%     depth grid.  Labels are automatically generated if not given in
%     DLABEL.
%
%     MANTLEPROFILE(...,'DLABEL',LABELS,...) sets the depth tick labels.
%     Can be numeric, char array, or a cellstr array.
%
%     MANTLEPROFILE(...,'DISTINC',INCREMENT,...) sets the distance sampling
%     increment in km.  Sets the same value as the GCINC option.  Not used
%     by default.
%
%     MANTLEPROFILE(...,'GCINC',INCREMENT,...) sets the distance sampling
%     increment in degrees.  Sets the same value as the DISTINC option.
%     Default is .25 degrees.
%
%     MANTLEPROFILE(...,'DEPINC',INCREMENT,...) sets the depth sampling
%     increment in km.  Default is 25km.
%
%     MANTLEPROFILE(...,'DISTPAD',KMPAD,...) extend the profile by this
%     many kilometers on either end.  Useful for raypath plots.
%
%     MANTLEPROFILE(...,'GCPAD',DEGPAD,...) extend the profile by this
%     many degrees on either end.  Useful for raypath plots.
%
%     MANTLEPROFILE(...,'DVRNG',DVLIMITS,...) sets the coloring limits.
%     Anything outside of this range is set to either the maximum or
%     minimum colors in the color map.  The default is set to the limits of
%     the data.  The units are in % dlnv.
%
%     MANTLEPROFILE(...,'COLORMAP',CMAP,...) alters the colormap used in
%     the profile.  The default colormap is 'seis'.  The colormap may be a
%     Nx3 RGB triplet array or a string that may be evaluated to a Nx3 RGB
%     triplet.
%
%     MANTLEPROFILE(...,'FGCOLOR',COLOR,...) specifies the foreground color
%     of the map.  The default is 'w'.  If BGCOLOR is specified and FGCOLOR
%     is not, then FGCOLOR will be set using INVERTCOLOR.
%
%     MANTLEPROFILE(...,'BGCOLOR',COLOR,...) specifies the background color
%     of the map.  The default is 'k'.  If FGCOLOR is specified and BGCOLOR
%     is not, then BGCOLOR will be set using INVERTCOLOR.
%
%     MANTLEPROFILE(...,'SHOWCOLORBAR',LOGICAL,...) turns on/off the
%     drawing of a colorbar.  The default is TRUE.
%
%     MANTLEPROFILE(...,'WEDGEOPT',{'OPT',VAL,...},...) passes additional
%     options to the WEDGE call.
%
%     MANTLEPROFILE(...,'AXIS',AX,...) sets the axes to draw in.  This is
%     useful for subplots, guis, etc.  The default draws the profile in a
%     new figure.
%
%     AX=MANTLEPROFILE(...) returns the axes handle for the profile.
%
%    Notes:
%
%    Examples:
%     % How can I do a greater circle (full 360)?  Specify a point+delaz:
%     mantleprofile('spnt',[0 0],'delaz',[360 -90]); % equatorial slice
%
%     % Fallaron slab:
%     mantleprofile('spnt',[40 -120],'epnt',[40 -60],'dvrng',[-1 1]);
%
%     % African Super Plume:
%     mantleprofile('spnt',[-60 -20],'epnt',[40 60],...
%                   'dvrng',[-2 2],'mo','tx2007');
%
%     % Hawaii Plume North/South & East/West slices:
%     mantleprofile('cpnt',[20 -155],'delaz',[30 0],...
%                   'dvrng',[-1 1],'mo','pri-p05');
%     mantleprofile('cpnt',[20 -155],'delaz',[30 90],...
%                   'dvrng',[-1 1],'mo','pri-p05');
%
%     % A raypath through the African Super Plume:
%     ax=mantleprofile('spnt',[-60 -20],'epnt',[40 60],...
%                      'dvrng',[-2 2],'mo','tx2007','gcpad',5);
%     hold(ax,'on');
%     p=tauppath('deg',95,'ph','P');
%     wedge(ax,p.path.distance,6371-p.path.depth,...
%           'linewidth',2,'color',[.5 .5 .5]);
%
%    See also: MANTLEMAP, MANTLEDV, AVAILABLE_3DMODELS, WEDGE

%     Version History:
%        Feb. 25, 2011 - initial version
%        May   2, 2012 - minor example fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   2, 2011 at 20:30 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('seizmo:mantleprofile:badNumInputs',...
        'Unpaired Option/Value!');
end

% option defaults
varargin=[{'m' 'S20RTS' 'deprng' [0 2891] 'gcinc' .25 'depinc' 25 ...
    'dvrng' [] 'cmap' 'seis' 'wo' [] 'fg' [] 'bg' [] 'cb' true ...
    'tt' true 'a' [] 'gcpad' 0 'dtick' [] 'dlabel' []} varargin];

% check options are strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:mantleprofile:badOption',...
        'All Options must be specified with a string!');
end

% check option/value pairs
cnt=false(4,1);
for i=1:2:numel(varargin)
    % skip empty by default (but still checking option exists)
    skip=false;
    if(isempty(varargin{i+1})); skip=true; end
    val=varargin{i+1};
    
    % check option is available
    switch lower(varargin{i})
        case {'model' 'mo' 'm'}
            if(skip); continue; end
            if(ischar(val))
                model=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'MODEL must be a string!');
            end
        case {'depthrange' 'deprng' 'dr'}
            if(skip); continue; end
            if(isreal(val) && isequal(size(val),[1 2]))
                deprng=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'DEPRNG must be [depth_lo depth_hi]!');
            end
        case {'startpoint' 'start' 'spnt'}
            if(skip)
                cnt(1)=false;
                spnt=[];
            elseif(isreal(val) && isequal(size(val),[1 2]))
                cnt(1)=true;
                spnt=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'SPNT must be [LAT LON]!');
            end
        case {'endpoint' 'end' 'epnt'}
            if(skip)
                cnt(2)=false;
                epnt=[];
            elseif(isreal(val) && isequal(size(val),[1 2]))
                cnt(2)=true;
                epnt=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'EPNT must be [LAT LON]!');
            end
        case {'centerpoint' 'center' 'cpnt'}
            if(skip)
                cnt(3)=false;
                cpnt=[];
            elseif(isreal(val) && isequal(size(val),[1 2]))
                cnt(3)=true;
                cpnt=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'CPNT must be [LAT LON]!');
            end
        case {'delaz'}
            if(skip)
                cnt(4)=false;
                delaz=[];
            elseif(isreal(val) && isequal(size(val),[1 2]))
                cnt(4)=true;
                delaz=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'DELAZ must be [DEGDIST AZIMUTH]!');
            end
        case {'depthtick' 'dtick'}
            if(skip)
                dtick=[];
            elseif(isreal(val) && isvector(val) && all(val>0))
                dtick=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'DTICK must be a vector of positive values in km!');
            end
        case {'depthticklabel' 'dticklabel' 'dlabel'}
            if(skip)
                dlabel={};
            elseif(ischar(val))
                dlabel=cellstr(val);
            elseif(iscellstr(val))
                dlabel=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'DLABEL must be a char/cellstr array of labels!');
            end
        case {'diststep' 'distinc'}
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>0)
                gcstep=val/6371*180/pi;
            else
                error('seizmo:mantleprofile:badInput',...
                    'DISTINC must be a positive scalar in km!');
            end
        case {'gcstep' 'gcinc'}
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>0)
                gcstep=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'GCINC must be a positive scalar in degrees!');
            end
        case {'depthstep' 'depthinc' 'depstep' 'depinc'}
            if(skip); continue; end
            if(isreal(val) && isscalar(val) && val>0)
                depstep=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'DEPTHSTEP must be a positive scalar!');
            end
        case 'distpad'
            if(skip); continue; end
            if(isreal(val) && isscalar(val))
                gcpad=val/6371*180/pi;
            else
                error('seizmo:mantleprofile:badInput',...
                    'DISTPAD must be a scalar in km!');
            end
        case 'gcpad'
            if(skip); continue; end
            if(isreal(val) && isscalar(val))
                gcpad=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'GCPAD must be a scalar in degrees!');
            end
        case {'clim' 'dvrng'}
            if(skip)
                dvrng=[];
            elseif(isreal(val) && isequal(size(val),[1 2]))
                dvrng=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'DVRNG must be [%%dv_lo %%dv_hi]!');
            end
        case {'colormap' 'cmap'}
            if(skip); continue; end
            if(isreal(val) && ndims(val)==2 ...
                    && size(val,2)==3 ...
                    && all(val(:)>=0 & val(:)<=1))
                cmap=val;
            elseif(ischar(val) && isvector(val) && size(val,1)==1)
                cmap=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    ['COLORMAP must be a colormap function\n'...
                    'string or a Nx3 RGB triplet array!']);
            end
        case {'fgcolor' 'fg'}
            if(skip)
                fg=[];
            elseif(ischar(val) ...
                    || (isreal(val) && isequal(size(val),[1 3])))
                fg=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'FGCOLOR must be a colorname or RGB triplet!');
            end
        case {'bgcolor' 'bg'}
            if(skip)
                bg=[];
            elseif(ischar(val) ...
                    || (isreal(val) && isequal(size(val),[1 3])))
                bg=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'BGCOLOR must be a colorname or RGB triplet!');
            end
        case {'showcolorbar' 'colorbar' 'cb'}
            if(skip); continue; end
            if((islogical(val) || isreal(val)) && isscalar(val))
                showcb=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'SHOWCOLORBAR must be TRUE or FALSE!');
            end
        case {'titletype' 'tt'}
            if(skip); continue; end
            if((islogical(val) || isreal(val)) && isscalar(val))
                titletype=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    'TITLETYPE must be a scalar value!');
            end
        case {'wedgeopt' 'wopt' 'wo'}
            if(skip)
                wopt={};
            elseif(iscell(val) && iscellstr(val(1:2:end)))
                wopt=val;
            else
                error('seizmo:mantleprofile:badInput',...
                    ['WEDGEOPT option must be a cell array of ' ...
                    '''property''/value pairs!']);
            end
        case {'axis' 'ax' 'a'}
            if(skip)
                ax=[];
            else
                ax=val;
            end
        otherwise
            error('seizmo:mantleprofile:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% check count
if(sum(cnt)~=2)
    error('seizmo:mantleprofile:badInput',...
        'MANTLEPROFILE requires 2 pnts or 1 pnt + distance & azimuth!');
end

% fix fg/bg colors
if(isempty(fg))
    if(isempty(bg))
        fg='w'; bg='k';
    else
        fg=invertcolor(bg,true);
    end
elseif(isempty(bg))
    bg=invertcolor(fg,true);
end

% convert colornames
if(ischar(fg)); fg=name2rgb(fg); end
if(ischar(bg)); bg=name2rgb(bg); end

% sort out start/end positions + distance & azimuth
if(cnt(1))
    [spnt(1),spnt(2)]=fixlatlon(spnt(1),spnt(2));
    if(cnt(2)) % start + end
        [epnt(1),epnt(2)]=fixlatlon(epnt(1),epnt(2));
        [gcarc,az]=sphericalinv(spnt(1),spnt(2),epnt(1),epnt(2));
    elseif(cnt(3)) % start + center
        [cpnt(1),cpnt(2)]=fixlatlon(cpnt(1),cpnt(2));
        [gcarc,az]=sphericalinv(spnt(1),spnt(2),cpnt(1),cpnt(2));
        gcarc=gcarc*2;
        [epnt(1),epnt(2)]=sphericalfwd(spnt(1),spnt(2),gcarc,az);
    else % start + delaz
        [gcarc,az]=deal(delaz(1),delaz(2));
        [epnt(1),epnt(2)]=sphericalfwd(spnt(1),spnt(2),gcarc,az);
    end
elseif(cnt(2))
    [epnt(1),epnt(2)]=fixlatlon(epnt(1),epnt(2));
    if(cnt(3)) % end + center
        [cpnt(1),cpnt(2)]=fixlatlon(cpnt(1),cpnt(2));
        [gcarc,baz]=sphericalinv(epnt(1),epnt(2),cpnt(1),cpnt(2));
        gcarc=gcarc*2;
        [spnt(1),spnt(2),az]=sphericalfwd(epnt(1),epnt(2),gcarc,baz);
    else % end + delaz
        [gcarc,baz]=deal(delaz(1),delaz(2));
        [spnt(1),spnt(2),az]=sphericalfwd(epnt(1),epnt(2),gcarc,baz);
    end
elseif(cnt(3))
    % center + delaz
    [cpnt(1),cpnt(2)]=fixlatlon(cpnt(1),cpnt(2));
    [gcarc,baz]=deal(delaz(1),delaz(2));
    baz=baz-180;
    [spnt(1),spnt(2),az]=sphericalfwd(cpnt(1),cpnt(2),gcarc,baz);
    gcarc=gcarc*2;
    [epnt(1),epnt(2)]=sphericalfwd(spnt(1),spnt(2),gcarc,az);
end

% now get points along profile at defined step size
gc=-gcpad:gcstep:gcarc+gcpad;
[lat,lon]=sphericalfwd(spnt(1),spnt(2),gc,az);
nll=numel(lat);

% depths
dep=(deprng(1):depstep:deprng(2))';
ndep=numel(dep);

% get %dv
[dv,wtype]=mantledv(model,lat(ones(ndep,1),:),lon(ones(ndep,1),:),...
    dep(:,ones(1,nll)));
dv=dv*100;

% setup axis
if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    % new figure
    fh=figure('color',bg);
    ax=axes('parent',fh,'color',bg);
else
    axes(ax);
    h=get(ax,'children'); delete(h);
    h=findobj(get(get(ax,'parent'),'children'),'peer',ax); delete(h);
end

% automatic depth ticks
if(isempty(dtick)); dtick=deprng; end
if(isempty(dlabel)); dlabel=dtick; end

% plot mantle profile
hold(ax,'on');
wedge(ax,gc,6371-dep,dv,...
    'azlim',[-gcpad gcarc+gcpad],'rlim',6371-fliplr(deprng),...
    'azoffset',-gcarc/2,'rticklabelunits','km',...
    'rtick',6371-flipud(dtick(:)),'rticklabel',flipud(dlabel(:)),...
    'backgroundcolor',bg,'azcolor',fg,'rcolor',fg,'gridcolor',fg,...
    wopt{:});

% pretty it up
colormap(ax,cmap);
if(~isempty(dvrng)); set(ax,'clim',dvrng); end
hold(ax,'off');

% colorbar & title
switch titletype
    case 1
        if(spnt(1)>=0); sns='N'; else sns='S'; end
        if(spnt(2)>=0); sew='E'; else sew='W'; end
        if(epnt(1)>=0); ens='N'; else ens='S'; end
        if(epnt(2)>=0); eew='E'; else eew='W'; end
        title(ax,{['MODEL: ' upper(model)] ...
            ['(' num2str(abs(spnt(1))) '^o' sns ', ' ...
            num2str(abs(spnt(2))) '^o' sew ...
            ')  to  (' num2str(abs(epnt(1))) '^o' ens ', ' ...
            num2str(abs(epnt(2))) '^o' eew ')'] ...
            ['DEGDIST: ' num2str(gcarc) '^o' ...
            ', AZIMUTH: ' num2str(az) '^o'] ''},'color',fg);
end
if(showcb)
    c=colorbar('southoutside','peer',ax,'xcolor',fg,'ycolor',fg);
    xlabel(c,['% \delta{}lnv_' wtype],'color',fg);
end

% export basic info
userdata.model=model;
userdata.spnt=spnt;
userdata.gcarc=gcarc;
userdata.az=az;
set(ax,'userdata',userdata);

% return axes handle
set(ax,'tag','mantleprofile');
if(nargout); varargout{1}=ax; end

end
