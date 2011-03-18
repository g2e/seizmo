function [bad,azlim,ddlim,ax]=delazcut(cpnt,pnts,azlim,ddlim,color,varargin)
%DELAZCUT    Interactively select degree-distance & azimuth range
%
%    Usage:    [bad,azlim,ddlim,ax]=delazcut(cpnt,pnts)
%              [...]=delazcut(cnpt,pnts,azlim,ddlim)
%              [...]=delazcut(cnpt,pnts,azlim,ddlim,color)
%              [...]=delazcut(cnpt,pnts,azlim,ddlim,color,'opt1',val1,...)
%
%    Description:
%     [BAD,AZLIM,DDLIM,AX]=DELAZCUT(CPNT,PNTS) presents a map in azimuthal
%     equidistant projection centered on CPNT with a distance and azimuth
%     grid drawn in blue overlain by the lat/lon positions in PNTS.  The
%     map is interactive and allows the user to select a geographic region
%     in terms of azimuth and distance.  CPNT & PNTS should be formated as
%     [LAT LON].  The output BAD indicates which points in PNTS are _not_
%     in the region selected.  AZLIM & DDLIM are the limits of the azimuth
%     and degree distance region.  AX is the handle of the axes of the map.
%
%     [...]=DELAZCUT(CNPT,PNTS,AZLIM,DDLIM) pre-sets the azimuth & degree
%     distance limits.  The user may alter these (the new values are
%     returned in the output AZLIM & DDLIM).
%
%     [...]=DELAZCUT(CNPT,PNTS,AZLIM,DDLIM,COLOR) sets the face color of
%     the points in PNTS.  This may be a valid NAME2RGB string or an array
%     of rgb triplets as a Nx3 array where N is the number of points in
%     PNTS.
%
%     [...]=DELAZCUT(CNPT,PNTS,AZLIM,DDLIM,COLOR,'OPT1',VAL1,...) passes
%     additional options to MMAP.
%
%    Notes:
%
%    Examples:
%     % This is useful for deleting records in a seizmo dataset based on
%     % their station's position from an event:
%     [ev,st]=getheader(data,'ev','st');
%     data(delazcut(ev(1,1:2),st(:,1:2)))=[];
%
%    See also: ARRCUT, AMPCUT, ERRCUT, SNRCUT, POPCUT

%     Version History:
%        Mar.  9, 2011 - initial version
%        Mar. 18, 2011 - added informative title
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 18, 2011 at 23:55 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% default inputs
if(nargin<3 || isempty(azlim)); azlim=[0 360]; end
if(nargin<4 || isempty(ddlim)); ddlim=[0 179.999]; end % 180 causes trouble
if(nargin<5 || isempty(color)); color=[]; end

% check inputs
if(~isreal(cpnt) || ~isequal(size(cpnt),[1 2]))
    error('seizmo:delazcut:badInput',...
        ['CPNT must be a real-valued vector giving the position of a\n' ...
        'single point as [LAT LON] for getting distance & azimuth!']);
elseif(~isreal(pnts) || size(pnts,2)~=2)
    error('seizmo:delazcut:badInput',...
        'PNTS must be real-valued array points as [LAT LON]!');
elseif(~isreal(azlim) || ~isequal(size(azlim),[1 2]) || diff(azlim)<0)
    error('seizmo:delazcut:badInput',...
        'AZLIM must be the default azimuth limits as [AZMIN AZMAX]!');
elseif(~isreal(ddlim) || ~isequal(size(ddlim),[1 2]) || diff(ddlim)<0 ...
        || any(ddlim<0 | ddlim>180))
    error('seizmo:delazcut:badInput',...
        'DDLIM must be the default distance limits as [DDMIN DDMAX]!');
end
if(~isempty(color))
    if(ischar(color))
        % keep 'none' or try name2rgb (it errors if not valid)
        if(~strcmpi(color,'none')); color=name2rgb(color); end
    elseif(~isreal(color) || ndims(color)~=2 ...
            || any(color(:)<0 | color(:)>1) || size(color,2)~=3 ...
            || ~any(size(color,1)~=[1 size(pnts,1)]))
        error('seizmo:delazcut:badInput',...
            'Numeric COLOR must be a valid rgb triplet!');
    end
end

% distance/azimuth
[dd,az]=sphericalinv(cpnt(1),cpnt(2),pnts(:,1),pnts(:,2));

% create map
ax=mmap('proj','Azimuthal Equidistant',...
    'po',{'lat',cpnt(1),'lon',cpnt(2),'rad',min(180,max(dd)*1.1)},...
    'st',pnts,'stms',15,'ev',cpnt,varargin{:});
mapeventgrid(ax,cpnt(1),cpnt(2));
if(~isempty(color) && isnumeric(color))
    set(findobj(ax,'tag','stations'),'cdata',color);
elseif(~isempty(color) && ischar(color))
    if(strcmpi(color,'none'))
        set(findobj(ax,'tag','stations'),'markerfacecolor',color);
    else
        set(findobj(ax,'tag','stations'),'cdata',name2rgb(color));
    end
end
drawnow;
movekids(findobj(ax,'tag','stations'),'front');
movekids(findobj(ax,'tag','events'),'front');

% retrieve coloring
fgcolor=get(findobj(ax,'tag','m_grid_box'),'color');

% instructive title
title(ax,{'LEFT CLICK:   Set Distance/Azimuth Limit' ...
    'MIDDLE CLICK:   Finalize Distance/Azimuth Limits' ...
    'RIGHT CLICK:   Undo last LEFT CLICK'},...
    'color',fgcolor);

% coloring
bgcolor=invertcolor(fgcolor,true);

% draw/get initial selection
if(~isequal(azlim,[0 360]) || ~isequal(ddlim,[0 179.999]))
    % fix azlim >360
    if(diff(azlim)>360); azlim(2)=azlim(1)+360; end
    
    % set variables appropriately
    none=false; done=true;
    ddrng=ddlim; azrng=azlim;
    
    % draw box
    hold(ax,'on');
    
    % draw the corner
    [lat,lon]=sphericalfwd(cpnt(1),cpnt(2),ddlim(1),azlim(1));
    [x,y]=m_ll2xy(lon,lat);
    pnth=plot(ax,x,y,'.','color',bgcolor,'markersize',20,'linewidth',2);
    
    % draw radial box lines
    [lat,lon]=sphericalfwd(cpnt(1),cpnt(2),...
        ddlim([1 1],:)',azlim([1 1],:));
    [x,y]=m_ll2xy(lon,lat);
    boxh(1:2)=plot(ax,x,y,'color',bgcolor,'linewidth',2);
    
    % draw arc box lines
    [lat,lon]=sphericalfwd(cpnt(1),cpnt(2),...
        ddlim([ones(30,1) 2*ones(30,1)]),...
        [linspace(azlim(1),azlim(2),30)' ...
        linspace(azlim(1),azlim(2),30)']);
    [x,y]=m_ll2xy(lon,lat);
    boxh(3:4)=plot(ax,x,y,'color',bgcolor,'linewidth',2);
    
    % finished box
    hold(ax,'off');
    
    % get bad/good
    aztmp=az;
    while(any(aztmp-azlim(2)>0))
        aztmp(aztmp>azlim(2))=aztmp(aztmp>azlim(2))-360;
    end
    while(any(aztmp-azlim(1)<0))
        aztmp(aztmp<azlim(1))=aztmp(aztmp<azlim(1))+360;
    end
    bad=(aztmp<azlim(1) | aztmp>azlim(2) | dd<ddlim(1) | dd>ddlim(2));
else
    % default to all good
    none=true; done=false;
    bad=false(size(pnts,1),1);
end

% make selection
while(true)
    % get user input
    try
        % force axes to front first
        axes(ax);
        [x,y,button]=ginput(1);
    catch
        % plot closed so finish
        return;
    end
    
    % require ax
    if(gca~=ax); continue; end
    
    % act by button
    switch button
        case 1 % left or ctrl+shift
            % 1st or 2nd corner?
            if(none) % 1st = show point
                % x/y 2 lat/lon
                [lon,lat]=m_xy2ll(x,y);
                
                % lat/lon 2 dd & az
                [ddrng(1),azrng(1)]=sphericalinv(cpnt(1),cpnt(2),lat,lon);
                
                % draw + on image
                hold(ax,'on');
                pnth=plot(ax,x,y,'.',...
                    'color',bgcolor,'markersize',20,'linewidth',2);
                hold(ax,'off');
                
                % increment
                none=false;
            else % 2nd = draw box
                % x/y 2 lat/lon
                [lon,lat]=m_xy2ll(x,y);
                
                % lat/lon 2 dd & az
                [ddrng(2),azrng(2)]=sphericalinv(cpnt(1),cpnt(2),lat,lon);
                
                % update output limits
                ddlim=ddrng; azlim=azrng;
                
                % force dd/az to increase
                if(diff(ddlim)<0); ddlim=fliplr(ddlim); end
                if(diff(azlim)<0); azlim(2)=azlim(2)+360; end
                
                % drawing new box using 1st pnt + current so delete last
                if(done); delete(boxh); end
                
                % draw box
                hold(ax,'on');
                
                % draw radial box lines
                [lat,lon]=sphericalfwd(cpnt(1),cpnt(2),...
                    ddlim([1 1],:)',azlim([1 1],:));
                [x,y]=m_ll2xy(lon,lat);
                boxh(1:2)=plot(ax,x,y,'color',bgcolor,'linewidth',2);
                
                % draw arc box lines
                [lat,lon]=sphericalfwd(cpnt(1),cpnt(2),...
                    ddlim([ones(30,1) 2*ones(30,1)]),...
                    [linspace(azlim(1),azlim(2),30)' ...
                    linspace(azlim(1),azlim(2),30)']);
                [x,y]=m_ll2xy(lon,lat);
                boxh(3:4)=plot(ax,x,y,'color',bgcolor,'linewidth',2);
                
                % finished box
                hold(ax,'off');
                
                % 2 points
                done=true;
                
                % get points inside/outside box
                aztmp=az;
                while(any(aztmp-azlim(2)>0))
                    aztmp(aztmp>azlim(2))=aztmp(aztmp>azlim(2))-360;
                end
                while(any(aztmp-azlim(1)<0))
                    aztmp(aztmp<azlim(1))=aztmp(aztmp<azlim(1))+360;
                end
                bad=(aztmp<azlim(1) | aztmp>azlim(2) ...
                    | dd<ddlim(1) | dd>ddlim(2));
            end
        case 2 % middle or shift
            return;
        case 3 % right or ctrl
            % undo last action
            if(none) % skip if nothing done
                continue;
            elseif(done) % delete bounding box
                % delete box
                delete(boxh);
                
                % back to all good
                ddlim=[0 180];
                azlim=[0 360];
                bad(:)=false;
                
                % there is only 1 pnt
                done=false;
            else % delete point & decrement
                delete(pnth);
                none=true;
            end
        otherwise
            key2zoompan(button,ax);
    end
end

end
