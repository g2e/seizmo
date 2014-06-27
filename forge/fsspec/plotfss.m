function [varargout]=plotfss(s,varargin)
%PLOTFSS    Plots frequency-slowness spectra
%
%    Usage:    plotfss(s)
%              plotfss(s,dblim)
%              plotfss(s,dblim,zerodb)
%              plotfss(s,dblim,zerodb,fgcolor,bgcolor)
%              plotfss(s,dblim,zerodb,fgcolor,bgcolor,ax)
%              ax=plotfss(...)
%
%    Description:
%     PLOTFSS(S) plots the frequency-slowness spectra S as power spectral
%     density scaled so that the maximum value corresponds to 0dB and any
%     spectra within 12dB are given different colors in the FIRE colormap.
%     See FSS for details on the struct S.  Note that the FSS spectra must
%     be scalar along the frequency dimension (use FSSAVG).  Note that this
%     function requires the spectra to be sampled in a regular grid which
%     should not be a problem for spectra from FSS.
%
%     PLOTFSS(S,DBLIM) sets the dB limits for coloring the spectra.
%     The default is [-12 0] for the default ZERODB (see next Usage form).
%     If ZERODB IS 'min' or 'median', the default DBLIM is [0 12].  DBLIM
%     must be a real-valued 2-element vector.
%
%     PLOTFSS(S,DBLIM,ZERODB) changes what 0dB corresponds to in the
%     plot.  The allowed values are 'min', 'max', 'median', & 'abs'.  The
%     default is 'max'.
%
%     PLOTFSS(S,DBLIM,ZERODB,FGCOLOR,BGCOLOR) specifies foreground and
%     background colors of the plot.  The default is 'w' for FGCOLOR & 'k'
%     for BGCOLOR.  Note that if one is specified and the other is not, an
%     opposing color is found using INVERTCOLOR.  For aesthetics, the
%     colormap is also changed so the noise clip gives BGCOLOR.
%
%     PLOTFSS(S,DBLIM,ZERODB,FGCOLOR,BGCOLOR,AX) sets the axes to draw
%     in to the handle AX.  This is useful for subplots, guis, etc.
%
%     AX=PLOTFSS(...) returns the axis handle of the plot.
%
%    Notes:
%     - Polar plots use the function WEDGE for drawing.  See that function
%       to alter the plots beyond PLOTFSS.
%
%    Examples:
%     % Show spectra for a dataset at about 50s periods:
%     plotfss(fssavg(fss(data,50,101,[1/51 1/49])));
%
%    See also: PLOTFSSARF, PLOTFSSUPDATE, FSSFREQSLIDE, FSSFRAMESLIDE,
%              FSS, FSSXC, FSSHORZ, FSSHORZXC

%     Version History:
%        May   4, 2010 - initial version
%        May  11, 2010 - updated for struct changes, got polar working,
%                        coloring options, option for setting axes handle
%        May  21, 2010 - display period rather than frequency
%        May  24, 2010 - labeling the top of colorbar is broken in r2009a
%        May  26, 2010 - added dblim & zerodb args, updated docs
%        June 16, 2010 - labels call correct axes, update see also section
%        July  6, 2010 - major update for new struct
%        Oct. 10, 2010 - all plotting functions use proper ax calls, tagged
%                        plots as 'fks'
%        Jan.  7, 2011 - delete commented mmpolar lines
%        Feb. 16, 2011 - fix pcolor offset in polar plots, color code fix
%        Feb. 23, 2011 - replace polar call with wedge, fix pcolor
%                        out-of-bounds pixels
%        Apr.  4, 2012 - minor doc update
%        Sep. 12, 2012 - adapted from plotfkmap
%        Sep. 27, 2012 - spectra to db simplified
%        Mar. 25, 2014 - scale unwhitened spectra appropriately
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2014 at 10:50 GMT

% todo:
% - pv passthrough
% - remove redundant code

% check nargin
error(nargchk(1,6,nargin));

% check fss struct
error(chkfss(s));

% require scalar struct
if(~isscalar(s))
    error('seizmo:plotfss:badInput',...
        'S must be a scalar fss struct');
end

% require spectra to be averaged or scalar in frequency domain
if(size(s.spectra,3)~=1)
    error('seizmo:plotfss:badInput',...
        'S needs to be reduced using FSSAVG!');
end

% plotting function call depends on polar
if(s.polar)
    [ax]=plotfsspolar(s,varargin{:});
else % cartesian
    [ax]=plotfsscart(s,varargin{:});
end
if(nargout); varargout{1}=ax; end

end



function ax=plotfsspolar(s,dblim,zerodb,fgcolor,bgcolor,ax)

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'med' 'abs'}))
    error('seizmo:plotfss:badZERODB',...
        'ZERODB must be ''min'' ''max'' ''median'' or ''abs''!');
end
zerodb=lower(zerodb);

% default/check dblim
if(nargin<2 || isempty(dblim))
    switch zerodb
        case {'min' 'median' 'med'}
            dblim=[0 12];
        case {'max' 'abs'}
            dblim=[-12 0];
    end
end
if(~isreal(dblim) || numel(dblim)~=2)
    error('seizmo:plotfss:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% scale factor if unwhitened
if(~s.whiten); s.spectra=2*s.delta*s.spectra/s.npts; end

% convert to dB
s.spectra=10*log10(abs(s.spectra));

% rescale spectra
switch zerodb
    case 'min'
        zdb=min(s.spectra(:));
        s.spectra=s.spectra-zdb;
    case 'max'
        zdb=max(s.spectra(:));
        s.spectra=s.spectra-zdb;
    case {'median' 'med'}
        zdb=nanmedian(s.spectra(:));
        s.spectra=s.spectra-zdb;
    case 'abs'
        zdb=0;
end

% check colors
if(nargin<4); fgcolor=[]; end
if(nargin<5); bgcolor=[]; end
if(nargin<6); ax=[]; end
if(isempty(fgcolor))
    if(isempty(bgcolor))
        if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
                || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
            fgcolor='w'; bgcolor='k';
        else
            bgcolor=get(get(ax,'parent'),'color');
            fgcolor=invertcolor(bgcolor,true);
        end
    else
        fgcolor=invertcolor(bgcolor,true);
    end
elseif(isempty(bgcolor))
    bgcolor=invertcolor(fgcolor,true);
end

% change char to something rgb
if(ischar(fgcolor)); fgcolor=name2rgb(fgcolor); end
if(ischar(bgcolor)); bgcolor=name2rgb(bgcolor); end

% check handle
if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')))
    fh=figure('color',bgcolor);
    ax=axes('parent',fh);
else
    axes(ax);
end

% plot polar image
wedge(ax,s.x,s.y,double(s.spectra),...
    'rcolor',fgcolor,'azcolor',fgcolor,...
    'backgroundcolor',bgcolor,'gridcolor',fgcolor);

% add title
s.butc=[s.butc(1:4) fix(s.butc(5)) fix(1000*mod(s.butc(5),1))];
s.eutc=[s.eutc(1:4) fix(s.eutc(5)) fix(1000*mod(s.eutc(5),1))];
title(ax,{['Stations: ' num2str(s.nsta)] ...
    ['Bgn Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.butc) ' UTC'] ...
    ['End Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.eutc) ' UTC'] ...
    ['Period: ' num2str(1/min(s.freq),'%3.3g') 's to ' ...
    num2str(1/max(s.freq),'%3.3g') 's'] ...
    ['0dB = ' num2str(zdb,'%3.3g') 'dB']},'color',fgcolor);

% modify
set(ax,'clim',dblim);
shading(ax,'flat');
if(strcmp(bgcolor,'w') || isequal(bgcolor,[1 1 1]))
    colormap(ax,flipud(fire));
elseif(strcmp(bgcolor,'k') || isequal(bgcolor,[0 0 0]))
    colormap(ax,fire);
else
    if(ischar(bgcolor))
        bgcolor=name2rgb(bgcolor);
    end
    hsv=rgb2hsv(bgcolor);
    colormap(ax,hsvcustom(hsv));
end
c=colorbar('eastoutside','peer',ax,...
    'xcolor',fgcolor,'ycolor',fgcolor);
%set(c,'xaxislocation','top');
xlabel(c,'dB','color',fgcolor);

% set zerodb in userdata
% - this is for plotfssupdate
userdata.zerodb=zerodb;
userdata.fgcolor=fgcolor;
userdata.bgcolor=bgcolor;
set(ax,'userdata',userdata);
setappdata(ax,'made_by_plotfss',true);

end



function ax=plotfsscart(s,dblim,zerodb,fgcolor,bgcolor,ax)

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'med' 'abs'}))
    error('seizmo:plotfss:badZERODB',...
        'ZERODB must be ''min'' ''max'' ''median'' or ''abs''!');
end
zerodb=lower(zerodb);

% default/check dblim
if(nargin<2 || isempty(dblim))
    switch zerodb
        case {'min' 'median' 'med'}
            dblim=[0 12];
        case {'max' 'abs'}
            dblim=[-12 0];
    end
end
if(~isreal(dblim) || numel(dblim)~=2)
    error('seizmo:plotfss:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% scale factor if unwhitened
if(~s.whiten); s.spectra=2*s.delta*s.spectra/s.npts; end

% convert to dB
s.spectra=10*log10(abs(s.spectra));

% rescale spectra
switch zerodb
    case 'min'
        zdb=min(s.spectra(:));
        s.spectra=s.spectra-zdb;
    case 'max'
        zdb=max(s.spectra(:));
        s.spectra=s.spectra-zdb;
    case {'median' 'med'}
        zdb=nanmedian(s.spectra(:));
        s.spectra=s.spectra-zdb;
    case 'abs'
        zdb=0;
end

% check colors
if(nargin<4); fgcolor=[]; end
if(nargin<5); bgcolor=[]; end
if(nargin<6); ax=[]; end
if(isempty(fgcolor))
    if(isempty(bgcolor))
        if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
                || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
            fgcolor='w'; bgcolor='k';
        else
            bgcolor=get(get(ax,'parent'),'color');
            fgcolor=invertcolor(bgcolor,true);
        end
    else
        fgcolor=invertcolor(bgcolor,true);
    end
elseif(isempty(bgcolor))
    bgcolor=invertcolor(fgcolor,true);
end

% change char to something rgb
if(ischar(fgcolor)); fgcolor=name2rgb(fgcolor); end
if(ischar(bgcolor)); bgcolor=name2rgb(bgcolor); end

% check handle
if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')))
    fh=figure('color',bgcolor);
    ax=axes('parent',fh);
else
    axes(ax);
end

% first plot the s
imagesc(s.x,s.y,s.spectra,'parent',ax);
set(ax,'xcolor',fgcolor,'ycolor',fgcolor,'ydir','normal',...
    'color',bgcolor,'clim',dblim);
hold(ax,'on');

% phase specific bullseye
% Phase:       Rg    Lg    Sn    Pn    Sdiff  Pdiff  PKPcdiff
% Vel (km/s):  3.0   4.0   4.5   8.15  13.3   25.1   54.0
% S (s/deg):   37.0  27.8  24.7  13.6  8.36   4.43   2.06
%if(smax>=37)
%    ph=[37 27.8 24.7 13.6 8.36 4.43];
%elseif(smax>=28)
%    ph=[27.8 24.7 13.6 8.36 4.43];
%elseif(smax>=25)
%    ph=[24.7 13.6 8.36 4.43];
%elseif(smax>=14)
%    ph=[13.6 8.36 4.43 2.06];
%elseif(smax>=8.5)
%    ph=[8.36 4.43 2.06];
%elseif(smax>=4.5)
%    ph=[4.43 2.06];
%else
%    ph=2.06;
%end

% pertinent info for making grid below
smax=max(max(abs(s.x)),max(abs(s.y)));

% regular rings (want only 3 to 5)
pot=[0.1 0.2 0.25 0.5 1 2 2.5 5 10 20 2.5 50 100];
rings=ceil(smax./pot);
idx=find(rings>=3 & rings<=5,1);
ph=(1:rings(idx))*pot(idx);

% plot the bull's eye
% first the radial lines
[x,y]=circle(0,12);
[x2,y2]=circle(ph(end),12);
plot(ax,[x; x2],[y; y2],'color',fgcolor,...
    'linewidth',1,'linestyle',':','tag','bullseye');
% second are the rings
for i=ph
    [x,y]=circle(i);
    plot(ax,x,y,'color',fgcolor,'linewidth',1,'linestyle',':',...
        'tag','bullseye');
end
hold(ax,'off');

% finally take care of labels/coloring/etc
s.butc=[s.butc(1:4) fix(s.butc(5)) fix(1000*mod(s.butc(5),1))];
s.eutc=[s.eutc(1:4) fix(s.eutc(5)) fix(1000*mod(s.eutc(5),1))];
title(ax,{['Stations: ' num2str(s.nsta)] ...
    ['Bgn Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.butc) ' UTC'] ...
    ['End Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.eutc) ' UTC'] ...
    ['Period: ' num2str(1/min(s.freq),'%3.3g') 's to ' ...
    num2str(1/max(s.freq),'%3.3g') 's'] ...
    ['0dB = ' num2str(zdb,'%3.3g') 'dB']},'color',fgcolor);
xlabel(ax,'East Slowness (sec/deg)','color',fgcolor);
ylabel(ax,'North Slowness (sec/deg)','color',fgcolor);
if(strcmp(bgcolor,'w') || isequal(bgcolor,[1 1 1]))
    colormap(ax,flipud(fire));
elseif(strcmp(bgcolor,'k') || isequal(bgcolor,[0 0 0]))
    colormap(ax,fire);
else
    if(ischar(bgcolor))
        bgcolor=name2rgb(bgcolor);
    end
    hsv=rgb2hsv(bgcolor);
    colormap(ax,hsvcustom(hsv));
end
c=colorbar('eastoutside','peer',ax,...
    'xcolor',fgcolor,'ycolor',fgcolor);
%set(c,'xaxislocation','top');
xlabel(c,'dB','color',fgcolor);
axis(ax,'equal','tight');

% set zerodb in userdata
% - this is for plotfssupdate
userdata.zerodb=zerodb;
userdata.fgcolor=fgcolor;
userdata.bgcolor=bgcolor;
set(ax,'userdata',userdata);
setappdata(ax,'made_by_plotfss',true);

end
