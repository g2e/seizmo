function [varargout]=plotfkazifreq(vol,srng,dblim,zerodb,fgcolor,bgcolor,ax)
%PLOTFKAZIFREQ    Plots beam intensity as a function of azimuth & frequency
%
%    Usage:    plotfkazifreq(vol)
%              plotfkazifreq(vol,srng)
%              plotfkazifreq(vol,srng,dblim)
%              plotfkazifreq(vol,srng,dblim,zerodb)
%              plotfkazifreq(vol,srng,dblim,zerodb,fgcolor,bgcolor)
%              plotfkazifreq(vol,srng,dblim,zerodb,fgcolor,bgcolor,ax)
%              ax=plotfkazifreq(...)
%
%    Description:
%     PLOTFKAZIFREQ(VOL) graphs the maximum beam strength across all
%     slownesses as a function of frequency and azimuth using the struct
%     VOL which is output from FKVOLUME or equivalent.  This is an
%     excellent way to summarize the direction strength with frequency.
%
%     PLOTFKAZIFREQ(VOL,SRNG) limits the horizontal slowness magnitudes
%     when searching for the maximum db value.  The default is [] and is
%     unrestricted.
%
%     PLOTFKAZIFREQ(VOL,SRNG,DBLIM) sets the dB limits for coloring the
%     beam info.  The default is [-12 0] for the default ZERODB (see next
%     Usage form).  If ZERODB IS 'min' or 'median', the default DBLIM is
%     [0 12].  DBLIM must be a real-valued 2-element vector.
%
%     PLOTFKAZIFREQ(VOL,SRNG,DBLIM,ZERODB) changes what 0dB corresponds to
%     in the plot.  The allowed values are 'min', 'max', 'median', & 'abs'.
%     The default is 'max'.
%
%     PLOTFKAZIFREQ(VOL,SRNG,DBLIM,ZERODB,FGCOLOR,BGCOLOR) specifies
%     foreground and background colors of the plot.  The default is 'w' for
%     FGCOLOR & 'k' for BGCOLOR.  Note that if one is specified and the
%     other is not, an opposing color is found using INVERTCOLOR.  The
%     color scale is also changed so the noise clip is at BGCOLOR.
%
%     PLOTFKAZIFREQ(VOL,SRNG,DBLIM,ZERODB,FGCOLOR,BGCOLOR,AX) sets the axes
%     to draw in.  This is useful for subplots, guis, etc.
%
%     AX=PLOTFKAZIFREQ(...) returns the axis handle of the plot.
%
%    Notes:
%
%    Examples:
%     % Limit to surface wave velocities:
%     plotfkazifreq(vol,[25 35]);
%
%    See also: FKVOLUME, FKXCVOLUME, FKFREQSLIDE, FKSUBVOL

%     Version History:
%        July 14, 2010 - initial version
%        Feb. 16, 2011 - color code fix
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 16:05 GMT

% todo:

% check nargin
error(nargchk(1,7,nargin));

% check fk struct
error(chkfkstruct(vol));

% don't allow array/volume
if(~isscalar(vol) || ~vol.volume)
    error('seizmo:plotfkazifreq:badInput',...
        'VOL must be a scalar fk struct and a volume!');
end

% if cartesian map, interpolate to polar
if(~vol.polar)
    vol=fkcart2pol(vol);
end

% check/default slowness range
if(nargin<2 || isempty(srng)); srng=[min(vol.y) max(vol.y)]; end
if(~isreal(srng) || numel(srng)~=2 || any(srng<0))
    error('seizmo:plotfkazifreq:badInput',...
        'SRNG must be [SLOWLO SLOWHI] in s/deg!');
end

% get pertinent slowness info
sidx=vol.y>=srng(1) & vol.y<=srng(2);
smin=min(vol.y(sidx));
smax=max(vol.y(sidx));

% reduce beam space by finding max across slowness
vol.beam=max(vol.beam(sidx,:,:));

% now permute into 2D with baz varying across rows, freq across columns
vol.beam=permute(vol.beam,[2 3 1]);

% default/check scaling type
if(nargin<4 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotfkazifreq:badSTYPE',...
        'STYPE must be ''min'' ''max'' ''median'' or ''abs''!');
end
zerodb=lower(zerodb);

% default/check dblim
if(nargin<3 || isempty(dblim))
    switch zerodb
        case {'min' 'median'}
            dblim=[0 12];
        case {'max' 'abs'}
            dblim=[-12 0];
    end
end
if(~isreal(dblim) || numel(dblim)~=2)
    error('seizmo:plotfkazifreq:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% rescale beam
switch zerodb
    case 'min'
        vol.normdb=vol.normdb+min(vol.beam(:));
        vol.beam=vol.beam-min(vol.beam(:));
    case 'max'
        vol.normdb=vol.normdb+max(vol.beam(:));
        vol.beam=vol.beam-max(vol.beam(:));
    case 'median'
        vol.normdb=vol.normdb+median(vol.beam(:));
        vol.beam=vol.beam-median(vol.beam(:));
    case 'abs'
        vol.beam=vol.beam+vol.normdb;
        vol.normdb=0;
end

% check colors
if(nargin<5);
    fgcolor='w'; bgcolor='k';
elseif(nargin<6)
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
else
    if(isempty(fgcolor))
        if(isempty(bgcolor))
            fgcolor='w'; bgcolor='k';
        else
            fgcolor=invertcolor(bgcolor,true);
        end
    elseif(isempty(bgcolor))
        bgcolor=invertcolor(fgcolor,true);
    end
end

% change char to something rgb
if(ischar(fgcolor)); fgcolor=name2rgb(fgcolor); end
if(ischar(bgcolor)); bgcolor=name2rgb(bgcolor); end

% check handle
if(nargin<7 || isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    figure('color',bgcolor);
    ax=gca;
else
    axes(ax);
end

% plot the map
imagesc(vol.freq,vol.x,vol.beam);
set(ax,'xcolor',fgcolor,'ycolor',fgcolor,'ydir','normal',...
    'color',bgcolor,'fontweight','bold','clim',dblim);

% label
title(ax,{['Number of Stations:  ' num2str(vol.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',vol.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',vol.eutc) ' UTC'] ...
    ['Slowness Range:  ' num2str(smin) ' to ' num2str(smax) 'sec/deg'] ...
    ['0 dB = ' num2str(vol.normdb) 'dB']},...
    'fontweight','bold','color',fgcolor);
xlabel(ax,'Frequency (Hz)','fontweight','bold','color',fgcolor);
ylabel(ax,'Azimuth (deg)','fontweight','bold','color',fgcolor);
if(strcmp(bgcolor,'w') || isequal(bgcolor,[1 1 1]))
    colormap(flipud(fire));
elseif(strcmp(bgcolor,'k') || isequal(bgcolor,[0 0 0]))
    colormap(fire);
else
    if(ischar(bgcolor))
        bgcolor=name2rgb(bgcolor);
    end
    hsv=rgb2hsv(bgcolor);
    colormap(hsvcustom(hsv));
end
c=colorbar('eastoutside','peer',ax,...
    'fontweight','bold','xcolor',fgcolor,'ycolor',fgcolor);
%set(c,'xaxislocation','top');
xlabel(c,'dB','fontweight','bold','color',fgcolor)

% handle output
if(nargout); varargout{1}=ax; end

end
