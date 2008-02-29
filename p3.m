function [fh]=p3(data,xlimits,ylimits,scale,style,label,fh,sfh) 
%P3   Plot evenly spaced record section
%
%    Description: Plots timeseries and xy records spaced out evenly in the 
%     same plot.  Other record types are ignored.  Optional inputs are the 
%     x/y plot limits ([low high] - y-axis is record number), the amplitude
%     scaling factor (in y-axis units), the style of normalization ('record'
%     or 'global' - scale each record independently or as a group), the 
%     label logical (y-axis ticks have record names for labels - default is
%     0 which is off) and the figure and subplot handles (to put p3 plot in
%     a particular figure and subplot). Default amplitude scaling is 1/10th
%     the number of records.  Default normalization style is 'record'.  
%     Output is the figure handle.
%
%    Usage:  [fh]=p3(data,xlim,ylim,scale,style,labelnames,fh,sfh)
%
%    Notes:  The use of ylim is there but it really isn't useful as you can
%            just select a specific range of data using data(a:b) to do the
%            same thing.
%
%    See also:  p1, p2, recsec

% check number of inputs
error(nargchk(1,8,nargin))

% check data structure
error(seischk(data,'x'))

% plotting style defaults
bgc='k';        % background color
fgc='w';        % axis/text color
ts=4;           % text size
tn='times';     % text name
tw='light';     % text weight
lw=1;           % record line width
%fw=3;           % flag width
cmap='hsv';     % record colormap
ticdir='out';   % tick direction
ticlen=[0 0];   % tick length [2D 3D]
gridit='off';   % grid parameter

% allow access to plot styling using a global structure
global SEISLAB
if(isfield(SEISLAB,'BGCOLOR')); bgc=SEISLAB.BGCOLOR; end
if(isfield(SEISLAB,'FGCOLOR')); fgc=SEISLAB.FGCOLOR; end
if(isfield(SEISLAB,'FONTSIZE')); ts=SEISLAB.FONTSIZE; end
if(isfield(SEISLAB,'FONTNAME')); tn=SEISLAB.FONTNAME; end
if(isfield(SEISLAB,'FONTWEIGHT')); tw=SEISLAB.FONTWEIGHT; end
if(isfield(SEISLAB,'TRACEWIDTH')); lw=SEISLAB.TRACEWIDTH; end
%if(isfield(SEISLAB,'BARWIDTH')); fw=SEISLAB.BARWIDTH; end
if(isfield(SEISLAB,'COLORMAP')); cmap=SEISLAB.COLORMAP; end
if(isfield(SEISLAB,'TICKDIR')); ticdir=SEISLAB.TICKDIR; end
if(isfield(SEISLAB,'TICKLEN')); ticlen=SEISLAB.TICKLEN; end
if(isfield(SEISLAB,'GRIDIT')); gridit=SEISLAB.GRIDIT; end

% initialize plot
if(nargin<7 || isempty(fh) || fh<1); fh=figure;
else figure(fh); if(nargin==8 && ~isempty(sfh)); subplot(sfh); end; end
whitebg(bgc);
set(gcf,'Name','P3 -- SAC Seismogram Plotting Utility', ...
    'NumberTitle','off','color',bgc,'Pointer','crosshair');
set(gca,'FontName',tn,'FontWeight',tw,'FontSize',ts,'box','on',...
    'xcolor',fgc,'ycolor',fgc,'TickDir',ticdir,'ticklength',ticlen);
grid(gridit);

% number of records
nrecs=length(data);

% check options
if(nargin<4 || isempty(scale)); scale=nrecs/10; end
if(nargin<5 || isempty(style)); style='record'; end
if(~any(strcmp(style,{'record' 'global'})))
    warning('seislab:p3:badInput','bad normalization style')
    style='record';
end

% record coloring
cmap=str2func(cmap);
colors=cmap(nrecs);

% header info
leven=glgc(data,'leven');
iftype=genum(data,'iftype');
[b,npts,delta]=gh(data,'b','npts','delta');

% check filetype
goodfiles=(strcmp(iftype,'itime') | strcmp(iftype,'ixy'));
indices=find(goodfiles).';

% find max amplitude for global style normalization
if(strcmp(style,'global'));
    ampmax=0;
    for i=indices
        tmax=max(abs(data(i).x));
        if(tmax>ampmax); ampmax=tmax; end
    end
end

% loop through each file
hold on
for i=indices
    % get record timing
    if(strcmp(leven(i),'true')); time=(b(i)+(0:npts(i)-1)*delta(i)).';
    else time=data(i).t; end
    
    % get max amplitude for record style normalization
    if(strcmp(style,'record')); ampmax=max(abs(data(i).x)); end
    
    % plot series
    plot(time,i+data(i).x*scale/ampmax,'color',colors(i,:),'linewidth',lw);
end
hold off

% zooming
axis tight;
if(nargin>1 && ~isempty(xlimits)); axis auto; xlim(xlimits); end
if(nargin>2 && ~isempty(ylimits)); ylim(ylimits); end

% label with record names
if(nargin>5 && label==1)
    if(isfield(data,'name'))
        set(gca,'ytick',1:nrecs);
        set(gca,'yticklabel',{data.name}.');
    end
end

% title (total plotted / total records)
title([num2str(length(indices)) '/' num2str(nrecs) ' Records']);

end
