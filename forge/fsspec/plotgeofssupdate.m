function [varargout]=plotgeofssupdate(s,ax)
%PLOTGEOFSSUPDATE    Updates an existing geofss plot with a new spectra
%
%    Usage:    plotgeofssupdate(s,ax)
%              ax=plotgeofssupdate(s)
%
%    Description:
%     PLOTGEOFSSUPDATE(S,AX) draws the new geofss spectra given by S in an
%     existing axes AX created by PLOTGEOFSS.  This is mainly intended for
%     exploring GEOFSS datasets by making movies in a faster fashion.
%
%     AX=PLOTGEOFSSUPDATE(S) is the same as calling PLOTGEOFSS(S) -- ie. a
%     new figure is drawn.
%
%    Notes:
%
%    Examples:
%     % Create a geofss plot and then change to a new spectra:
%     [lat,lon]=meshgrid(-89:2:89,-179:2:179);
%     s=geofss(d,[lat(:) lon(:)],27:33,[1/50 1/20],'center');
%     s1=geofsssub(s,[1/50 1/40]);
%     s2=geofsssub(s,[1/40 1/30]);
%     ax=plotgeofss(geofssavg(s1));
%     pause(3);
%     plotgeofssupdate(geofssavg(s2),ax);
%
%    See also: PLOTGEOFSS, GEOFSS, GEOFSSXC, GEOFSSSUB, GEOFSSAVG,
%              GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Dec.  8, 2010 - improved degrees symbol usage (^o)
%        Apr.  4, 2012 - minor doc update
%        June  8, 2012 - adapt from updategeofkmap
%        Sep. 29, 2012 - update for struct changes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check struct
error(chkgeofss(s));

% require scalar struct
if(~isscalar(s))
    error('seizmo:plotgeofssupdate:badInput',...
        'S must be a scalar geofss struct');
end

% require scalar freq & slow dimensions
if(size(s.spectra,2)~=1 || size(s.spectra,3)~=1)
    error('seizmo:plotgeofssupdate:badInput',...
        'S needs to be reduced using GEOFSSAVG!');
end

% just replot if ax isn't an axes handle
if(nargin<2)
    ax=plotgeofss(s);
    if(nargout); varargout{1}=ax; end
    return;
elseif(~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')) ...
        || ~isappdata(ax,'made_by_plotgeofss'))
    error('seizmo:plotgeofssupdate:badInput',...
        'AX is not a valid PLOTGEOFSS axes!');
end

% get zerodb/dblim
userdata=get(ax,'userdata');
if(isempty(userdata) || ~isstruct(userdata) ...
        || any(~isfield(userdata,{'zerodb'})))
    zerodb='max';
else
    zerodb=userdata.zerodb;
end

% scale factor if unwhitened
if(~s.whiten); s.spectra=s.spectra/(2*pi); end

% convert to dB
s.spectra=10*log10(abs(s.spectra));

% rescale response
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

% reshape spectra for plotting & to account for pcolor
nlat=numel(unique(s.latlon(:,1)));
nlon=numel(unique(s.latlon(:,2)));
s.latlon=reshape(s.latlon,[nlon nlat 2]); % this breaks if non-meshgrid
s.spectra=reshape(s.spectra,[nlon nlat]); % this breaks if non-meshgrid
latstep=s.latlon(1,2,1)-s.latlon(1,1,1);
lonstep=s.latlon(2,1,2)-s.latlon(1,1,2);
s.latlon(:,:,1)=s.latlon(:,:,1)-latstep/2;
s.latlon(:,:,2)=s.latlon(:,:,2)-lonstep/2;
s.latlon=s.latlon([1:end end],[1:end end],:);
s.latlon(:,end,1)=s.latlon(:,end,1)+latstep;
s.latlon(end,:,2)=s.latlon(end,:,2)+lonstep;
s.spectra=s.spectra([1:end end],[1:end end]);

% convert to map coordinates
[s.latlon(:,:,2),s.latlon(:,:,1)]=m_ll2xy(...
    s.latlon(:,:,2),s.latlon(:,:,1),'clip','patch');

% find previous
pc=findobj(ax,'tag','m_pcolor');

% slip in new data
set(pc(1),...
    'xdata',s.latlon(:,:,2),'ydata',s.latlon(:,:,1),...
    'zdata',0*s.latlon(:,:,2),'cdata',double(s.spectra));

% adjust title
fmin=min(s.freq); fmax=max(s.freq);
smin=min(s.slow); smax=max(s.slow);
s.butc=[s.butc(1:4) fix(s.butc(5)) fix(1000*mod(s.butc(5),1))];
s.eutc=[s.eutc(1:4) fix(s.eutc(5)) fix(1000*mod(s.eutc(5),1))];
set(get(ax,'Title'),'string',...
    {[] ['Stations: ' num2str(s.nsta)] ...
    ['Bgn Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.butc) ' UTC'] ...
    ['End Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.eutc) ' UTC'] ...
    ['Period: ' num2str(1/fmax) ' to ' num2str(1/fmin) ' s'] ...
    ['Slowness: ' num2str(smin) ' to ' num2str(smax) ' s/^o'] ...
    ['0dB = ' num2str(zdb) 'dB'] []});

if(nargout); varargout{1}=ax; end

end
