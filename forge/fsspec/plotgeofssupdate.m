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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  8, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check struct
error(chkgeofss(s));

% don't allow array/volume
if(~isscalar(s) || any(s.vector))
    error('seizmo:plotgeofssupdate:badInput',...
        'S must be a scalar geofss struct and be averaged w/ GEOFSSAVG!');
end

% just replot if ax isn't an axes handle
if(nargin<2)
    ax=plotgeofss(s);
    if(nargout); varargout{1}=ax; end
    return;
elseif(~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')) ...
        || ~strcmp(get(ax,'createfcn'),'plotgeofss'))
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

% convert to dB
switch s.method
    case 'coarray'
        s.spectra=10*log10(pos(real(s.spectra)));
    otherwise
        s.spectra=10*log10(s.spectra);
end

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
set(get(ax,'Title'),'string',...
    {[] ['Number of Stations:  ' num2str(s.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',s.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',s.eutc) ' UTC'] ...
    ['Period    :  ' num2str(1/fmax) ' to ' num2str(1/fmin) ' s'] ...
    ['Slowness  :  ' num2str(smin) ' to ' num2str(smax) ' s/^o'] ...
    ['0 dB = ' num2str(zdb) 'dB'] []});

if(nargout); varargout{1}=ax; end

end
