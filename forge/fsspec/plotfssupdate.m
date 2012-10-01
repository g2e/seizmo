function [varargout]=plotfssupdate(s,ax)
%PLOTFSSUPDATE    Updates an existing fss plot with a new spectra
%
%    Usage:    plotfssupdate(s,ax)
%              ax=plotfssupdate(s)
%
%    Description:
%     PLOTFSSUPDATE(S,AX) draws an new fss spectra given by S in an
%     existing axes AX created by PLOTFSS.  This is mainly intended for
%     exploring FSS datasets by making movies in a faster fashion.
%
%     AX=PLOTFSSUPDATE(S) is the same as calling PLOTFSS(S) -- ie. a new
%     figure is drawn.
%
%    Notes:
%
%    Examples:
%     % Make & slide through a few fss spectra:
%     s(1)=fssavg(fss(data,25,201,[1/50 1/45],true));
%     s(2)=fssavg(fss(data,25,201,[1/45 1/40],true));
%     s(3)=fssavg(fss(data,25,201,[1/40 1/35],true));
%     s(4)=fssavg(fss(data,25,201,[1/35 1/30],true));
%     s(5)=fssavg(fss(data,25,201,[1/30 1/25],true));
%     s(6)=fssavg(fss(data,25,201,[1/25 1/20],true));
%     ax=plotfss(s(1));
%     for i=2:6
%         pause(1);
%         plotfssupdate(s(i),ax);
%     end
%
%    See also: PLOTFSS, FSSFREQSLIDE, FSSFRAMESLIDE, FSS, FSSXC, FSSHORZ,
%              FSSHORZXC, FSSARF

%     Version History:
%        May  11, 2010 - initial version
%        May  21, 2010 - display period rather than frequency
%        May  26, 2010 - updated for new plotfkmap args (requires passing
%                        info through userdata)
%        June 16, 2010 - better see also section
%        July  6, 2010 - major update for new struct
%        Apr. 13, 2011 - better axis handle usage
%        Apr.  4, 2012 - minor doc update
%        Sep. 13, 2012 - adapted from updatefkmap
%        Sep. 27, 2012 - spectra to db simplified
%        Sep. 29, 2012 - support whiten field, minor title change
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 29, 2012 at 16:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check fss struct
error(chkfss(s));

% require scalar struct
if(~isscalar(s))
    error('seizmo:plotfssupdate:badInput',...
        'S must be a scalar fss struct');
end

% require spectra to be averaged or scalar in frequency domain
if(size(s.spectra,3)~=1)
    error('seizmo:plotfssupdate:badInput',...
        'S needs to be reduced using FSSAVG!');
end

% just replot if ax isn't an axes handle
if(nargin<2)
    ax=plotfss(s);
    if(nargout); varargout{1}=ax; end
    return;
elseif(~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')) ...
        || ~isappdata(ax,'made_by_plotfss'))
    error('seizmo:plotfssupdate:badInput',...
        'AX is not a valid PLOTFSS axes!');
end

% plotting function call depends on polar
if(s.polar)
    plotfssupdatepolar(s,ax);
else % cartesian
    plotfssupdatecart(s,ax);
end
if(nargout); varargout{1}=ax; end

end


function plotfssupdatepolar(s,ax)

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

% convert spectra to dB
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

axes(ax);
h=findobj(ax,'type','surface');
delete(h);
hold(ax,'on');
wedge(ax,s.x,s.y,double(s.spectra));
hold(ax,'off');
s.butc=[s.butc(1:4) fix(s.butc(5)) fix(1000*mod(s.butc(5),1))];
s.eutc=[s.eutc(1:4) fix(s.eutc(5)) fix(1000*mod(s.eutc(5),1))];
set(get(ax,'Title'),'string',...
    {['Stations: ' num2str(s.nsta)] ...
    ['Bgn Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.butc) ' UTC'] ...
    ['End Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.eutc) ' UTC'] ...
    ['Period: ' num2str(1/min(s.freq),'%3.3g') 's to ' ...
    num2str(1/max(s.freq),'%3.3g') 's'] ...
    ['0dB = ' num2str(zdb,'%3.3g') 'dB']});

end


function plotfssupdatecart(s,ax)

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

% convert spectra to dB
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

axes(ax);
h=findobj(ax,'type','image');
delete(h);
hold(ax,'on');
imagesc(s.x,s.y,s.spectra,'parent',ax);
childs=get(ax,'children');
set(ax,'children',[childs(2:end); childs(1)]);
hold(ax,'off');
s.butc=[s.butc(1:4) fix(s.butc(5)) fix(1000*mod(s.butc(5),1))];
s.eutc=[s.eutc(1:4) fix(s.eutc(5)) fix(1000*mod(s.eutc(5),1))];
set(get(ax,'Title'),'string',...
    {['Stations: ' num2str(s.nsta)] ...
    ['Bgn Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.butc) ' UTC'] ...
    ['End Time: ' sprintf('%d.%03d %02d:%02d:%02d.%03d',s.eutc) ' UTC'] ...
    ['Period: ' num2str(1/min(s.freq),'%3.3g') 's to ' ...
    num2str(1/max(s.freq),'%3.3g') 's'] ...
    ['0dB = ' num2str(zdb,'%3.3g') 'dB']});

end
