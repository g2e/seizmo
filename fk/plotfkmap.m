function [varargout]=plotfkmap(map)
%PLOTFKMAP    Plots the frequency-wavenumber output from FKMAP
%
%    Usage:    h=plotfkmap(map)
%
%    Description: H=PLOTFKMAP(MAP) plots a slowness map using the struct
%     MAP which was output from FKMAP.  See FKMAP for details on the
%     struct.  This is mainly so you can save the results and replot them
%     later (because FKMAP is quite slow).
%
%    Notes:
%     - The circles of the bull's eye in the plot correspond to several
%       surface/diffracted/head seismic phases (which appear depends on the
%       plot's maximum slowness):
%        Phase:       Rg    Lg    Sn    Pn    Sdiff  Pdiff  PKPcdiff
%        Vel (km/s):  3.0   4.0   4.5   8.15  13.3   25.1   54.0
%        S (s/deg):   37.0  27.8  24.7  13.6  8.36   4.43   2.06
%       The radial lines correspond to 30deg steps in backazimuth.
%
%    Examples:
%     Show slowness map for a dataset at about 50s periods:
%      map=fkmap(data,50,201,[1/51 1/49]);
%      plotfkmap(map);
%
%    See also: FKMAP, FKARF

%     Version History:
%        May   4, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   4, 2010 at 11:00 GMT

% todo:
% - increase checks for new fields
% - plotting for rose style
% - verify pcolor works for polar by using very pixelated versions

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check map is proper struct
fields={'map' 'nsta' 'stla' 'stlo' 'stel' 'stdp' 'butc' 'eutc' ...
    'npts' 'delta' 'smax' 'spts' 'freqs' 'rose' 'center' 'normdb'};
if(~isstruct(map) || ~isscalar(map) ...
        || ~all(ismember(fields,fieldnames(map))))
    error('seizmo:plotfkmap:badInput',...
        ['MAP must be a struct with the fields:\n' ...
        sprintf('%s ',fields{:})]);
end

% consistency check
if(~isreal(map.map) || ~isequal(size(map.map),[map.spts map.spts]) ...
        || ~isreal(map.nsta) || ~isscalar(map.nsta) ...
        || map.nsta~=fix(map.nsta) || map.nsta<2 ...
        || ~isequal(size(map.stla),[map.nsta 1]) || ~isreal(map.stla) ...
        || ~isequal(size(map.stlo),[map.nsta 1]) || ~isreal(map.stlo) ...
        || ~isequal(size(map.stel),[map.nsta 1]) || ~isreal(map.stel) ...
        || ~isequal(size(map.stdp),[map.nsta 1]) || ~isreal(map.stdp) ...
        || ~isequal(size(map.butc),[1 5]) || ~isreal(map.butc) ...
        || ~isequal(size(map.eutc),[1 5]) || ~isreal(map.eutc) ...
        || ~isreal(map.npts) || ~isscalar(map.npts) ...
        || map.npts~=fix(map.npts) || map.npts<0 ...
        || ~isreal(map.delta) || ~isscalar(map.delta) || map.delta<=0 ...
        || ~isreal(map.smax) || ~isscalar(map.smax) || map.smax<=0 ...
        || ~isreal(map.spts) || ~isscalar(map.spts) ...
        || map.spts~=fix(map.spts) || map.spts<0 ...
        || ~isreal(map.freqs) || numel(map.freqs)<0 ...
        || ~isvector(map.freqs) || any(map.freqs)<=0)
    warning('seizmo:plotfkmap:badStruct',...
        'MAP appears to be corrupt!');
end

% get pertinent info
rose=map.rose;
smax=map.smax;
spts=size(map.map,1);
sx=-smax:2*smax/(spts-1):smax;
fmin=min(map.freqs);
fmax=max(map.freqs);

% plotting slowness space
if(rose)
    bazpts=map.bazpts;
    figure('color','default');
    defaultaxescolor=get(0,'defaultaxescolor');
    defaultaxesxcolor=get(0,'defaultaxesxcolor');
    defaultaxesycolor=get(0,'defaultaxesycolor');
    set(0,'defaultaxescolor','w');
    set(0,'defaultaxesxcolor','w')
    set(0,'defaultaxesycolor','w')
    ph=polar([0 2*pi],[0 smax]);
    axis('ij');
    %ph=mmpolar([0 2*pi],[0 smax],'style','compass',...
    %    'backgroundcolor','k','bordercolor','w');
    set(0,'defaultaxescolor',defaultaxescolor);
    set(0,'defaultaxesxcolor',defaultaxesxcolor)
    set(0,'defaultaxesycolor',defaultaxesycolor)
    delete(ph);
    hold on;
    view([-90 90]);
    smag=(0:spts-1)/(spts-1)*smax;
    smag=smag(ones(bazpts,1),:)';
    baz=(0:bazpts-1)/(bazpts-1)*360*pi/180;
    baz=baz(ones(spts,1),:);
    [x,y]=pol2cart(baz,smag);
    pcolor(x,y,map.map);
    colormap(ritz);
    set(gca,'clim',[-12 0]);
    shading flat;
    colorbar;
    die
end

% first plot the map
h=figure;
imagesc(sx,fliplr(sx),map.map);
set(gca,'ydir','normal');
set(gca,'clim',[-12 0]);
hold on

% next plot the bull's eye
% Phase:       Rg    Lg    Sn    Pn    Sdiff  Pdiff  PKPcdiff
% Vel (km/s):  3.0   4.0   4.5   8.15  13.3   25.1   54.0
% S (s/deg):   37.0  27.8  24.7  13.6  8.36   4.43   2.06
if(smax>=37)
    ph=[37 27.8 24.7 13.6 8.36 4.43];
elseif(smax>=28)
    ph=[27.8 24.7 13.6 8.36 4.43];
elseif(smax>=25)
    ph=[24.7 13.6 8.36 4.43];
elseif(smax>=14)
    ph=[13.6 8.36 4.43 2.06];
elseif(smax>=8.5)
    ph=[8.36 4.43 2.06];
elseif(smax>=4.5)
    ph=[4.43 2.06];
else
    ph=2.06;
end
[x,y]=circle(ph(1),12);
[x2,y2]=circle(ph(end),12);
plot([x; x2],[y; y2],'w','linewidth',1,'tag','bullseye');
for i=ph
    [x,y]=circle(i);
    plot(x,y,'w','linewidth',1,'tag','bullseye');
end
hold off

% finally take care of labels/coloring/etc
title({...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Frequency Range:    ' num2str(fmin) 'Hz to ' num2str(fmax) 'Hz']},...
    'fontweight','bold');
xlabel('East/West Slowness (s/deg)','fontweight','bold');
ylabel('North/South Slowness (s/deg)','fontweight','bold');
colormap(ritz);
c=colorbar('eastoutside','fontweight','bold');
set(c,'xaxislocation','top');
set(gca,'fontweight','bold');
xlabel(c,'dB')
axis equal tight
if(nargout); varargout{1}=h; end

end
