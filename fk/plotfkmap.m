function []=plotfkmap(map)

% get pertinent info
smax=map.smax;
spts=size(map.map,1);
sx=-smax:2*smax/(spts-1):smax;
fmin=min(map.freqs);
fmax=max(map.freqs);

% plotting slowness space
% Phase:       Rg    Lg    Sn    Pn    Sdiff  Pdiff  PKPcdiff
% Vel (km/s):  3.0   4.0   4.5   8.15  13.3   25.1   54.0
% S (s/deg):   37.0  27.8  24.7  13.6  8.36   4.43   2.06
figure;
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
[x,y]=circle2(ph(1),12);
[x2,y2]=circle2(ph(end),12);
plot([x; x2],[y; y2],'k','linewidth',2);
hold on
for i=ph
    [x,y]=circle(i);
    plot(x,y,'k','linewidth',2);
end
imagesc(sx,fliplr(sx),map.map);
set(gca,'ydir','normal');
set(gca,'clim',[-12 0]);
set(gca,'fontweight','bold');
alpha(0.9)
hold off
title(['Array Slowness Map @ ' ...
    num2str(fmin) 'Hz to ' num2str(fmax) 'Hz'],'fontweight','bold');
xlabel('East/West Slowness (s/deg)','fontweight','bold');
ylabel('North/South Slowness (s/deg)','fontweight','bold');
colormap(ritz);
c=colorbar('eastoutside','fontweight','bold');
set(c,'xaxislocation','top');
xlabel(c,'dB')
axis equal tight

end

function [cx,cy]=circle(r)
    ang=0:0.002:pi*2;
    cx=sin(ang)*r;
    cy=cos(ang)*r;
end

function [cx,cy]=circle2(r,steps)
    ang=0:2*pi/steps:pi*2;
    cx=sin(ang)*r;
    cy=cos(ang)*r;
end

