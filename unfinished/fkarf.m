function [arf]=fkarf(stla,stlo,smax,spts,s,baz,f,arf)
%ARF    Returns the fk array response function for a SEIZMO dataset

% add array response function to existing one or not
if(nargin<8); arf=zeros(spts,spts); end
if(~isequal(size(arf),[spts spts]))
    error('seizmo:fkarf:badARF',...
        'ARF stacking must use the same slowness grid!');
end

% angular frequency
w=2*pi*f;

% number of stations
nst=numel(stla);

% calculate array center
[cla,clo]=arraycenter(stla,stlo);

% calculate location vectors
% - temporary measure: get dist, az and get r from those
[dist,az]=vincentyinv(cla,clo,stla,stlo);
x=dist.*sin(az*pi/180);
y=dist.*cos(az*pi/180);

% break down slowness/baz to sx/sy
sx0=s*sin(baz*pi/180);
sy0=s*cos(baz*pi/180);

% make wavenumber arrays
sx=-smax:2*smax/(spts-1):smax;
kk0x=w*(sx(ones(spts,1),:)-sx0)/(6371*pi/180);
kk0y=w*(fliplr(sx(ones(spts,1),:))'-sy0)/(6371*pi/180);

% get array response
%  2*pi*i*(k-k0)r
% e
for i=1:nst
    arf=arf+real(exp(-0.5i*(kk0x*x(i)+kk0y*y(i)))).^2/nst;
end

% plot it
figure;
hold on
for i=[3 4 8 25]
    [x,y]=circle(1/i*6371*pi/180);
    plot(x,y,'k','linewidth',2);
end
[x,y]=circle2(1/3*6371*pi/180,12);
[x2,y2]=circle2(1/25*6371*pi/180,12);
plot([x; x2],[y; y2],'k','linewidth',2);
imagesc(sx,fliplr(sx),10*log10(arf)); set(gca,'ydir','normal');
set(gca,'fontweight','bold');
alpha(0.9)
hold off
title(['Array Response Function @ ' ...
    num2str(f) 'Hz (' num2str(1/f) 's)'],'fontweight','bold');
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
