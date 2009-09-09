function [x,y,p,a]=sensitivity(f,amp,vel)
% makes a surface wave sensitivity kernel near the receiver

% kernel map limits (make bigger than your region size)
w=1500; d=10;
x=-w:d:w;
n=numel(x);

% normalize amps
amp=amp/sum(amp);

% make sensitivity kernel
k=2*pi*6371*f/vel;
y=x(ones(n,1),:);
x=y.';
r=sqrt(x.^2+y.^2);
r(floor(n*n/2)+1)=d;
c=1./(4*r.^2*sqrt(2*pi*abs(sin(r./6371)))./(d^2));
p=zeros(n); a=zeros(n);
for i=1:numel(f)
    p=p-c.*amp(i).*k(i).^1.5.*sin(k(i).*(x+r)./6371+pi/4);
    a=a-c.*amp(i).*k(i).^1.5.*cos(k(i).*(x+r)./6371+pi/4);
end

figure;
surface(x,y,p);
colormap(jet(1024));
shading interp;
hold on
title('Phase Kernel')
figure;
surface(x,y,a);
colormap(jet(1024));
shading interp;
hold on
title('Amplitude Kernel')

end