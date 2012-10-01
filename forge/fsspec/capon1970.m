function [data]=capon1970()
%CAPON1970    Artificial data for HiRes fk testing
%
%    Usage:    data=capon1970()
%
%    Description:
%     DATA=CAPON1970() returns a SEIZMO dataset with waveforms roughly
%     equal to the artificial data published by Capon in 1970.  This is for
%     testing and comparing slowness spectra.
%
%    Notes:
%     - Reference:
%        Capon 1970, Analysis of Rayleigh-Wave Multipath Propagation at
%        LASA, BSSA, Vol. 60, No. 5, pp. 1701-1731
%
%    Examples:
%     % Capon estimation of the 40s wave:
%     plotfss(fssavg(fss(capon1970,50,101,[1/45 1/35],'m','capon')));
%
%    See also: LASA, FSS, PLOTFSS

% station locations
[st,nm]=lasa;
[e,n]=geographic2enu(st(:,1),st(:,2),0,st(1,1),st(1,2),0);

% planewave backazimuths
baz=[310 340];
vph=3.7;

% integer delays
D1=-round((n*cosd(baz(1))+e*sind(baz(1)))/vph);
D2=-round((n*cosd(baz(2))+e*sind(baz(2)))/vph);

% dispersive wave train
k=1:1000; T=1;
g=nan(1000,21);
for i=1:21
    g(:,i)=f(k-100-D1(i),T)+f(k-300-D2(i),T);
end

% covert to a seizmo dataset
[data{1,1:21}]=deal((1:1000)');
data(2,:)=mat2cell(g,1000,ones(1,21));
data=bseizmo(data{:});
data=ch(data,'kstnm',nm,'stla',st(:,1),'stlo',st(:,2));

end

function [v]=f(k,T)
f0=.025; f1=.050;
v=sin(2*pi*f0*k*T+(f1-f0)*pi/600*(k*T).^2);
v(k<1 | k>600)=0;
end
