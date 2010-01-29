function bb(fm, centerX, centerY, diam, ta, color)
%function bb(fm, centerX, centerY, diam, ta, color)
%
%	Draws beachball diagram of earthquake double-couple focal mechanism(s). S1, D1, and
%		R1, the strike, dip and rake of one of the focal planes, can be vectors of
%		multiple focal mechanisms.
%	fm - focal mechanism that is either number of mechnisms (NM) by 3 (strike, dip, and rake)
%		or NM x 6 (mxx, myy, mzz, mxy, mxz, myz - the six independent components of
%		the moment tensor). The strike is of the first plane, clockwise relative to north.
%		The dip is of the first plane, defined clockwise and perpedicular to strike, 
%		relative to horizontal such that 0 is horizontal and 90 is vertical. The rake is 
%		of the first focal plane solution. 90 moves the hanging wall up-dip (thrust),
%		0 moves it in the strike direction (left-lateral), -90 moves it down-dip
%		(normal), and 180 moves it opposite to strike (right-lateral).
%	centerX - place beachball(s) at position centerX
%	centerY - place beachball(s) at position centerY
%	diam - draw with this diameter.  If diam is zero, beachball
%		is drawn on stereonet.
%	ta - type of axis. If 0, this is a normal axis, if 1 it is a map axis. In case
%		of the latter, centerX and centerY are Lons and Lats, respectively.
%	color - color to use for quadrants of tension; can be a string, e.g. 'r'
%		'b' or three component color vector, [R G B].
%


[ne,n] = size(fm);
if n == 6
	for j = 1:ne
		[s1(j),d1(j),r1(j)] = mij2sdr(fm(j,1),fm(j,2),fm(j,3),fm(j,4),fm(j,5),fm(j,6));
	end
else
	s1 = fm(:,1);
	d1 = fm(:,2);
	r1 = fm(:,3);
end

r2d = 180/pi;
d2r = pi/180;
ampy = cos(mean(centerY)*d2r);

if ne > 1
	[ds,i] = sort(diam,1,'descend');
	diam = diam(i);
	s1 = s1(i);
	d1 = d1(i);
	r1 = r1(i);
	centerX = centerX(i);
	centerY = centerY(i);
end

mech = zeros(ne,1);
j = find(r1 > 180);
r1(j) = r1(j) - 180;
mech(j) = 1;
j = find(r1 < 0);
r1(j) = r1(j) + 180;
mech(j) = 1;

% Get azimuth and dip of second plane
[s2,d2,r2] = AuxPlane(s1,d1,r1);

if diam(1) > 0
	hold on
end
for ev = 1:ne
	S1 = s1(ev);
	D1 = d1(ev);
	S2 = s2(ev);
	D2 = d2(ev);
	P = r1(ev);
	CX = centerX(ev);
	CY = centerY(ev);
	D = diam(ev);
	M = mech(ev);

if M > 0
   P = 2;
else
   P = 1;
end

if D1 >= 90
   D1 = 89.9999;
end
if D2 >= 90
   D2 = 89.9999;
end

phi = 0:.01:pi;
d = 90 - D1;
m = 90;
l1 = sqrt(d^2./(sin(phi).^2 + cos(phi).^2 * d^2/m^2));

d = 90 - D2;
m = 90;
l2 = sqrt(d^2./(sin(phi).^2 + cos(phi).^2 * d^2/m^2));

if D == 0
   stereo(phi+S1*d2r,l1,'k')
   hold on
   stereo(phi+S2*d2r,l2,'k')
end

inc = 1;
[X1,Y1] = pol2cart(phi+S1*d2r,l1);
if P == 1
   lo = S1 - 180;
   hi = S2;
   if lo > hi
      inc = -inc;
   end
   th1 = S1-180:inc:S2;
   [Xs1,Ys1] = pol2cart(th1*d2r,90*ones(1,length(th1)));
   [X2,Y2] = pol2cart(phi+S2*d2r,l2);
   th2 = S2+180:-inc:S1;
else
   hi = S1 - 180;
   lo = S2 - 180;
   if lo > hi
      inc = -inc;
   end
   th1 = hi:-inc:lo;
   [Xs1,Ys1] = pol2cart(th1*d2r,90*ones(1,length(th1)));
   [X2,Y2] = pol2cart(phi+S2*d2r,l2);
   X2 = fliplr(X2);
   Y2 = fliplr(Y2);
   th2 = S2:inc:S1;
end
[Xs2,Ys2] = pol2cart(th2*d2r,90*ones(1,length(th2)));

X = cat(2,X1,Xs1,X2,Xs2);
Y = cat(2,Y1,Ys1,Y2,Ys2);

if D > 0
   X = ampy*X * D/90 + CY;
   Y = Y * D/90 + CX;
   phid = 0:.01:2*pi;
   [x,y] = pol2cart(phid,90);
   xx = x*D/90 + CX;
   yy = ampy*y*D/90 + CY;
   if ta == 0
      fill(xx,yy,'w')
      fill(Y,X,color)
      line(xx,yy,'color','k','linewidth',0.2);
   else
      fillm(yy,xx,'w')
      fillm(X,Y,color)
      linem(yy,xx,'color','k','linewidth',0.2);
   end
else
   if ta == 0
      fill(X,Y,color)
   else
      fillm(Y,X,color)
   end
   view(90,-90)
end

end

function [strike, dip, rake] = AuxPlane(s1,d1,r1);
%function [strike, dip, rake] = AuxPlane(s1,d1,r1);
% Get Strike and dip of second plane, adapted from Andy Michael bothplanes.c
r2d = 180/pi;

z = (s1+90)/r2d;
z2 = d1/r2d;
z3 = r1/r2d;
%/* slick vector in plane 1 */
sl1 = -cos(z3).*cos(z)-sin(z3).*sin(z).*cos(z2);
sl2 = cos(z3).*sin(z)-sin(z3).*cos(z).*cos(z2);
sl3 = sin(z3).*sin(z2);
[strike, dip] = strikedip(sl2,sl1,sl3);

n1 = sin(z).*sin(z2);  %/* normal vector to plane 1 */
n2 = cos(z).*sin(z2);
n3 = cos(z2);
h1 = -sl2; %/* strike vector of plane 2 */
h2 = sl1;
%/* note h3=0 always so we leave it out */

z = h1.*n1 + h2.*n2;
z = z./sqrt(h1.*h1 + h2.*h2);
z = acos(z);

rake = zeros(size(strike));
j = find(sl3 > 0);
rake(j) = z(j)*r2d;
j = find(sl3 <= 0);
rake(j) = -z(j)*r2d;

function [strike, dip] = strikedip(n, e, u)
%function [strike, dip] = strikedip(n, e, u)
%       Finds strike and dip of plane given normal vector having components n, e, and u
%

% Adapted from Andy Michaels stridip.c

r2d = 180/pi;

j = find(u < 0);
n(j) = -n(j);
e(j) = -e(j);
u(j) = -u(j);

strike = atan2(e,n)*r2d;
strike = strike - 90;
while strike >= 360
        strike = strike - 360;
end
while strike < 0
        strike = strike + 360;
end

x = sqrt(n.^2 + e.^2);
dip = atan2(x,u)*r2d;

function hpol = stereo(theta,rho,line_style)
%function hpol = stereo(theta,rho,line_style)

if nargin < 1
    error('Requires 2 or 3 input arguments.')
elseif nargin == 2 
    if isstr(rho)
        line_style = rho;
        rho = theta;
        [mr,nr] = size(rho);
        if mr == 1
            theta = 1:nr;
        else
            th = (1:mr)';
            theta = th(:,ones(1,nr));
        end
    else
        line_style = 'auto';
    end
elseif nargin == 1
    line_style = 'auto';
    rho = theta;
    [mr,nr] = size(rho);
    if mr == 1
        theta = 1:nr;
    else
        th = (1:mr)';
        theta = th(:,ones(1,nr));
    end
end
if isstr(theta) | isstr(rho)
    error('Input arguments must be numeric.');
end
if ~isequal(size(theta),size(rho))
    error('THETA and RHO must be the same size.');
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold on;
    maxrho = max(abs(rho(:)));
    hhh=plot([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho]);
    set(gca,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    set(gca,'xlim',[-90 90])
    set(gca,'ylim',[-90 90])
    v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
%    if rticks > 5   % see if we can reduce the number
%        if rem(rticks,2) == 0
%            rticks = rticks/2;
%        elseif rem(rticks,3) == 0
%            rticks = rticks/3;
%        end
%    end
	 rticks = 1;

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~isstr(get(cax,'color')),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end

% draw radial circles
    c82 = cos(82*pi/180);
    s82 = sin(82*pi/180);
    rinc = (rmax-rmin)/rticks;
    for i=(rmin+rinc):rinc:rmax
        hhh = plot(xunit*i,yunit*i,ls,'color',tc,'linewidth',1,...
                   'handlevisibility','off');
%        text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
%            ['  ' num2str(i)],'verticalalignment','bottom',...
%            'handlevisibility','off')
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = (1:6)*2*pi/12;
    cst = cos(th); snt = sin(th);
%    cs = [-cst; cst];
%    sn = [-snt; snt];
%    plot(rmax*cs,rmax*sn,ls,'color',tc,'linewidth',1,...
%         'handlevisibility','off')

% annotate spokes in degrees
    rt = 1.1*rmax;
    for i = 1:length(th)
        text(rt*cst(i),rt*snt(i),int2str(i*30),...
             'horizontalalignment','center',...
             'handlevisibility','off');
        if i == length(th)
            loc = int2str(0);
        else
            loc = int2str(180+i*30);
        end
        text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
             'handlevisibility','off')
    end

% set view to 2-D
    view(2);
% set axis limits
    axis(rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

% transform data to Cartesian coordinates.
xx = rho.*cos(theta);
yy = rho.*sin(theta);

% plot data on top of grid
if strcmp(line_style,'auto')
    q = plot(xx,yy);
else
    q = plot(xx,yy,line_style);
end
if nargout > 0
    hpol = q;
end
if ~hold_state
    set(gca,'dataaspectratio',[1 1 1]), axis off; set(cax,'NextPlot',next);
end
set(get(gca,'xlabel'),'visible','on')
set(get(gca,'ylabel'),'visible','on')


function [str,dip,rake] = mij2sdr(mxx,myy,mzz,mxy,mxz,myz)
%function [str,dip,rake] = mij2sdr(mxx,myy,mzz,mxy,mxz,myz)
%
%   INPUT
%	mij - siz independent components of the moment tensor
%
%   OUTPUT
%	str - strike of first focal plane (degrees)
%	dip - dip of first focal plane (degrees)
%	rake - rake of first focal plane (degrees)
%
%Adapted from code, mij2d.f, created by Chen Ji and given to me by Gaven Hayes.
%

a = [mxx mxy mxz; mxy myy myz; mxz myz mzz];
[V,d] = eig(a);

D = [d(3,3) d(1,1) d(2,2)];
V(2:3,1:3) = -V(2:3,1:3);
V = [V(2,3) V(2,1) V(2,2); V(3,3) V(3,1) V(3,2); V(1,3) V(1,1) V(1,2)];

IMAX = find(D == max(D));
IMIN = find(D == min(D));
AE = (V(:,IMAX)+V(:,IMIN))/sqrt(2.0);
AN = (V(:,IMAX)-V(:,IMIN))/sqrt(2.0);
AER = sqrt(AE(1)^2+AE(2)^2+AE(3)^2);
ANR = sqrt(AN(1)^2+AN(2)^2+AN(3)^2);
AE = AE/AER;
AN = AN/ANR;
if (AN(3) <= 0.)
	AN1 = AN;
	AE1 = AE;
else
	AN1 = -AN;
	AE1 = -AE;
end 
[ft,fd,fl] = TDL(AN1,AE1);
str = 360 - ft;
dip = fd;
rake = 180 - fl;

function [FT,FD,FL] = TDL(AN,BN)
XN=AN(1);
YN=AN(2);
ZN=AN(3);
XE=BN(1);
YE=BN(2);
ZE=BN(3);
AAA=1.0E-06;
CON=57.2957795;
if (abs(ZN) < AAA)
	FD=90.;
	AXN=abs(XN);
	if (AXN > 1.0) 
		AXN=1.0;
	end
	FT=asin(AXN)*CON;
	ST=-XN;
	CT=YN;
	if (ST >= 0. & CT < 0) 
		FT=180.-FT;
	end
	if (ST < 0. & CT <= 0) 
		FT=180.+FT;
	end
	if (ST < 0. & CT > 0) 
		FT=360.-FT;
	end
	FL=asin(abs(ZE))*CON;
	SL=-ZE;
	if (abs(XN) < AAA) THEN
		CL=XE/YN;
	else
		CL=-YE/XN;
	end 
	if (SL >= 0. & CL < 0) 
		FL=180.-FL;
	end
	if (SL < 0. & CL <= 0) 
		FL=FL-180.;
	end
	if (SL < 0. & CL > 0) 
		FL=-FL;
	end
else
	if (-ZN > 1.0) 
		ZN=-1.0;
	end
	FDH=acos(-ZN);
	FD=FDH*CON;
	SD=sin(FDH);
	if  (SD == 0)
		return;
	end 
	ST=-XN/SD;
	CT=YN/SD;
	SX=abs(ST);
	if (SX > 1.0) 
		SX=1.0;
	end
	FT=asin(SX)*CON;
	if (ST >= 0. & CT < 0) 
		FT=180.-FT;
	end
	if (ST < 0. & CT <= 0) 
		FT=180.+FT;
	end
	if (ST < 0. & CT > 0) 
		FT=360.-FT;
	end
	SL=-ZE/SD;
	SX=abs(SL);
	if (SX > 1.0) 
		SX=1.0;
	end
	FL=asin(SX)*CON;
	if (ST == 0) THEN
		CL=XE/CT;
	else
		XXX=YN*ZN*ZE/SD/SD+YE;
		CL=-SD*XXX/XN;
		if (CT == 0) 
			CL=YE/ST;
		end
	end 
	if (SL >= 0. & CL < 0) 
		FL=180.-FL;
	end
	if (SL < 0. & CL <= 0) 
		FL=FL-180.;
	end
	if (SL < 0. & CL > 0) 
		FL=-FL;
	end
end 
