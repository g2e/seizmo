% Function anaseis - Compute analytical seismograms for different source types.
%
%   call: anaseis(stype) 
%
%   stype : source type. One of: 'DC_XY' 'DC_XZ' 'DC_YZ' 'EXPL'
%           default: 'DC_YZ'
%
%   See equations 4.32-3, page 81, Aki and Richards, 1st edition.
%
%   Note:
%     The R-component in the seismogram plots is amplified by
%     sqrt(3)^3 in order to make P-waves more visible
%
%   Changes:
%     Switched x- and y-seismogram components for obtaining plausible P-seismograms
%
%   Known bugs:
%     (1) For vertical incidence, the PHI component is wrong (division by zero?)
%
% G. Jahnke, 31. Oct. 2004 - 04. Nov. 2004
function anaseis(stype)
  % Set up model and source parameters
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vp=5.;                % P velocity
  vs=5./sqrt(3);        % S velocity
  rho=2.;               % density
  sLoc=[0.,0.,0.];      % Earthquake location (km)
  f0=2;                 % mean frequency
  dt=.01;               % time increment
  tmax=6.;              % seismogram length
  nt=ceil(tmax/dt);     % number of samples

  if(~exist('stype'))
    stype='DC_YZ';
    fprintf(1,'\n  * Using default source type ''%s''\n\n',stype);
  end
  
  [rLoc1,ang1]=defineReceiverRing(1);   % 1: xz-plane
  [rLoc2,ang2]=defineReceiverRing(2);   % 2: xy-plane
  [rLoc3,ang3]=defineReceiverRing(3);   % 3: yz-plane   

  rLoc=[rLoc1 rLoc2 rLoc3];            % concatenate receiver rings
  ang=[ang1 ang2 ang3];                % and angles
  no_rings=3;                          % number of receiver rings (for plotting)

  checkValidSrcRecPos(sLoc,rLoc);       % check if sLoc<>rLoc
  plotSrcRec(sLoc,rLoc,no_rings,stype); % plot source and receivers
  plotRadPat(stype,2.0); % plot P and S radiation pattern
 
  % Note: 'Seis' is a struct array! (Seis.x Seis.y Seis.z ...)
  %
  Seis=computeGreensFunctions(sLoc,rLoc,nt,dt,rho,vp,vs,stype);
  Seis=addSourceTimeFunction(Seis,f0,dt); % perform convolution
  plotSeismograms(Seis,dt,f0,no_rings);   % plot seismograms
  return
  
% Create receiver ring.
%
% mode=1 : ring in x/z plane (y=0)
%      2 : ring in x/y plane (z=0)
%      3 : ring in y/z plane (x=0)
%
function [rLoc,ang]=defineReceiverRing(mode)
  nr =16;                 % number of receivers
  rad=10;                 % radius of the ring [km]
  ang=[0:nr-1]*2*pi/nr;   % angle in selected plane
  v1=rad*sin(ang);        % 1st coord in ring-plane
  v2=rad*cos(ang);        % 2nd coord in ring-plane
  v3=zeros(1,nr);         % the constant coordinate

  switch(mode)
    case 1                   % ring in xz plane (y=0)
      rLoc=[v1', v3', v2']';
    case 2                   % ring in xy plane (z=0)
      rLoc=[v1', v2', v3']';
    case 3                   % ring in yz plane (x=0)
      rLoc=[v3', v2', v1']';
    otherwise
      error(sprintf(' defineReceiverRing() : Unknown mode %d (sorry)\n',mode));
  end
return

% Compute the greens functions for all receivers.
%
function Seis= computeGreensFunctions(sLoc,rLoc,nt,dt,rho,vp,vs,stype);
  cp=1/(4*pi*rho*vp^3); % derived constants
  cs=1/(4*pi*rho*vs^3); % Note: cs/cp is also the S/P amplitude ratio (3*sqrt(3))

  nr=size(rLoc,2); % number of receivers

  xr=rLoc(1,:); yr=rLoc(2,:); zr=rLoc(3,:);
  xs=sLoc(1);   ys=sLoc(2);   zs=sLoc(3);

  % allocate output seismograms
  sx=zeros(nr,nt);  sy=zeros(nr,nt);  sz=zeros(nr,nt);
  sr=zeros(nr,nt);  st=zeros(nr,nt);  sp=zeros(nr,nt);

  % for each receiver  
  for i=1:nr
    phi  =atan2( (xr(i)-xs),(yr(i)-ys));           % angle in xy plane ( old: atan(val) )
    rr   =norm([ xr(i)-xs , yr(i)-ys, zr(i)-zs]);  % src-receiver distance

    theta=acos((zr(i)-zs)/rr);                     % vertical angle (???)

    % unit vectors in source-receiver system
    ur=[ sin(theta)*cos(phi), sin(theta)*sin(phi),  cos(theta)]'; %'
    ut=[ cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)]'; %'
    up=[ -sin(phi)          , cos(phi)           ,           0]'; %'

    % TODO: check why the wrong components are assigned!!!
    % Probably this is due to the wrong order in the atan2() term above.
    % The next line correct this by flipping x- and y-comp.
    ur=ur([2,1,3]); % correct switched x- and y-components
    
    % P & S wave amplitudes
    switch(stype)
      case 'DC_XY'
        Ap =cp/rr*  sin(theta)*sin(2*phi);
        Ast=cs/rr*  sin(theta)*cos(2*phi);
        Asp=cs/rr*(-sin(theta)*sin(2*phi));
      case 'DC_XZ'
        phi=mod(phi+pi/2,2*pi);
	Ap =cp/rr*  sin(2*theta)*cos(phi);
        Ast=cs/rr*  cos(2*theta)*cos(phi);
        Asp=cs/rr*(-cos(  theta)*sin(phi));
      case 'DC_YZ'
        Ap =cp/rr*  sin(2*theta)*cos(phi);
        Ast=cs/rr*  cos(2*theta)*cos(phi);
        Asp=cs/rr*(-cos(  theta)*sin(phi));
      case 'EXPL'
        Ap =cp/rr; % CHECK THIS!
        Ast=0;
        Asp=0;
      otherwise
        error(sprintf('\n   * Unknown source type ''%s''!\n',stype));
    end

    % arrival times
    ttp=rr/vp; ntp=round(ttp/dt);
    tts=rr/vs; nts=round(tts/dt);

    A =Ap*ur;          % project P amplitude to xyz components
    sr(i,ntp)=Ap;      % radial component
    sx(i,ntp)=A(1);    % Ax=Ap*ur*ux; sx(i,ntp)=Ax; %% with ux=[1; 0; 0];
    sy(i,ntp)=A(2);    % Ay=Ap*ur*uy; sy(i,ntp)=Ay; %% with uy=[0; 1; 0];
    sz(i,ntp)=A(3);    % Az=Ap*ur*uz; sz(i,ntp)=Az; %% with uz=[0; 0; 1];

    A=Ast*ut;          % project S_theta amplitude to xyz components
    st(i,ntp)=Ast;     % theta component
    sx(i,nts)=A(1);    % Ax=Ast*ut*ux; sx(i,nts)=Ax;
    sy(i,nts)=A(2);    % Ay=Ast*ut*uy; sy(i,nts)=Ay;
    sz(i,nts)=A(3);    % Az=Ast*ut*uz; sz(i,nts)=Az;

    A=Asp*up;                     % project S_phi amplitude to xyz components
    sp(i,nts)=Asp;                % phi component
    sx(i,nts)=sx(i,nts)+A(1);     % Ax=Asp*up*ux; sx(i,nts)=sx(i,nts)+Ax;
    sy(i,nts)=sy(i,nts)+A(2);     % Ay=Asp*up*uy; sy(i,nts)=sy(i,nts)+Ay;
    sz(i,nts)=sz(i,nts)+A(3);     % Az=Asp*up*uz; sz(i,nts)=sz(i,nts)+Az;
    Seis.ttp(i)=ttp;
    Seis.tts(i)=tts; % preserve P & S travel times 
  end

  Seis.x=sx;  Seis.y=sy;  Seis.z=sz;
  Seis.r=sr;  Seis.t=st;  Seis.p=sp;
return

% Check if no receiver is at source position
%
function checkValidSrcRecPos(sLoc,rLoc)
  % create src-rec difference vectors
  a=[rLoc(1,:)'-sLoc(1) [rLoc(2,:)'-sLoc(2)] rLoc(3,:)'-sLoc(3)]; %'

  % for all src-rec vectors
  for i=1:size(a,1) 
    if( norm(abs(a(i,:)))==0) % invalid receiver position found
      error(sprintf('\n   * Not allowed: Receiver %d at source position!\n',i));
    end
  end
return

% plot P and S radiation pattern
%
function plotRadPat(stype,sc)
  radpat(stype,sc);
return


% Create 3d-plot receivers and source in figure 1.
%
function plotSrcRec(sLoc,rLoc,no_rings,stype)

  xs=sLoc(1);   ys=sLoc(2);   zs=sLoc(3);
  xr=rLoc(1,:); yr=rLoc(2,:); zr=rLoc(3,:);

  figure(1),clg; % create figure for source-receiver configuration
  set(gca,'FontSize',16);
  plot3(xs,ys,zs,'*'); hold on         % plot source
  % label with receiver numbers
  for i=1:length(xr)
    if(plotlabel(i,length(xr),no_rings)==1) % avoid multiple labels at same position
      text(xr(i),yr(i),zr(i),sprintf(' %02d',i),'FontSize',16);
    end
  end

  if(no_rings==3)
    is=1; ie=length(xr)/3;
    plot3([xr(is:ie) xr(is)],[yr(is:ie) yr(is)],[zr(is:ie) zr(is)],'r-','LineWidth',2);
    is=  length(xr)/3+1;  ie=2*length(xr)/3;
    plot3([xr(is:ie) xr(is)],[yr(is:ie) yr(is)],[zr(is:ie) zr(is)],'k-','LineWidth',2);
    is=2*length(xr)/3+1;  ie=length(xr);
    plot3([xr(is:ie) xr(is)],[yr(is:ie) yr(is)],[zr(is:ie) zr(is)],'b-','LineWidth',2);
  elseif(no_rings==1)
    plot3(xr,yr,zr,'k-'), hold on
  end
  plot3(xr,yr,zr,'kx','MarkerSize',10,'LineWidth',2), hold on

  title(sprintf(' Source Type: %s ',strrep(stype,'_','\_')));
  xlabel('x'),  ylabel('y'),  zlabel('z')
  grid,  axis equal,  hold off
return

% Quick Hack
% works only for 12 or 16 receivers per  receiver ring!!!
function flag=plotlabel(i,len,no_rings)
  flag=1;
  if( len==36 & no_rings==3 & (i==4 | i==10 | i==19 | i==25 | i== 28 | i==34))
    flag=0;
  end  
  if( len==48 & no_rings==3 & (i==5 | i==13 | i==25 | i==33| i==37 | i==45))
    flag=0;
  end  
return

% Do convolution of Greens-functions with source time function.
%
function Seis=addSourceTimeFunction(Seis,f0,dt)

  sx=Seis.x;  sy=Seis.y; sz=Seis.z;
  sr=Seis.r;  st=Seis.t; sp=Seis.p;
  nr=size(sx,1);

  stf=ricker(f0,dt);                      % create source time function

  ssx=zeros(nr,length(sx(1,:))+length(stf)-1);
  ssy=ssx;  ssz=ssx;
  ssr=ssx;  sst=ssx;  ssp=ssx;

  for i=1:nr
     ssx(i,:)=conv(sx(i,:),stf');
     ssy(i,:)=conv(sy(i,:),stf');
     ssz(i,:)=conv(sz(i,:),stf');
     ssr(i,:)=conv(sr(i,:),stf');
     sst(i,:)=conv(st(i,:),stf');
     ssp(i,:)=conv(sp(i,:),stf');
  end
  
  Seis.x=ssx;  Seis.y=ssy;  Seis.z=ssz;
  Seis.r=ssr;  Seis.t=sst;  Seis.p=ssp;
return

% Plot seismograms from struct array 'Seis' in figure 2.
%
function plotSeismograms(Seis,dt,f0,no_rings)
  nr=size(Seis.x,1); % number of receivers
  nt=size(Seis.x,2); % number of samples
  amp=1;             % seismogram amplification factor
  taxis=[0:nt-1]*dt; % create time axis
  col=getcol(nr,no_rings); % define colors for all traces

  titstr=strvcat('x-comp','y-comp','z-comp','radial','theta','phi');

  % amplify P by S/P amplitude ratio which is: cs/cp=sqrt(3)*3
  Seis.r=Seis.r*sqrt(3)^3; 

  % 'ce': norm preserving ratios between components;
  % 'ci':  norm each component individually
  Seis=normSeis(Seis,'ce');

  % store seismograms in [NRxNTxNC] array, for convenient plotting
  % ( ts := ts(irec,ismp,icomp) )
  ts=seis2tensor(Seis);
  ttp=Seis.ttp; tts=Seis.tts;
  
  figure(2);      % open seismogram window
  %PCOMP=[1 2 3]; % plot x/y/z component
  %PCOMP=[4 5 6]; % plot r/t/p component
  PCOMP=[1:6];    % plot x/y/z/r/t/p component
  ik=0; % increment index

  % find maximum non-zero sample-no. (ilen) in tensor 'ts'
  %
  [inr,ilen,icmp]=ind2sub(size(ts),max(find(ts~=0)));
  ilen=min(ilen+round(f0/4/dt),size(Seis.x,2)); % include samples at the end of the trace
  
  for k=PCOMP % for each selected component
    ik=ik+1;
    subplot(100+length(PCOMP)*10++ik); % e.g. subplot(131)
    for i=1:nr
      len=find(max(ts(i,:,k))~=0);
     plot(taxis(1:ilen),amp*ts(i,1:ilen,k)+i,col(i)), hold on
    end

    plot(ttp,[1:nr],'b:'); plot(tts,[1:nr],'r:'); % mark P- and S-traveltime
    title(titstr(k,:));
    set(gca,'YTick',1:nr);axis([taxis(1) taxis(ilen) 0 nr+1])
    hold off
  end
return

function col=getcol(nr,no_rings)
  switch(no_rings)
    case 1
      col(1:nr)='b';
    case 3
      col(1:nr)='b';
      col(1:end/3)='r';
      col(end/3+1:2*end/3)='k';
  end
return



% Store seismograms in 3-dim array, for convenient plotting.
%
% 1st index: receiver no.
% 2nd index: sample index
% 3rd index: component [x,y,z,r,t]
%
function ts=seis2tensor(Seis);
  ts(:,:,1)=Seis.x;  ts(:,:,2)=Seis.y;  ts(:,:,3)=Seis.z;
  ts(:,:,4)=Seis.r;  ts(:,:,5)=Seis.t;  ts(:,:,6)=Seis.p;
return

% Create ricker wavelet for given parameters,
% taken from Lars Ceranna''s matlab library.
%
function [ft] = ricker(fs,dt,time,ts)
  % ricker.........ricker-wavelet-like source-signal
  %
  %  call: [ft]=ricker(fs,dt,time,ts)
  %    ft - source excitation [sec]
  %         [length(time)+1 x 1]
  %
  %    fs - source frequency [Hz]
  %         [scalar]
  %    dt - sampling intervall [sec]
  %         [scalar]
  %  time - time-vector [sec]
  %         [1 x M] or [M x 1]
  %    ts - origin time [sec]
  %         [scalar]
  %
  %  d/dx^2(exp{-x^2})=-(2-4x^2)exp{-x^2}; x=wt
  %  Ricker-Wavelet: -d/dx^2(exp{-x^2})
  %
  %  !length(output)=length(time)+1
  %
  %  lc, 01.10.1998
  %
  t=[0:dt:3/(2*fs)]'; %'
  if(nargin == 2)
    ts=0;
    time=[0:dt:6/fs+dt]'; %'
  elseif nargin == 3
    ts=0;
  end
  fq=2*(1-2*pi^2*fs^2*t.^2).*exp(-pi^2*fs^2*t.^2);
  if (length([t; t]) > length(time))
    error('source frequency is too small, use a longer time-window')
  end
  fqq=[flipud(fq); fq];
  ft=[zeros(round(ts/dt)+1,1); fqq; zeros(length(time)-length(fqq)-round(ts/dt),1)];
return

% Normalize traces in struct-array Seis.
%
% flag='ci' : norm components individually
% flag='ce' : norm components equally
%      (preserves ratios between components)
%
function Seis=normSeis(Seis,flg)

  switch(flg)
    case 'ce'
      % global maximum of all traces
      gmax=max(max(abs([Seis.x Seis.y Seis.z])));
      Seis.x=normalize(Seis.x,gmax);
      Seis.y=normalize(Seis.y,gmax);
      Seis.z=normalize(Seis.z,gmax);
      gmax=max(max(abs([Seis.r Seis.t Seis.p])));
      Seis.r=normalize(Seis.r,gmax);
      Seis.t=normalize(Seis.t,gmax);
      Seis.p=normalize(Seis.p,gmax);
    case 'ci'
      Seis.x=normalize(Seis.x); Seis.y=normalize(Seis.y);
      Seis.z=normalize(Seis.z); Seis.r=normalize(Seis.r);
      Seis.t=normalize(Seis.t); Seis.p=normalize(Seis.p);
    otherwise
      error(sprintf(' normSeis() : Unknown mode %s (sorry)\n',flg));
  end
return

% normalize time series
%
% normalize(mat)      : set maximum of 'mat' to one
% normalize(mat,gmax) : divide 'mat' by 'gmax'
%
%
function mat=normalize(mat,gmax)

  mmax=max(max(abs(mat)));

  if(exist('gmax') & gmax>0)
    mat=mat/gmax;
  elseif(mmax>0)
    mat=mat/mmax;
  else
    fprintf(1,'  * normalize() : Can not normalize flat trace.\n');
  end
return

% Function: radpat(type,sc)
%
% Plots double-couple or explosion radiation pattern
%
% type: determines the type and orientation.
%       possible values:
%       'DC_XZ' 'DC_XY' 'DC_YZ' or 'EXPL'
%
% sc: scale factor, controls total plot size
%     default: sc=1
%
% For a double couple source, the p-radiation is blue
%   and the s-radiation is red.
%
% G. Jahnke, 04. Nov. 2004

% ------------------------------------------------------------------------
% Note for different notations in Matlab and e.g. Aki&Richards:
%
% Matlab definition of spherical coordinates [R,THETA,PHI]:
%  "THETHA is the counterclockwise angle in the xy plane measured from the
%    positive x axis.  PHI is the elevation angle from the xy plane.
%
% Common (e.g. Aki&Richards) definition:
%
%  "phi is the clockwise angle in the xy plane measured from the
%    positive x axis.  theta is the vertical angle starting from
%    the y axis"
%
% Matlab > Aki&Richards: [r,theta,phi]= [R,phi,pi/2-THETA]
%
%
function radpat(stype,sc)
  if(~exist('sc'))sc=1;end
  RES=.2;              % increase RES for faster plotting
  range=[-2:RES:2]*sc; % 3D volume to process

  [x,y,z]=meshgrid(range);

  % transform coordinate system
  switch(stype)
    case 'DC_XZ'
     [phi,theta,r] = cart2sph(x,y,z); % (see man cart2sph for order of return values)
     theta=pi/2-theta;                % correct theta (see comment in header)
     % no phi corrections necessary
   case 'DC_YZ'
     [phi,theta,r] = cart2sph(x,y,z);
     theta=pi/2-theta;
     phi=phi+pi/2;  % rotate 90deg in xy-plane
   case 'DC_XY'
     [phi,theta,r] = cart2sph(x,z,y);
     theta=pi/2-theta;
     %     phi=phi+pi/8;  % rotate 22.5deg in xy-plane
   case 'EXPL'
     [phi,theta,r] = cart2sph(x,z,y);
     theta=pi/2-theta;
    otherwise
      error('not yet implemented (sorry)');
  end

  % Compute radiation pattern
  %
  if(strcmp(stype(1:2),'DC'))
    Ap =1./r.*  sin(2*theta).*cos(phi);    % P radiation
    Ast=1./r.*  cos(2*theta).*cos(phi);    % S in theta-dir
    Asp=1./r.*(-cos(  theta).*sin(phi));   % S in phi-dir
  elseif(strcmp(stype,'EXPL'))
    Ap =0.5/r;    % P radiation
    Ast=r.*0;    % no S in theta-dir
    Asp=r.*0;    % no S in phi-dir
  end
  As=sqrt(Ast.^2+ Asp.^2);               % absolute S-amplitude

  p= patch(isosurface(x, y, z, abs(Ap),+.5/sc));
  s= patch(isosurface(x, y, z, As,+.5/sc)); alpha(s,.75); % activate alpha-blending

  isonormals(x,y,z,As,s);  isonormals(x,y,z,Ap,p);

  set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
  set(s, 'FaceColor', 'red', 'EdgeColor', 'none');

  daspect([1 1 1]);  view(3);
  camlight; lighting phong;

  xlabel('X'), ylabel('Y'), zlabel('Z');
return
