function [G,m_synth,tt,s,e]=easytomo(nrays,nblocks,vmean,vvar)
%EASYTOMO    Generates a random, simple tomography problem
%
%    Usage:    [G,m_synth,tt,s,e]=easytomo(nrays,model)
%              [G,m_synth,tt,s,e]=easytomo(nrays,nblocks,vmean,vvar)
%
%    Description:
%     [G,M_SYNTH,TT,S,E]=EASYTOMO(NRAYS,MODEL) draws NRAYS number of
%     raypaths through the velocity model MODEL.  Their total travel time
%     through the region and individual path lengths within each block are
%     then calculated and returned.  MODEL should have velocities scaled to
%     the block size (so if the blocks are 5x5km, then a velocity of 1
%     would mean 5kmps).  Blocks are assumed to be square.  The model does
%     not have to be square (ie the matrix can be MxN).  The outputs are
%     the kernel matrix G (contains the path length in each block for each
%     raypath with one raypath per row - this defines your data/model
%     relationship and is the hard part to work out), m_synth (1/velocity
%     aka the slowness for each block - this is what you want to resolve),
%     tt (the total travel time for each raypath through the region - this
%     is your data), and the starting and ending positions of the raypaths
%     in [x y] format, S and E.
%
%     [G,M_SYNTH,TT,S,E]=EASYTOMO(NRAYS,NBLOCKS,VMEAN,VVAR) creates a
%     region of NBLOCKS by NBLOCKS (blocks are square) with each block
%     having a random velocity within the range vmean+/-vvar blocks/sec.
%     NBLOCKS may be a 2 element vector to specify a non-square region.
%
%    Notes/Caveats:
%     - in the plot, green squares indicate the starting position of a ray
%       path and the red diamonds give the ending position
%     - no receivers or sources - its just random raypaths
%     - no ray path bending - just straight lines
%     - rays always start and end on the boundary of the region
%     - no error is added to the travel times
%
%    Examples:
%     % 15 rays through a 3x3 grid with velocity 10+/-1:
%     [G,m_synth,tt,s,e]=easytomo(15,3,10,1);
%
%    See also:

%     Version History:
%        Fall     2004 - initial version
%        Fall     2006 - made into a function, doc update
%        Fall     2008 - major code cleanup, renamed to EASYTOMO
%        Mar.  6, 2009 - added to SEIZMO
%        Sep.  6, 2009 - updated documentation
%        Mar. 24, 2012 - minor doc update, plot calls use handles now
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 02:50 GMT

% todo:

% act by number of inputs
if(nargin==4)
    % random velocity model
    vel=vmean+vvar*2*(rand(nblocks)-.5);
elseif(nargin==2)
    % user model - remember velocity must be normed to block size!
    vel=nblocks;
else
    error('bad number of inputs!')
end

% check input model
[m,n]=size(vel);

% get model parameters
m_synth=1./vel(:);

% plot velocity model
fh=figure;
ax=axes('parent',fh);
pcolor(ax,0:n,0:m,[vel zeros(m,1); zeros(1,n) 0]);
shading(ax,'flat');
axis(ax,'ij');
caxis(ax,[min(vel(:))-2*eps max(vel(:))+2*eps]);
colormap(ax,flipud(jet(1024)));
a=colorbar('peer',ax);
ylabel(a,'Velocity (blocks/sec)','fontsize',14,'fontweight','bold');
title(ax,'Velocity Grid and Raypaths','fontsize',14,'fontweight','bold');
set(ax,'fontsize',14,'fontweight','bold','xaxislocation','top');
hold(ax,'on');

% random position on side (0 to 1)
pos=rand(nrays,2);
%pos

% random side (no same side)
%
%        1
%      _____
%     |     |
%  3  |     |  4
%     |_____|
%        
%        2
%
s1=ceil(rand(nrays,1)*4);
s2=mod(s1+ceil(rand(nrays,1)*3)-1,4)+1;
%[s1 s2]

% start/end position
% (y,x) axis agrees with (row,column) indices in matrix:
%
%   (1,1) (1,2) (1,3)
%   (2,1) (2,2) (2,3)
%   (3,1) (3,2) (3,3)
%
%  -y
%
%   ^
%   |
%   |
%   +---->  +x
%  
s11=s1==1; s12=s1==2; s13=s1==3; s14=s1==4;
s21=s2==1; s22=s2==2; s23=s2==3; s24=s2==4;
sx=zeros(nrays,1); sy=sx; ex=sx; ey=sx;
sx(s11)=pos(s11,1)*n;
sx(s12)=pos(s12,1)*n; sy(s12)=m;
                      sy(s13)=pos(s13,1)*m;
sx(s14)=n;            sy(s14)=pos(s14,1)*m;
ex(s21)=pos(s21,2)*n;
ex(s22)=pos(s22,2)*n; ey(s22)=m;
                      ey(s23)=pos(s23,2)*m;
ex(s24)=n;            ey(s24)=pos(s24,2)*m;
s=[sx sy];
e=[ex ey];

% plot ray paths
plot(ax,[sx.'; ex.'],[sy.'; ey.'],'k','linewidth',2);
plot(ax,sx,sy,'s',...
    'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',12);
plot(ax,ex,ey,'d',...
    'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',7);
hold(ax,'off');

% geometry
dx=ex-sx;            % x length
dy=ey-sy;            % y length
phi=atan2(dx,dy);    % angle from down, increases ccw

% starting block (takes care of starting on corner too)
sbx=ceil(sx+eps*(sx==fix(sx) & dx>0));
sby=ceil(sy+eps*(sy==fix(sy) & dy>0));

% create straight line length matrix
G=zeros(nrays,m*n);
for i=1:nrays
    % plan of action:
    %  get the distance to each boundary
    %  sort the distances keeping the info on which type (x vs y)
    %  use the type info to get blocks along path
    %  use the distance info for length in each block along path
    
    % intercept points at block boundaries
    if(dx(i)>0)
        xbx=(ceil(sx(i)):ex(i))-sx(i);
    else
        xbx=sx(i)-(floor(sx(i)):-1:ex(i));
    end
    if(dy(i)>0)
        yby=(ceil(sy(i)):ey(i))-sy(i);
    else
        yby=sy(i)-(floor(sy(i)):-1:ey(i));
    end
    
    % distance to intercept points
    xbd=abs(xbx./sin(phi(i)));
    ybd=abs(yby./cos(phi(i)));
    
    % drop 0 distance
    xbd(xbd==0)=[];
    ybd(ybd==0)=[];
    
    % create distance-step matrix
    nx=numel(xbd); ny=numel(ybd);
    ds=[xbd.' ones(nx,1) zeros(nx,1); ybd.' zeros(ny,1) ones(ny,1)];
    ds=sortrows(ds);
    
    % lengths in blocks
    d=diff([0; ds(:,1)]);
    
    % steps
    ixs=[0; ds(:,2)]*sign(dx(i));
    iys=[0; ds(:,3)]*sign(dy(i));
    
    % crunch and combine repeats (happens at corner)
    r=find(d==0);
    if(~isempty(r))
        ixs(r+1)=ixs(r)+ixs(r+1);
        iys(r+1)=iys(r)+iys(r+1);
        d(r)=[]; ixs(r)=[]; iys(r)=[];
        disp('corner!'); % haven't actually seen this yet
    end
    
    % drop final block step (steps off grid)
    ixs(end)=[]; iys(end)=[];
    
    % indices of blocks
    ix=sbx(i)+cumsum(ixs);
    iy=sby(i)+cumsum(iys);
    
    % lengths in blocks
    %len=zeros(m,n);
    %len(iy+m*(ix-1))=d
    G(i,iy+m*(ix-1))=d;
end

% calculating the travel times for each ray
tt=G*m_synth;

end
