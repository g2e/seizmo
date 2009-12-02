function [h]=plot_mlocate_ellipsoids_in_cartesian(varargin)
%PLOT_MLOCATE_ELLIPSOIDS_IN_CARTESIAN  Plots MLOCATE confidence ellipsoids
%
% Usage: h=plot_mlocate_ellipsoids(files)

% parse input
files=onefilelist(varargin);

% colormap
cint=256;
repeatint=100;
cmap=hsv(cint);

% wireframe parameters
wflatlinelatint=15;
wflonlinelonint=15;
%wfdepint=100;
wflatlinelonint=1;
wflonlinelatint=1;

% make wireframe
wflatlinelats=-90:wflatlinelatint:90;
wflatlinelons=-180:wflatlinelonint:180;
wflatlinelat=wflatlinelats(ones(numel(wflatlinelons),1),:);
wflatlinelon=wflatlinelons(ones(numel(wflatlinelats),1),:).';
wflonlinelats=-90:wflonlinelatint:90;
wflonlinelons=-180:wflonlinelonint:180;
wflonlinelat=wflonlinelats(ones(numel(wflonlinelons),1),:).';
wflonlinelon=wflonlinelons(ones(numel(wflonlinelats),1),:);

% convert wireframe to xyz
[wfx1,wfy1,wfz1]=geographic2xyz(wflatlinelat,wflatlinelon,zeros(size(wflatlinelat)));
[wfx2,wfy2,wfz2]=geographic2xyz(wflonlinelat,wflonlinelon,zeros(size(wflonlinelat)));

% plot wireframe
h=figure;
set(0,'DefaultAxesColorOrder',[0 0 0])
plot3(wfx1,wfy1,wfz1);
hold on
plot3(wfx2,wfy2,wfz2);

% how many mlocate files to read in?
zmin=nan; zmax=nan;
for n=1:length(files)
    % read in info
    [lat,lon,dep]=textread(files(n).name,'%10.4f%10.4f%10.4f',1,'headerlines',1);
    [years,months,days,hours,minutes,mag,...
        latdev,londev,depdev,...
        xx,xy,xz,xdist,...
        yx,yy,yz,ydist,...
        zx,zy,zz,zdist]=textread(files(n).name,['%4d-%2d-%2d-%2d:%2d-%5.2f\n'...
                                           '%8.2f%8.2f%8.2f\n'...
                                           '   %10.5f%10.5f%10.5f%10.5f\n'...
                                           '   %10.5f%10.5f%10.5f%10.5f\n'...
                                           '   %10.5f%10.5f%10.5f%10.5f\n'],...
                                           'headerlines',7);
    
    % loop through each ellipsoid
    nellips=length(xx);
    for i=1:nellips
        % rotation matrix (to orient ellipsoid)
        v=inv([xx(i) xy(i) xz(i);
               yx(i) yy(i) yz(i);
               zx(i) zy(i) zz(i)]);
        
        % ellipsoid (unoriented)
        [x,y,z]=ellipsoid(0,0,0,xdist(i),ydist(i),zdist(i));
        
        % rotate ellipsoid
        temp=[x(:) y(:) z(:)]*v;
        
        % translate ellipsoid (to centroid reference frame)
        x(:)=temp(:,1)+latdev(i);
        y(:)=temp(:,2)+londev(i);
        z(:)=temp(:,3)+depdev(i);
        clear temp
        
        % translate ellipsoid (to WGS-84 ellipsoid Earth reference frame)
        % km north and east ==> distance and azimuth (this could cause error)
        % lat, lon, distance, bearing ==> lat, lon
        [x(:),y(:)]=vincentyfwd(lat(ones(441,1)),lon(ones(441,1)),...
            sqrt(x(:).^2+y(:).^2),-atan2(x(:),y(:))*180/pi+90);
        z=dep+z;
        
        % check for minimum/maximum
        zmin=min([zmin; z(:)]);
        zmax=max([zmax; z(:)]);
        
        % translate geographic position to cartesian
        [x(:),y(:),z(:)]=geographic2xyz(x(:),y(:),z(:));
        
        % draw ellipsoid (cycle ellipsoid coloring every 100km)
        figure(h);
        surf(x,y,z,'facecolor',cmap(1+mod(round(cint*(dep+depdev(i))/repeatint),cint),:),'edgecolor','none');
    end
end

% setup colormap used
zmin=round(zmin); zmax=round(zmax);
zminr=fix(zmin/repeatint); zmaxr=fix(zmax/repeatint);
fmini=1+mod(round(cint*zmin/repeatint),cint);
fmaxi=1+mod(round(cint*zmax/repeatint),cint);
cmap=[cmap(fmini:end,:); repmat(cmap,[zmaxr-zminr 1]); cmap(1:fmaxi,:)];
colormap(cmap)

% universal plotting parameters
figure(h)
colorbar
box off
axis vis3d
axis equal
axis off
lighting phong
material shiny
light('Position',[ 1 -1  0],'Style','infinite')
light('Position',[ 1  1  0],'Style','infinite')
light('Position',[-1  0  0],'Style','infinite')
hold off
rotate3d

end


