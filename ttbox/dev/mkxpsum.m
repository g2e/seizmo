function [dist,segx,segz]=mkxpsum(p,v,z,zmin,zmax,novertex);
% mkxpsum.........integrates distance over several layers IN FLAT EARTH
%
% call: [dist,segx,segz]=mkxpsum(p,v,z,zmin,zmax,novertex);
%
%             p: ray parameter [s/km]
%             v: velocities at depths given in z [km/s]
%             z: depths, at which velocities are given [km]
%          zmin: depth at which ray starts [km]
%          zmax: depth at which ray ends [km]
%      novertex: set novertex=1 if you use MKXPSUM to calculate
%                a focal depth correction.
%
% all parameters have to be in FLAT EARTH coordinates/dimensions!
%
% result: dist:  horizontal distance covered while traveling from depth zmin to depth zmax [km]
%                NaN in case of any kind of problem
%                inf if ray goes through a layer with zero velocity.
%                (as S does in the earth's outer core)
%                -1 for layers above ray segment (used espcially for K-leg of SKS and the like)
%         segx:  cumulative sum of horizontal distance segments traveled in
%                each layer [km] IN FLAT EARTH
%                together with SEGZ, this gives the ray geometry
%         segz:  vertical distance segments traveled in each layer [km]
%                these are more or less the layer interface depths.
%                first element is ZMIN, last element is ZMAX, elements inbetween
%                are the interface depths between ZMIN and ZMAX.
%
% This function computes the horizontal distance covered of by ray that
% travels from depth ZMIN to depth ZMAX. Nothing else. You have to care 
% for changes from P to S, for handling of reflections, for doubliung
% the time after the vertex etc.
%
% uses MKXP to compute distance travelled within a layer.
%
% Martin Knapmeyer, 30.04.2002, 13.09.2004, 16.02.2006, 13.07.2006,
%                   11.10.2006
%
% BUGS: an error occurs if a vertex is in the lowermost layer, due to sloppy array
%       indexing. This only happens if yur ray parameter is too small.


% init result
dist=NaN;
segx=NaN;
segz=NaN;
res=[];
resz=[];

%%% va and z have to be row vectors
v=v(:)';
z=z(:)';

%%% largest ray parameter for which zero epocentral distance is allowed
pepsilon=1e-10; % s/km

%%% number of layers
anz=length(v)-1; 

%%% go through all layers which are not below ZMAX
done=0;
cnt=1;
while ~done
   %disp(['MKXPSUM: cnt=' int2str(cnt) ' ' num2str(z(cnt)) ' ' num2str(zmin) ' ' num2str(zmax)]);
   %disp([int2str((z(cnt)>=zmax)) ' ' int2str(((zmin>z(cnt))&(zmax>z(cnt))))]);
   if (z(cnt)>=zmax)
      %%% layer is below deepest allowed point of ray - and we're done.
      %disp(['MKXPSUM: layer too deep - cancel loop (cnt=' int2str(cnt) ')']);
      done=1;
      res(cnt)=-1;
      resz(cnt)=-1; % MK13092004
   else
      if (zmin>z(cnt))&(zmax>z(cnt))
         %disp(['MKXPSUM: layer too shallow - return -1 (cnt=' int2str(cnt) ')']);
         res(cnt)=-1;
         resz(cnt)=-1; % MK13092004
      else
         % ray goes through this layer
         layerv=[v(cnt) v(cnt+1)]; % velocities at top and bottom of current layer
         layerz=[z(cnt) z(cnt+1)]; % depth of top and bottom of current layer
         if novertex==1
            res(cnt)=mkxp(p,layerv,layerz,zmin,zmax,1);
            %resz(cnt)=max(layerz); % MK13092004
            resz(cnt)=min(max(layerz),zmax); % MK15022006
         else
            res(cnt)=mkxp(p,layerv,layerz,zmin,zmax);
            %resz(cnt)=max(layerz); % MK13092004
            resz(cnt)=min(max(layerz),zmax); % MK15022006
            %disp(['MKXPSUM: cnt=' int2str(cnt) ' res=' num2str(res(cnt)) ' ' num2str(layerz) ' ' num2str(zmin) ' ' num2str(zmax)]);
         end; % if novertex
      end; % if (zmin>z(cnt))&(zmax>z(cnt))
   end; % if layerz(1)
   cnt=cnt+1; % number of next layer
end; % while ~done

%%% find elements of res that contribute to the distance:
%%% these are all positive elements before the first NaN element.
%%% and all elements that are not -1.
indies=find(res~=-1); % non-(-1) elements only mk 14052002
%disp(['size: ' int2str(size(indies)) ' isempty: ' int2str(isempty(indies))]);
if ~isempty(indies)
   %disp('MKXPSUM: this is RES:')
   %dmy=res;
   %disp('MKXPSUM: this is INDIES:');
   %disp(int2str(indies))
   %disp(num2str(dmy(indies)))
   res=res(indies);
   resz=resz(indies);
else
   %%% what if there are only -1 values?
   %%% return a dummy ray of length zero to avoid indexing problems in MKX4P
   if ~isempty(res)
      if isempty(find(res~=-1))
         res=[0 0]; %% 0: layer does not contribute anything!
         resz=[0 0]; % MK13092004
      end; % if isempty(find())
   end; % if ~isempty(res)
end; % if ~isempty(indies)
indies=find(isnan(res));
if ~isempty(indies) 
   indies=1:(indies(1)-1);
else
   indies=1:length(res);
end; % if ~isnan(indies) else



%%% return result
if isempty(indies)
   %%% this means that there is no non-NaN element in RES
   dist=NaN;
   segx=NaN;
   segz=NaN;
else
   dist=sum(res(indies));  % sum of distances of all layers
   if sum(res)==0
      %%% ray does not cover any horizontal distance. This is allowed only
      %%% if it is very near to vertical. MK11102006
      if p<pepsilon
         segx=([0 res(indies)]); % made of individual horizontal distances for each layer
         segz=[zmin resz(indies)]; % all z between zmin and zmax MK11102006
      else
          segx=res;
          segz=res;
      end; % if p<pepsilon
   else
      segx=([0 res(indies)]); % made of individual horizontal distances for each layer
      %segz=[zmin z(find((z>zmin)&(z<zmax))) zmax]; % all z between zmin and zmax MK30042002
      segz=[zmin resz(indies)]; % all z between zmin and zmax MK13092004
   end; % if prod(res==[0 0])~=0
end; % if isempty(indies)

%%% if the source is in a Low Velocity Zone, rays for certain ray
%%% parameters cannot leave the LVZ and the focal depth correction will
%%% fail to construct a surface-to-focus ray. In this case, the ray
%%% constructed here and defined by SEGX and SEGZ does not reach the focal
%%% depth and is thus invalid.
%%% Note that it is still possible to compute (meaningless) traveltimes, the
%%% recognition of such invalid rays can be done via ray geometry only!
%%% MK13072006
if max(segz)<zmax
   %%% ray does not reach desired depth, so throw it away!
   dist=NaN;
   segx=NaN;
   segz=NaN;
end;

% figure(3);
% subset=length(segx);
% hold on
% plot(p,segx(end),'.'); %sqrt((segx(subset).*segx(subset)+segz(subset).*segz(subset))),'.');
% hold off
% figure(1);