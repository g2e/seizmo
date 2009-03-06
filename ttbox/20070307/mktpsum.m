function ttime=mktpsum(p,v,z,zmin,zmax,novertex);
% mktpsum.........integrates traveltime over several layers IN FLAT EARTH
%
% call: ttime=mktpsum(p,v,z,zmin,zmax,novertex);
%
%             p: ray parameter [s/km]
%             v: velocities at depths given in z [km/s]
%             z: depths, at which velocities are given [km]
%          zmin: depth at which ray starts [km]
%          zmax: depth at which ray ends [km]
%      novertex: set novertex=1 if you use MKTPSUM to calculate
%                a focal depth correction.
%
% all parameters have to be in FLAT EARTH coordinates/dimensions!
%
% result: ttime: time needed to travel from depth zmin to depth zmax [s]
%                NaN in case of any kind of problem
%                inf if ray goes through layer with zero velocity
%                (as S does in earth's outer core)
%                -1 for layers above ray segment (used especially for K-leg of SKS and the like)
%
% This function computes the traveltime of a ray from depth ZMIN
% to depth ZMAX. Nothing else. You have to care for changes from P to S,
% for handling of reflections, for doubliung the time after the vertex
% etc.
%
% uses MKTP to compute times spent within a layer.
%
% Martin Knapmeyer, 24.04.2002

% init result
ttime=NaN;
res=[];


%%% va and z have to be row vectors
v=v(:)';
z=z(:)';

%%% number of layers
anz=length(v)-1; 


%%% go through all layers which are not below ZMAX
done=0;
cnt=1;
while ~done
   if z(cnt)>=zmax
      % layer is below deepest allowed point of ray - and we're done.
      %disp(['MKTPSUM: layer too deep - abort.']);
      done=1;
   else
      if (zmin>z(cnt))&(zmax>z(cnt))
         %disp(['MKTPSUM: layer too shallow - return -1 (cnt=' int2str(cnt) ')']);
         res(cnt)=-1;
      else
         % ray goes through this layer
         layerv=[v(cnt) v(cnt+1)]; % velocities at top and bottom of current layer
         layerz=[z(cnt) z(cnt+1)]; % depth of top and bottom of current layer
         if novertex==1
            res(cnt)=mktp(p,layerv,layerz,zmin,zmax,1);
         else
            res(cnt)=mktp(p,layerv,layerz,zmin,zmax);
         end; % if novertex
%figure(2); hold on; plot(layerz,[1 1]*res(cnt),'.-'); hold off; drawnow
      end; % if (zmin>z(cnt))&(zmax>z(cnt))
   end; % if layerz(1)
   cnt=cnt+1; % number of next layer
end; % while ~done
%figure(2); hold on; plot(z(1:(cnt-2)),res,'ro');


%%% if RES contains NaN elemtents, the ray is invalid! MK26102006
if sum(isnan(res))>0
   res=NaN;
   return;
end; % if sum(isnan(res))>0

%%% find elements of res that contribute to the travel time:
%%% these are all elements before the first NaN element.
%%% and all elements that are not -1
indies=find(res>-1); % non-(-1) elements only mk 14052002
if ~isempty(indies)
   res=res(indies);
else
   %%% what if there are only -1 values?
   %%% return a dummy ray of length zero to avoid indexing problems in MKX4P
   if ~isempty(res)
      if isempty(find(res~=-1))
         res=0; %% 0: layer does not contribute anything!
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
ttime=sum(res(indies));  % sum of times of all layers

