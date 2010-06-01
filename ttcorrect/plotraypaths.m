function []=plotraypaths(paths)
%PLOTRAYPATHS    3D plot of the Globe with raypaths

% check nargin
error(nargchk(1,1,nargin));

% check raypath struct
test=tauppath('ph','P','deg',10);
if(~isstruct(paths) || any(~ismember(fieldnames(paths),fieldnames(test))))
    error('seizmo:plotraypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember(fieldnames(paths(1).path),fieldnames(test(1).path))))
    error('seizmo:plotraypaths:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
end

% number of raypaths
nrp=numel(paths);

% constants
Re=6371;

% colormap
cmap=hsv(nrp);
figure;

% loop over each raypath
rph=nan(1,nrp);
for i=1:nrp
    % convert lat/lon/depth to x/y/z
    [x,y,z]=geocentric2xyz(paths(i).path.latitude,...
        paths(i).path.longitude,paths(i).path.depth,Re);
    
    % plot path
    rph(i)=plot3(-x,y,z,'color',cmap(i,:));
    if(i==1); hold on; end
end

% tranparent globe with continents
%cla reset;
load topo;
[x y z] = sphere(45);
s = surface(Re*x,Re*y,Re*z,'facecolor','texturemap','cdata',topo);
set(s,'edgecolor','none','facealpha','texture','alphadata',topo);
set(s,'backfacelighting','unlit');
colormap(topomap1);
alpha('direct');
alphamap([.1;1])
axis off vis3d;
campos([2 13 10]);
camlight;
lighting gouraud;

% need to plot cmb icb

% option to plot others?

% alpha settings?
% color settings?



end
