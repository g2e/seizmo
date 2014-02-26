function [cu]=read_cu_mod(file)
%READ_CU_MOD    Reads CU_SRT1.0 or CU_SDT1.0 style model
%
%    Usage:    model=read_cu_mod(file)
%
%    Description:
%     MODEL=READ_CU_MOD(FILE) reads in the Colorado surface wave tomography
%     model CU_SRT1.0 or CU_SDT1.0.  The output is a struct MODEL
%     containing several fields:
%      .name                    - model name
%      .reference               - reference article for model
%      .lat                     - latitudes in grid
%      .lon                     - longitudes in grid
%      .depth                   - depths in grid
%      .vs  .vsmin  .vsmax      - isotropic velocities & model limits
%      .vsv .vsvmin .vsvmax     - anisotropic vert velo & model limits
%      .vsh .vshmin .vshmax     - anisotropic horz velo & model limits
%
%    Notes:
%     - The Colorado CU_SRT1.0 & CU_SDT1.0 models are quite large.
%       Unfortunately reading the model requires ~6 gigabytes of memory.
%       I would suggest splitting the file into sections (I split it into 4
%       pieces of 23, 22, 22, & 22 latitude stripes).  Then concatenate the
%       models together.  
%
%    Examples:
%     % Graphically select model file, read it, and then plot a few
%     % layers to see if things look right:
%     model=read_cu_mod;
%     figure; imagesc(model.vs(:,:,3));
%     figure; imagesc(model.vsh(:,:,3));
%     figure; imagesc(model.vsv(:,:,3));
%
%    See also: READCRUST2, READCRUST10

%     Version History:
%        Jan. 21, 2011 - initial version
%        Jan. 26, 2014 - abs path exist fix
%        Feb.  9, 2014 - use readtxt
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2014 at 10:35 GMT

% todo

% check nargin
error(nargchk(0,1,nargin));

% read in file
if(nargin<1); file=[]; end
filterspec={
    '*.mod;*.MOD' 'MOD Files (*.mod,*.MOD)';
    '*.txt;*.TXT' 'TXT Files (*.txt,*.TXT)';
    '*.*' 'All Files (*.*)'};
[v,file]=readtxt(file,filterspec);
v=getwords(v,sprintf('\n'));

% initialize struct
cu=struct('name',[],'reference',[],...
    'lat',[],'lon',[],'depth',[],...
    'vs',[],'vsv',[],'vsh',[],...
    'vsmin',[],'vsvmin',[],'vshmin',[],...
    'vsmax',[],'vsvmax',[],'vshmax',[]);

% model name and reference
[path,name]=fileparts(file);
cu.name=name;
cu.reference='Ritzwoller et al 2002';

% parameters (change if the model parameters change)
nlat=89;
nlon=180;
ndep=100;

% check nlines
if(numel(v)~=nlat*nlon*(ndep+1))
    error('seizmo:read_cu_mod:badCUMOD',...
        'CU MOD file is not formatted correctly!');
end

% extract lat/lon positions
latlon=single(str2num(char(v(1:101:end))));
%row=(90-latlon(:,1))/2; % b/c we know (change if model parameters change)
%col=latlon(:,2)/2+1;    % b/c we know (change if model parameters change)
v(1:101:end)=[];
cu.lat=flipud(unique(latlon(:,1)));
cu.lon=unique(latlon(:,2)).';

% get data values
v=single(str2num(char(v)));

% extract depths
cu.depth=v(1:ndep,1);
v(:,1)=[];

% this is the super fast way (change if model parameters change)
% [depth lon lat] => [lat lon depth]
cu.vs=flipdim(permute(reshape(v(:,1),[ndep nlon nlat]),[3 2 1]),1);
v(:,1)=[];
cu.vsv=flipdim(permute(reshape(v(:,1),[ndep nlon nlat]),[3 2 1]),1);
v(:,1)=[];
cu.vsh=flipdim(permute(reshape(v(:,1),[ndep nlon nlat]),[3 2 1]),1);
v(:,1)=[];
cu.vsmin=flipdim(permute(reshape(v(:,1),[ndep nlon nlat]),[3 2 1]),1);
v(:,1)=[];
cu.vsvmin=flipdim(permute(reshape(v(:,1),[ndep nlon nlat]),[3 2 1]),1);
v(:,1)=[];
cu.vshmin=flipdim(permute(reshape(v(:,1),[ndep nlon nlat]),[3 2 1]),1);
v(:,1)=[];
cu.vsmax=flipdim(permute(reshape(v(:,1),[ndep nlon nlat]),[3 2 1]),1);
v(:,1)=[];
cu.vsvmax=flipdim(permute(reshape(v(:,1),[ndep nlon nlat]),[3 2 1]),1);
v(:,1)=[];
cu.vshmax=flipdim(permute(reshape(v(:,1),[ndep nlon nlat]),[3 2 1]),1);

end
