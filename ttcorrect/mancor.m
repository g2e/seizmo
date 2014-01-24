function [corr,tt1d,tt3d]=mancor(paths,mod3d)
%MANCOR    Returns mantle travel time corrections for a set of raypaths
%
%    Usage:    corr=mancor(paths,mod3d)
%              [corr,tt1d,tt3d]=mancor(...)
%
%    Description:
%     CORR=MANCOR(PATHS,MOD3D) calculates travel time corrections for phase
%     raypaths given by PATHS through the 3D mantle model MOD3D.  MOD3D
%     should be one of the models listed by AVAILABLE_3DMODELS.  PATHS are
%     expected to follow the format of the struct output from TAUPPATH.
%     Please note that paths are not checked for segments that leave the
%     mantle and so you should "prepare" the paths by running them through
%     CRUSTLESS_RAYPATHS & TRIM_DEPTHS_RAYPATHS to remove the crust
%     segments and core segments.  CORR is in seconds and gives
%     TT3D=TT1D+CORR.  Combining MANCOR with ELLCOR and CRUCOR will provide
%     a more complete 3D correction.
%
%     [CORR,TT1D,TT3D]=MANCOR(...) returns the total travel time for the
%     paths in the 3D model (TT3D) and in the 1D model (TT1D).
%
%    Notes:
%     - Latitudes in PATH are assumed to be geocentric.  You should always
%       pass geocentric latitudes to TAUPPATH if you use it.  The function
%       GETRAYPATHS does this for you.
%
%    Examples:
%     % An example of mantle correction workflow:
%     paths=tauppath('ev',[31.5 140.07],'st',[2.389 9.834],...
%                    'ph','ttp+','mod','prem');
%     cmb=2891; % prem based cmb depth
%     gpaths=crustless_raypaths(trim_depths_raypaths(paths,[0 cmb]));
%     corr=mancor(gpaths,'hmsl06p')
%
%    See also: AVAILABLE_3DMODELS, GETRAYPATHS, TRIM_DEPTHS_RAYPATHS,
%              CRUSTLESS_RAYPATHS, EXTRACT_UPSWING_RAYPATHS, PLOTRAYPATHS,
%              TAUPPATH, INSERT_DEPTHS_IN_RAYPATHS

%     Version History:
%        June  2, 2010 - initial version
%        June  3, 2010 - Fermat's principle formulation
%        June  4, 2010 - drop amp input
%        Jan. 14, 2011 - improved verbose message
%        Feb. 27, 2012 - update for tauppath changes
%        Aug.  6, 2012 - doc update
%        Jan. 23, 2014 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 02:45 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check inputs
sm3d=size(mod3d);
test=tauppath('ph','P','ev',[0 0],'st',[0 10]);
if(~isstruct(paths) || any(~ismember(fieldnames(paths),fieldnames(test))))
    error('seizmo:mancor:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember(fieldnames(paths(1).path),fieldnames(test(1).path))))
    error('seizmo:mancor:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember({'latitude' 'longitude'},fieldnames(paths(1).path))))
    error('seizmo:mancor:badStruct',...
        'Latitude & Longitude are required path fields!');
elseif(~ischar(mod3d) || numel(sm3d)~=2 || sm3d(1)~=1)
    error('seizmo:mancor:badInput',...
        'MOD3D must be a string!');
end

% number of ray paths
nrp=numel(paths);

% verbosity
verbose=seizmoverbose;
if(verbose)
    disp(['Getting Mantle Correction(s) For ' upper(mod3d)]);
    print_time_left(0,nrp);
end

% loop over each raypath
corr=zeros(size(paths)); tt1d=corr;
for i=1:nrp
    % look out for nan depths
    npts=numel(paths(i).path.depth);
    nn=~isnan(paths(i).path.depth);
    gseg=find(nn(1:npts-1) & nn(2:npts))';
    
    % get dlnv at each point
    dlnv=nan(npts,1);
    dlnv(nn)=mantledv(mod3d,paths(i).path.latitude(nn),...
        paths(i).path.longitude(nn),paths(i).path.depth(nn));
    
    % get travel time for each good segment
    tt=paths(i).path.time(gseg+1)-paths(i).path.time(gseg);
    
    % get correction
    % - using formulation based on Fermat's principle
    % - average perturbation effect of dlnv on each end of every segment
    corr(i)=-sum(tt.*(dlnv(gseg)+dlnv(gseg+1)))/2;
    
    % 1D total travel time (if desired)
    if(nargout>1); tt1d(i)=sum(tt); end
    
    % detail message
    if(verbose); print_time_left(i,nrp); end
end

% 3D total travel time (if desired)
if(nargout==3); tt3d=tt1d+corr; end

end
