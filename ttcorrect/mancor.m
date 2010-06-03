function [corr,tt3d,tt1d]=mancor(paths,mod3d,amp)
%MANCOR    Returns mantle travel time corrections for a set of raypaths
%
%    Usage:    corr=mancor(paths,mod3d)
%              corr=mancor(paths,mod3d,amp)
%              [corr,tt3d,tt1d]=mancor(...)
%
%    Description: CORR=MANCOR(PATHS,MOD3D) calculates travel time
%     corrections for phase raypaths given by PATHS through the 3D mantle
%     model MOD3D.  MOD3D should be one of the models listed by
%     AVAILABLE_3DMODELS.  PATHS are expected to follow the format of the
%     struct output from TAUPPATH.  Please note that paths are not checked
%     for segments that leave the mantle and so you should "prepare" the
%     paths by running them through CRUST2LESS_RAYPATHS &
%     TRIM_DEPTHS_RAYPATHS to remove the crust segments and core segments.
%     CORR is in seconds and gives TT3D=TT1D+CORR.  Combining MANCOR with
%     ELLCOR and CRUCOR will provide a more complete 3D correction.
%
%     CORR=MANCOR(PATHS,MOD3D,AMP) amplifies MOD3D by AMP when calculating
%     the travel time corrections.  Setting AMP to 2 will double the
%     strength of the velocity deviations in the model while setting amp to
%     1/3 would reduce the amplitude of the anomalies by 1/3.
%
%     [CORR,TT3D,TT1D]=MANCOR(...) returns the total travel time for the
%     paths in the 3D model (TT3D) and in the 1D model (TT1D).
%
%    Notes:
%
%    Examples:
%     An example mantle correction workflow:
%      paths=tauppath('ev',[31.5 140.07],'st',[2.389 9.834],...
%                     'ph','ttp+','mod','prem');
%      cmb=2891; % prem based cmb depth
%      gpaths=crust2less_raypaths(trim_depths_raypaths(paths,[0 cmb]));
%      corr=mancor(gpaths,'hmsl06p')
%
%    See also: AVAILABLE_3DMODELS, GETRAYPATHS, TRIM_DEPTHS_RAYPATHS,
%              CRUST2LESS_RAYPATHS, GET_UPSWING_RAYPATHS

%     Version History:
%        June  2, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  2, 2010 at 23:15 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% check inputs
sm3d=size(mod3d);
test=tauppath('ph','P','deg',10);
if(nargin==2 || isempty(amp)); amp=1; end
if(~isstruct(paths) || any(~ismember(fieldnames(paths),fieldnames(test))))
    error('seizmo:mancor:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(any(~ismember(fieldnames(paths(1).path),fieldnames(test(1).path))))
    error('seizmo:mancor:badStruct',...
        'PATHS does not appear to be a valid raypath struct!');
elseif(~ischar(mod3d) || numel(sm3d)~=2 || sm3d(1)~=1)
    error('seizmo:mancor:badInput',...
        'MOD3D must be a string!');
elseif(~isreal(amp) || ~isscalar(amp) || amp<=0)
    error('seizmo:mancor:badInput',...
        'AMP must be a positive real-valued scalar!');
end

% number of ray paths
nrp=numel(paths);

% verbosity
verbose=seizmoverbose;
if(verbose)
    disp('Getting Mantle Correction(s)');
    print_time_left(0,nrp);
end

% loop over each raypath
corr=zeros(size(paths)); tt3d=corr; tt1d=corr;
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
    
    % get total travel times
    % - 3d uses average of travel times based on starting and ending
    %   dlnv for each segment to account for velocity gradients
    tt1d(i)=sum(tt);
    tt3d(i)=sum(tt./(1+amp*dlnv(gseg))+tt./(1+amp*dlnv(gseg+1)))/2;
    
    % detail message
    if(verbose); print_time_left(i,nrp); end
    
    % debug
    %disp(['PHASE: ' paths(i).phase ' TAUP: ' num2str(paths(i).time) ...
    %      's TT1D: ' num2str(tt1d(i)) 's TT3D: ' num2str(tt3d(i)) 's']);
end

% corrections
corr=tt3d-tt1d;

end
