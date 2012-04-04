function [model]=flatten_1dmodel(model)
%FLATTEN_1DMODEL    Flattens a 1D Earth model
%
%    Usage:    model=flatten_1dmodel(model)
%
%    Description:
%     MODEL=FLATTEN_1DMODEL(MODEL) takes a spherical 1D Earth model and
%     performs an Earth-flattening transformation on it (converts it into a
%     half-space).  Depth points at the center of the Earth will be at
%     infinite depth.  For details see reference in Notes.
%
%    Notes:
%     - See Kennett 1983, Seismic Wave Propagation in a Stratified Medium,
%       page 20 for details on Earth flattening.
%
%    Examples:
%     % Compare a flattened Earth model to a non-flat one:
%     mod=prem;
%     fmod=flatten_1dmodel(mod);
%     figure; plot(mod.depth(1:50),mod.vp(1:50),...
%                 fmod.depth(1:50),fmod.vp(1:50))
%
%    See also: PREM, AK135, IASP91

%     Version History:
%        May  24, 2010 - initial version
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 14:45 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check model
chk1dmodel(model);

% already flattened?
if(model.flattened)
    warning('seizmo:flatten_1dmodel:alreadyFlat',...
        'MODEL is already flattened!');
    return;
end

% flatten model
Re=6371;
r=(Re-model.depth)/Re;
model.rho=model.rho.*r;
model.vp=model.vp./r;
model.vs=model.vs./r;
model.depth=-Re*log(r);

end
