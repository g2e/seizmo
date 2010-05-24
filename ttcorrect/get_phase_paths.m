function [paths]=get_phase_paths(evla,evlo,evdp,stla,stlo,mod,phase)
%GET_PHASE_PATHS    Returns seismic phase paths for a set of stations
%
%    Usage:    paths=get_phase_paths(evla,evlo,evdp,stla,stlo,mod,phase)
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also: SPLIT_PHASE_PATHS, GET_UPSWING_PATHS

%     Version History:
%        May  21, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  21, 2010 at 12:25 GMT

% todo:
% - check inputs
% - output struct

% check nargin
error(nargchk(7,7,nargin));

% check inputs


% expand scalars
[evla,evlo,evdp,stla,stlo,mod,phase]=expandscalars(...
    evla,evlo,evdp,stla,stlo,mod,phase);
nph=numel(evla);

% geographic to geocentric lat
evla=geographic2geocentriclat(evla);
stla=geographic2geocentriclat(stla);

% loop over each phase
for i=1:nph
    
end

end
