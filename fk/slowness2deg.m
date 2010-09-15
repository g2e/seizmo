function [deg]=slowness2deg(ph,slow)
%SLOWNESS2DEG    Converts horizontal slowness (s/deg) to distance (deg)
%
%    Usage:    deg=slowness2deg(phase,slowness)
%
%    Description:
%     DEG=SLOWNESS2DEG(PHASE,SLOWNESS) converts a horizontal slowness given
%     in seconds per degree for a specific seismic phase to the originating
%     distance in degrees using the travel time model AK135.  Allowed
%     seismic phase strings are 'P', 'PP', 'PKPab', 'PKPbc', & 'PKPdf'.
%     SLOWNESS may be an array.
%
%    Notes:
%
%    Examples:
%     % Where and when is the slowness of a P-wave 6 sec/deg?
%     deg=slowness2deg('P',6)
%     tauptime('ph','P','mod','ak135','deg',deg)
%
%    See also: TAUPCURVE, TAUPTIME

%     Version History:
%        Aug. 30, 2010 - initial version
%        Sep. 13, 2010 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 13, 2010 at 10:35 GMT

% todo

% check nargin
error(nargchk(2,2,nargin));

% check phase
if(~isstring(ph))
    error('seizmo:slowness2deg:badInput',...
        'PHASE must be a string!');
end

% check slowness
if(~isreal(slow))
    error('seizmo:slowness2deg:badInput',...
        'SLOWNESS must be a real valued array!');
end

% load model
dp=load('ak135slow');

% act by slowness type
switch ph
    case 'P'
        deg=interp1(dp.P(:,2),dp.P(:,1),slow);
    case 'PP'
        deg=interp1(dp.PP(:,2),dp.PP(:,1),slow);
    case 'PKPab'
        deg=interp1(dp.PKPab(:,2),dp.PKPab(:,1),slow);
    case 'PKPbc'
        deg=interp1(dp.PKPbc(:,2),dp.PKPbc(:,1),slow);
    case {'PKIKP' 'PKPdf'}
        deg=interp1(dp.PKPdf(:,2),dp.PKPdf(:,1),slow);
    otherwise
        error('seizmo:slowness2deg:badPhase',...
            'Unknown phase: %s',ph);
end

end
