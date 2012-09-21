function [deg]=slowness2deg(ph,slow)
%SLOWNESS2DEG    Converts horizontal slowness (s/deg) to distance (deg)
%
%    Usage:    deg=slowness2deg(phase,slowness)
%
%    Description:
%     DEG=SLOWNESS2DEG(PHASE,SLOWNESS) converts a horizontal slowness given
%     in seconds per degree for a specific seismic phase to the originating
%     distance in degrees using the travel time model AK135.  Allowed
%     seismic phase strings are 'P', 'PP', 'PKPab', 'PKPbc', 'PKPdf'
%     (aka PKIKP), 'PKiKP', 'S', 'SS', 'SKS', 'SKiKS', and 'SKIKS'.  NOTE
%     THAT THE S-WAVE CORE-PHASES ARE FROM IASP91.  This is because AK135
%     SKS is a bit strange.  SLOWNESS may be an array.
%
%    Notes:
%     - A slowness will give a specific distance but a distance may not
%       give a specific slowness due to triplications.  Please note that
%       the companion function DEG2SLOWNESS only works at lowermantle
%       distances for P, PP, S, & SS.  Use TAUPTIME to get multiple
%       slownesses for a phase at a specific distance.
%
%    Examples:
%     % Where and when is the slowness of a P-wave 6 sec/deg?
%     deg=slowness2deg('P',6)
%     tauptime('ph','P','mod','ak135','deg',deg)
%
%    See also: DEG2SLOWNESS, TAUPCURVE, TAUPTIME

%     Version History:
%        Aug. 30, 2010 - initial version
%        Sep. 13, 2010 - doc update
%        Sep. 21, 2010 - several more phases
%        Apr.  3, 2012 - minor doc update
%        May  18, 2012 - forgot to mention PKiKP in documentation
%        Sep.  6, 2012 - user *lowermantle for truncated curves
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  6, 2012 at 10:35 GMT

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

% act by phase
switch ph
    case 'P'
        deg=interp1(dp.P(:,2),dp.P(:,1),slow);
    case 'PP'
        deg=interp1(dp.PP(:,2),dp.PP(:,1),slow);
    case 'Plowermantle'
        deg=interp1(dp.Plowermantle(:,2),dp.Plowermantle(:,1),slow);
    case 'PPlowermantle'
        deg=interp1(dp.PPlowermantle(:,2),dp.PPlowermantle(:,1),slow);
    case 'PKPab'
        deg=interp1(dp.PKPab(:,2),dp.PKPab(:,1),slow);
    case 'PKPbc'
        deg=interp1(dp.PKPbc(:,2),dp.PKPbc(:,1),slow);
    case 'PKiKP'
        deg=interp1(dp.PKiKP(:,2),dp.PKiKP(:,1),slow);
    case {'PKIKP' 'PKPdf'}
        deg=interp1(dp.PKPdf(:,2),dp.PKPdf(:,1),slow);
    case 'S'
        deg=interp1(dp.S(:,2),dp.S(:,1),slow);
    case 'SS'
        deg=interp1(dp.SS(:,2),dp.SS(:,1),slow);
    case 'Slowermantle'
        deg=interp1(dp.Slowermantle(:,2),dp.Slowermantle(:,1),slow);
    case 'SSlowermantle'
        deg=interp1(dp.SSlowermantle(:,2),dp.SSlowermantle(:,1),slow);
    case 'SKS'
        % iasp91
        deg=interp1(dp.SKS(:,2),dp.SKS(:,1),slow);
    case 'SKiKS'
        % iasp91
        deg=interp1(dp.SKiKS(:,2),dp.SKiKS(:,1),slow);
    case 'SKIKS'
        % iasp91
        deg=interp1(dp.SKIKS(:,2),dp.SKIKS(:,1),slow);
    otherwise
        error('seizmo:slowness2deg:badPhase',...
            'Unknown phase: %s',ph);
end

end
