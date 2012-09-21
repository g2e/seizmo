function [slow]=deg2slowness(ph,deg)
%DEG2SLOWNESS    Converts distance (deg) to horizontal slowness (s/deg)
%
%    Usage:    slow=deg2slowness(phase,deg)
%
%    Description:
%     SLOW=DEG2SLOWNESS(PHASE,DEG) converts degree distance DEG to
%     horizontal slowness given in seconds per degree for a specific
%     seismic phase using the travel time model AK135.  Allowed seismic
%     phase strings are 'P', 'PP', 'PKPab', 'PKPbc', 'PKPdf' (aka PKIKP),
%     'S', 'SS', 'SKS', 'SKiKS', and 'SKIKS'.  NOTE THAT THE S-WAVE CORE
%     PHASES ARE FROM IASP91.  This is because AK135 SKS has some funkyness
%     going on.  SLOWNESS may be an array.
%
%    Notes:
%     - Only the lowermantle portion of the P,PP,S,SS phases can be used
%       because of triplications (and thus 3+ slownesses for some
%       distances).  If you really need all slownesses for a given distance
%       then use TAUPTIME and/or TAUPPATH.
%
%    Examples:
%     % What is the slowness of a P-wave at 60 deg?
%     slow=deg2slowness('P',60)
%     tauptime('ph','P','mod','ak135','deg',60)
%
%    See also: SLOWNESS2DEG, TAUPCURVE, TAUPTIME

%     Version History:
%        Apr.  3, 2012 - initial version
%        Sep.  6, 2012 - use *lowermantle for truncated curves but skip
%                        actually using them due to interpolation problem
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  6, 2012 at 10:35 GMT

% todo

% check nargin
error(nargchk(2,2,nargin));

% check phase
if(~isstring(ph))
    error('seizmo:deg2slowness:badInput',...
        'PHASE must be a string!');
end

% check distance
if(~isreal(deg))
    error('seizmo:deg2slowness:badInput',...
        'DEG must be a real valued array!');
end

% load model
dp=load('ak135slow');

% act by phase
switch ph
    case {'P' 'Plowermantle'}
        slow=interp1(dp.Plowermantle(:,1),dp.Plowermantle(:,2),deg);
    case {'PP' 'PPlowermantle'}
        slow=interp1(dp.PPlowermantle(:,1),dp.PPlowermantle(:,2),deg);
    case 'PKPab'
        slow=interp1(dp.PKPab(:,1),dp.PKPab(:,2),deg);
    case 'PKPbc'
        slow=interp1(dp.PKPbc(:,1),dp.PKPbc(:,2),deg);
    case 'PKiKP'
        slow=interp1(dp.PKiKP(:,1),dp.PKiKP(:,2),deg);
    case {'PKIKP' 'PKPdf'}
        slow=interp1(dp.PKPdf(:,1),dp.PKPdf(:,2),deg);
    case {'S' 'Slowermantle'}
        slow=interp1(dp.Slowermantle(:,1),dp.Slowermantle(:,2),deg);
    case {'SS' 'SSlowermantle'}
        slow=interp1(dp.SSlowermantle(:,1),dp.SSlowermantle(:,2),deg);
    case 'SKS'
        % iasp91
        slow=interp1(dp.SKS(:,1),dp.SKS(:,2),deg);
    case 'SKiKS'
        % iasp91
        slow=interp1(dp.SKiKS(:,1),dp.SKiKS(:,2),deg);
    case 'SKIKS'
        % iasp91
        slow=interp1(dp.SKIKS(:,1),dp.SKIKS(:,2),deg);
    otherwise
        error('seizmo:deg2slowness:badPhase',...
            'Unknown phase: %s',ph);
end

end
