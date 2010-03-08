function [normal,slip]=sdr2ns(strike,dip,rake)
%SDR2ND    Returns normal & slip vectors (NEU) for a given strike-dip-rake
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also:

%     Version History:
%        Mar.  8, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2010 at 13:50 GMT

% todo:

% for radians/degrees
R2D=180/pi;

if(nargin==1)
    
elseif(nargin==3)
    % check inputs
    
    % expand scalars
    
    % convert to radians
    strike=strike/R2D;
    dip=dip/R2D;
    rake=rake/R2D;
    
    % normal vector
    normal=[-sin(dip).*sin(strike) sin(dip).*cos(strike) cos(dip)];
    
    % slip vector
    slip=[cos(rake).*cos(strike)+sin(rake).*cos(dip).*sin(strike) ...
          cos(rake).*sin(strike)-sin(rake).*cos(dip).*cos(strike) ...
          sin(rake).*sin(dip)];
else
    error('seizmo:sdr2ns:badNumInputs',...
        'Incorrect number of inputs!');
end

end
