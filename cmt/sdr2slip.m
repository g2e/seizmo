function [n,e,u]=sdr2slip(strike,dip,rake)
%SDR2SLIP    Calculates the slip vector for a given strike, dip and rake
%
%    Usage:    neu=sdr2slip(sdr)
%              neu=sdr2slip(strike,dip,rake)
%              [n,e,u]=sdr2slip(...)
%
%    Description:
%     NEU=SDR2SLIP(SDR) calculates the slip vector NEU for a given fault
%     and rake defined by SDR.  SDR is expected to be a Nx3 array formatted
%     as [strike dip rake] where N allows for multiple planes & rakes to be
%     processed simultaneously.  Strike is positive clockwise from North,
%     dip is positive downward from the horizontal and rake is positive
%     counter-clockwise in the fault plane from the strike direction.  Note
%     that the strike must be such that when you look along the direction
%     of the strike the fault dips to your right.  NEU is returned as a Nx3
%     array formatted as [North East Up].  The slip vector defines the
%     motion of the hanging wall relative to the foot wall.
%
%     NEU=SDR2SLIP(STRIKE,DIP,RAKE) allows strike, dip & rake to be given
%     separately.
%
%     [N,E,U]=SDR2SLIP(...) returns North, East & Up separately.
%
%    Notes:
%
%    Examples:
%     % Slip vector of a thrust fault:
%     sdr2slip(0,30,90)
%
%    See also: SDR2NULL, STRIKEDIP2NORM, NORM2STRIKEDIP, AUXPLANE, SDR2TPB,
%              TPB2SDR, NORMSLIP2SDR, NODALLINES

%     Version History:
%        Mar. 15, 2013 - initial version
%        Mar. 18, 2013 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 18, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% one or both inputs
switch nargin
    case 1
        if(size(strike,2)~=3 || ndims(strike)>2)
            error('seizmo:sdr2slip:badInput',...
                'SDR must be a Nx3 array as [STRIKE DIP RAKE] !');
        elseif(~isnumeric(strike) || ~isreal(strike))
            error('seizmo:sdr2slip:badInput',...
                'SDR must be a real-valued Nx3 array!');
        end
        [strike,dip,rake]=deal(strike(:,1),strike(:,2),strike(:,3));
    case 3
        if(~isnumeric(strike) || ~isreal(strike) ...
                || ~isnumeric(dip) || ~isreal(dip) ...
                || ~isnumeric(rake) || ~isreal(rake))
            error('seizmo:sdr2slip:badInput',...
                'STRIKE/DIP/RAKE must be real-valued arrays!');
        end
        [strike,dip,rake]=expandscalars(strike,dip,rake);
    otherwise
        error('seizmo:sdr2slip:badNumInputs',...
            'Incorrect number of inputs (only 1 or 3)!');
end

% get slip vector
n=cosd(rake).*cosd(strike)+sind(rake).*cosd(dip).*sind(strike);
e=cosd(rake).*sind(strike)-sind(rake).*cosd(dip).*cosd(strike);
u=sind(rake).*sind(dip);

% combine if only one output
if(nargout<=1); n=[n(:) e(:) u(:)]; end

end
