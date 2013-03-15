function [n,e,u]=sdr2null(strike,dip,rake)
%SDR2NULL    Calculates the null vector for a given strike, dip and rake
%
%    Usage:    neu=sdr2null(sdr)
%              neu=sdr2null(strike,dip,rake)
%              [n,e,u]=sdr2null(...)
%
%    Description:
%     NEU=SDR2NULL(SDR) calculates the null vector NEU for a given fault
%     and rake defined by SDR.  SDR is expected to be a Nx3 array formatted
%     as [strike dip rake] where N allows for multiple planes & rakes to be
%     processed simultaneously.  Strike is clockwise from North, dip is
%     positive downward from the horizontal and rake is counter-clockwise
%     in the fault plane from the strike direction.  Note that the strike
%     must be such that when you look along the direction of the strike the
%     fault dips to your right.  NEU is returned as a Nx3 array formatted
%     as [North East Up].
%
%     NEU=SDR2NULL(STRIKE,DIP,RAKE) allows strike, dip & rake to be given
%     separately.
%
%     [N,E,U]=SDR2NULL(...) returns North, East & Up separately.
%
%    Notes:
%
%    Examples:
%     % Null vector of a thrust fault:
%     sdr2null(0,30,90)
%
%    See also: SDR2SLIP, STRIKEDIP2NORM, NORM2STRIKEDIP, AUXPLANE, SDR2TPB,
%              TPB2SDR, NORMSLIP2SDR, NODALLINES

%     Version History:
%        Mar. 15, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% one or both inputs
switch nargin
    case 1
        if(size(strike,2)~=3 || ndims(strike)>2)
            error('seizmo:sdr2null:badInput',...
                'SDR must be a Nx3 array as [STRIKE DIP RAKE] !');
        elseif(~isnumeric(strike) || ~isreal(strike))
            error('seizmo:sdr2null:badInput',...
                'SDR must be a real-valued Nx3 array!');
        end
        [strike,dip,rake]=deal(strike(:,1),strike(:,2),strike(:,3));
    case 3
        if(~isnumeric(strike) || ~isreal(strike) ...
                || ~isnumeric(dip) || ~isreal(dip) ...
                || ~isnumeric(rake) || ~isreal(rake))
            error('seizmo:sdr2null:badInput',...
                'STRIKE/DIP/RAKE must be real-valued arrays!');
        end
        [strike,dip,rake]=expandscalars(strike,dip,rake);
    otherwise
        error('seizmo:sdr2null:badNumInputs',...
            'Incorrect number of inputs (only 1 or 3)!');
end

% get null vector
n=-sind(rake).*cosd(strike)+cosd(rake).*cosd(dip).*sind(strike);
e=-sind(rake).*sind(strike)-cosd(rake).*cosd(dip).*cosd(strike);
u=cosd(rake).*sind(dip);

% combine if only one output
if(nargout<=1); n=[n(:) e(:) u(:)]; end

end
