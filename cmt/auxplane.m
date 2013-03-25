function [strike,dip,rake]=auxplane(strike,dip,rake)
%AUXPLANE    Returns strike-dip-slip of 2nd (auxiliary) focal plane
%
%    Usage:  sdr=auxplane(sdr)
%            sdr=auxplane(strike,dip,rake)
%            [strike,dip,rake]=auxplane(...)
%
%    Description:
%     SDR=AUXPLANE(SDR) finds the auxiliary plane(s) and rake(s) that match
%     the strain(s) produced by the input plane(s) and rake(s).  SDR must
%     be given as an Nx3 array arranged as [strike dip rake] where N allows
%     for multiple planes & rakes to be processed simultaneously.  Strike
%     is clockwise from North, dip is positive downward from the horizontal
%     and rake is counter-clockwise in the fault plane from the strike
%     direction.  Note that the strike must be such that when you look
%     along the direction of the strike the fault dips to your right.
%
%     SDR=AUXPLANE(STRIKE,DIP,RAKE) allows strike, dip & rake to be given
%     separately.
%
%     [STRIKE,DIP,RAKE]=AUXPLANE(...) returns strike, dip & rake
%     separately.
%
%    Notes:
%     - The auxiliary plane is perpendicular to the input fault plane.
%     - The normal & slip vector of the auxiliary plane is the slip &
%       normal vector of the fault plane (there may be a sign change
%       needed so that the normal vector points upwards).
%     - The plane perpendicular to both the fault plane and the auxiliary
%       plane has its normal along the null axis of the focal mechanism.
%
%    Examples:
%     % Corresponding auxiliary fault plane and rake for
%     % dip-slip on a normal fault that strikes North:
%     auxplane(0,45,-90)
%
%    See also: NORM2STRIKEDIP, STRIKEDIP2NORM, SDR2SLIP, SDR2NULL, SDR2TPB,
%              TPB2SDR, NORMSLIP2SDR, NODALLINES

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 22, 2010 - added docs
%        June  1, 2011 - updated docs
%        Mar. 14, 2013 - rewrite
%        Mar. 15, 2013 - more cleanup
%        Mar. 18, 2013 - sign change required for upwards normal
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
            error('seizmo:auxplane:badInput',...
                'SDR must be a Nx3 array as [STRIKE DIP RAKE] !');
        elseif(~isnumeric(strike) || ~isreal(strike))
            error('seizmo:auxplane:badInput',...
                'SDR must be a real-valued Nx3 array!');
        end
        [strike,dip,rake]=deal(strike(:,1),strike(:,2),strike(:,3));
    case 3
        if(~isnumeric(strike) || ~isreal(strike) ...
                || ~isnumeric(dip) || ~isreal(dip) ...
                || ~isnumeric(rake) || ~isreal(rake))
            error('seizmo:auxplane:badInput',...
                'STRIKE/DIP/RAKE must be real-valued arrays!');
        end
        [strike,dip,rake]=expandscalars(strike,dip,rake);
    otherwise
        error('seizmo:auxplane:badNumInputs',...
            'Incorrect number of inputs (only 1 or 3)!');
end

% size to reshape to
sz=size(dip);

% get strike/dip/rake of auxiliary plane
% -> sdr to normal, slip
% -> switch normal & slip
%    -> taking care to keep the normal pointed upwards so that it points to
%       the hanging wall (also have to flip the slip so it gives the motion
%       of the hanging wall relative to the foot wall)
%    -> note that this includes a trick such that +val=>1, 0=>1, -val=>-1
% -> convert back to sdr
normal=strikedip2norm(strike,dip);
slip=sdr2slip(strike,dip,rake);
[strike,dip,rake]=normslip2sdr((1-2*(slip(:,[3 3 3])<0)).*slip,...
    (1-2*(slip(:,[3 3 3])<0)).*normal);

% reshape
strike=reshape(strike,sz);
dip=reshape(dip,sz);
rake=reshape(rake,sz);

% combine if only one output
if(nargout<=1); strike=[strike(:) dip(:) rake(:)]; end

end
