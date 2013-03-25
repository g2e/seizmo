function [mt]=sdr2mt(strike,dip,rake)
%SDR2MT    Convert strike-dip-rake to moment tensor
%
%    Usage:    mt=sdr2mt(sdr)
%              mt=sdr2mt(strike,dip,rake)
%
%    Description:
%     MT=SDR2MT(SDR) converts a double-couple (aka focal mechanism) given
%     in strike-dip-rake form (Nx3) as [strike dip rake] to Harvard moment
%     tensor form (Nx6) as [Mrr Mtt Mpp Mrt Mrp Mtp] which uses the Up,
%     South, East coordinate system.  N allows for multiple focal
%     mechanisms to be converted simultaneously.  The main purpose of this
%     function is to allow plotting of the focal mechanism with PLOTMT.
%     Strike is positive clockwise from North, dip is positive downward
%     from the horizontal and rake is positive counter-clockwise in the
%     fault plane from the strike direction.  Note that the strike must be
%     such that when you look along the direction of the strike the fault
%     dips to your right.
%
%     MT=SDR2MT(STRIKE,DIP,RAKE) allows inputing the strike, dip & rake
%     separately.  Note that the inputs should all be Nx1 column vectors or
%     scalars.
%
%    Notes:
%     - Tested OK against George Helffrich's website:
%        http://www1.gly.bris.ac.uk/~george/focmec.html
%     - See Aki & Richards (2002) Figure 4.20 & Box 4.4 for details.
%
%    Examples:
%     % Plot a dip-slip fault with EW strike and 45 deg dip:
%     plotmt(0,0,sdr2mt(90,45,-90))
%     axis tight equal off
%
%     % Validate SDR2MT & MT2SDR using GlobalCMT catalog:
%     cmts=findcmts;
%     sdr1=[cmts.strike1 cmts.dip1 cmts.rake1];
%     sdr2=[cmts.strike2 cmts.dip2 cmts.rake2];
%     sdr3=mt2sdr(sdr2mt(sdr1));
%     max(min(abs(azdiff(sdr3,sdr1)),abs(azdiff(sdr3,sdr2))))
%     sdr4=mt2sdr(sdr2mt(sdr2));
%     max(min(abs(azdiff(sdr4,sdr1)),abs(azdiff(sdr4,sdr2))))
%
%    See also: MT2SDR, AUXPLANE, PLOTMT, MT_DECOMP, MT_DIAG, MT_UNDIAG

%     Version History:
%        Mar.  8, 2010 - initial version
%        June  1, 2011 - harvard output, now with docs
%        Mar. 21, 2013 - doc update, clean up code, no mt cmp output
%        Mar. 25, 2013 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% one or both inputs
switch nargin
    case 1
        if(size(strike,2)~=3 || ndims(strike)>2)
            error('seizmo:sdr2mt:badInput',...
                'SDR must be a Nx3 array as [STRIKE DIP RAKE] !');
        elseif(~isnumeric(strike) || ~isreal(strike))
            error('seizmo:sdr2mt:badInput',...
                'SDR must be a real-valued Nx3 array!');
        end
        [strike,dip,rake]=deal(strike(:,1),strike(:,2),strike(:,3));
    case 3
        if(~isnumeric(strike) || ~isreal(strike) ...
                || ~isnumeric(dip) || ~isreal(dip) ...
                || ~isnumeric(rake) || ~isreal(rake))
            error('seizmo:sdr2mt:badInput',...
                'STRIKE/DIP/RAKE must be real-valued arrays!');
        end
        [strike,dip,rake]=expandscalars(strike(:),dip(:),rake(:));
    otherwise
        error('seizmo:sdr2mt:badNumInputs',...
            'Incorrect number of inputs (only 1 or 3)!');
end

% convert strike-dip-rake to HRV mt
% f = strike, d = dip, l = rake
%1 Mrr =  Mzz =  Mo sin2d sinl
%2 Mtt =  Mxx = -Mo(sind cosl sin2f +     sin2d sinl (sinf)^2 )
%3 Mpp =  Myy =  Mo(sind cosl sin2f -     sin2d sinl (cosf)^2 )
%4 Mrt =  Mxz = -Mo(cosd cosl cosf  +     cos2d sinl sinf )
%5 Mrp = -Myz =  Mo(cosd cosl sinf  -     cos2d sinl cosf )
%6 Mtp = -Mxy = -Mo(sind cosl cos2f + 0.5 sin2d sinl sin2f )
mt=nan(numel(strike),6);
mt(:,1)=sind(2*dip).*sind(rake);
mt(:,2)=sind(dip).*cosd(rake).*sind(2*strike) ...
    + sind(2*dip).*sind(rake).*sind(strike).^2;
mt(:,3)=sind(dip).*cosd(rake).*sind(2*strike) ...
    - sind(2*dip).*sind(rake).*cosd(strike).^2;
mt(:,4)=cosd(dip).*cosd(rake).*cosd(strike) ...
    + cosd(2*dip).*sind(rake).*sind(strike);
mt(:,5)=cosd(dip).*cosd(rake).*sind(strike) ...
    - cosd(2*dip).*sind(rake).*cosd(strike);
mt(:,6)=sind(dip).*cosd(rake).*cosd(2*strike) ...
    + 0.5.*sind(2*dip).*sind(rake).*sind(2*strike);
mt(:,2:2:6)=-mt(:,2:2:6);

end
