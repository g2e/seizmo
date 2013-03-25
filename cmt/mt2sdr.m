function [strike,dip,rake]=mt2sdr(mt)
%MT2SDR    Converts moment tensors to strike-dip-rake using decomposition
%
%    Usage:    sdr=mt2sdr(mt)
%              [strike,dip,rake]=mt2sdr(mt)
%
%    Description:
%     SDR=MT2SDR(MT) decomposes the moment tensors in MT to get the maximum
%     double-couple component and then converts that to strike, dip & rake.
%     MT must be a Nx6 or 3x3xN array where N is the number of moment
%     tensors in MT.  Moment tensors in MT are expected to be in Harvard
%     convention (Up, South, East).  The output SDR is a Nx3 array
%     formatted as [strike dip rake].  Use the function AUXPLANE to get the
%     strike, dip, & rake of the auxiliary plane.  Strike is positive
%     clockwise from North, dip is positive downward from the horizontal
%     and rake is positive counter-clockwise in the fault plane from the
%     strike direction.  Note that the strike must be such that when you
%     look along the direction of the strike the fault dips to your right.
%
%     [STRIKE,DIP,RAKE]=MT2SDR(MT) allows strike, dip, and rake to be
%     returned separately.
%
%    Notes:
%
%    Examples:
%     % Comparison to published GlobalCMT fault solutions
%     % shows that the output is accurate to within +/-0.5deg:
%     cmts=findcmts;
%     max([cmts.strike1 cmts.dip1 cmts.rake1]-mt2sdr(mt_s2v(cmts)))
%
%    See also: SDR2MT, AUXPLANE, PLOTMT, TPB2MT, MT2TPB, MT_DECOMP,
%              MT_DIAG, MT_UNDIAG, TPB2SDR, SDR2TPB

%     Version History:
%        Mar. 20, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 20, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check tensor format
mtsz=size(mt);
if(~isnumeric(mt) || ~isreal(mt))
    error('seizmo:mt2sdr:badInput',...
        'MT must be a real-valued numeric array!');
elseif(~(isequal(mtsz(1:2),[3 3]) && any(numel(mtsz)==[2 3])) ...
        && ~(mtsz(2)==6 && numel(mtsz)==2))
    error('seizmo:mt2sdr:badInput',...
        'MT must be a harvard moment tensor array as 3x3xN or Nx6!');
end

% 1) decompose into maximum double-couple + clvd
% 2) convert to principal axes in vpa
% 3) convert principal axes to strike/dip/rake
[dblcpl,clvd,vec]=mt_decomp(mt,'maxdc');
[t,p,b]=mt2tpb(mt_undiag(dblcpl,vec));
[strike,dip,rake]=tpb2sdr(t,p,b);

% combine outputs is desired
if(nargout<=1)
    strike=[strike(:) dip(:) rake(:)];
end

end
