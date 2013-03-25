function [strike,dip,rake]=mt2sdr(varargin)
%MT2SDR    Converts moment tensors to strike-dip-rake using decomposition
%
%    Usage:    sdr=mt2sdr(mt)
%              sdr=mt2sdr(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
%              [strike,dip,rake]=mt2sdr(...)
%
%    Description:
%     SDR=MT2SDR(MT) decomposes the moment tensors in MT to get the maximum
%     double-couple component and then converts that to strike, dip & rake.
%     MT must be a scalar struct as output by FINDCMT/FINDCMTS, a Nx6
%     array, or a 3x3xN array where N is the number of moment tensors in
%     MT.  Moment tensors in MT are expected to be in Harvard convention
%     (Up, South, East).  The output SDR is a Nx3 array formatted as
%     [strike dip rake].  Use the function AUXPLANE to get the strike, dip,
%     & rake of the auxiliary plane.  Strike is positive clockwise from
%     North, dip is positive downward from the horizontal, and rake is
%     positive counter-clockwise in the fault plane from the strike
%     direction.  Note that the strike is given such that when you look
%     along the direction of the strike the fault dips to your right.
%
%     SDR=MT2SDR(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) allows specifying the moment
%     tensor components individually (Harvard system only).
%
%     [STRIKE,DIP,RAKE]=MT2SDR(...) returns strike, dip, and rake
%     separately.
%
%    Notes:
%
%    Examples:
%     % Comparison to published GlobalCMT fault solutions
%     % shows that the output is accurate to within +/-0.5deg:
%     cmts=findcmts;
%     max([cmts.strike1 cmts.dip1 cmts.rake1]-mt2sdr(cmts))
%
%    See also: SDR2MT, AUXPLANE, PLOTMT, TPB2MT, MT2TPB, MT_DECOMP,
%              MT_DIAG, MT_UNDIAG, TPB2SDR, SDR2TPB

%     Version History:
%        Mar. 20, 2013 - initial version
%        Mar. 25, 2013 - update for mt_check/mt_change
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check tensor format
error(mt_check(varargin{:}));
mt=mt_change('g',varargin{:});

% 1) decompose into maximum double-couple + clvd
% 2) convert to principal axes in vpa
% 3) convert principal axes to strike/dip/rake
[dblcpl,clvd,vec]=mt_decomp(mt,'maxdc');
[t,p,b]=mt2tpb(mt_undiag(dblcpl,vec));
[strike,dip,rake]=tpb2sdr(t,p,b);

% combine outputs if desired
if(nargout<=1)
    strike=[strike(:) dip(:) rake(:)];
end

end
