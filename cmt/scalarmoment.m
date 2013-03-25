function [mo]=scalarmoment(varargin)
%SCALARMOMENT    Returns the scalar moment of a moment tensor
%
%    Usage:    mo=scalarmoment(mt)
%              mo=scalarmoment(m1,m2,m3,m4,m5,m6)
%
%    Description:
%     MO=SCALARMOMENT(MT) returns the scalar moment of the moment tensor
%     MT in dyne*cm (a measure of energy).  MT must be in a format
%     recognized by MT_CHECK.  MO is a Nx1 column vector.
%
%     MO=SCALARMOMENT(M1,M2,M3,M4,M5,M6) allows inputting the moment tensor
%     components individually (Harvard or Aki & Richards system).
%
%    Notes:
%     - Scalar moment or seismic moment is used to define the size/energy
%       of an earthquake.  The moment is related to the earthquake by the
%       equation Mo = G*A*D where G is the shear modulus of the involved
%       rock (G=Vs^2*rho), A is the area of the fault that slipped and D
%       is the magnitude of the slip on the fault.
%
%    Examples:
%     % Compare scalar moments in the CMT catalog with those calculated
%     % directly from the tensor values (note the bias due to truncation):
%     cmts=findcmts;
%     mo=scalarmoment(cmts);
%     hist(cmts.scalarmoment-mo,100);
%
%     % How does the principal axes compare?
%     t=[cmts.eigval1 cmts.plunge1 cmts.azimuth1];
%     p=[cmts.eigval3 cmts.plunge3 cmts.azimuth3];
%     b=[cmts.eigval2 cmts.plunge2 cmts.azimuth2];
%     mo2=scalarmoment(tpb2mt(t,p,b));
%     hist(cmts.scalarmoment-mo2,100);
%
%     % How to include the exponent:
%     mo=scalarmoment(cmts).*10.^cmts.exponent;
%
%    See also: FINDCMTS, FINDCMT, MOMENTMAG, MO2HD

%     Version History:
%        Mar. 11, 2011 - initial version
%        June  7, 2011 - made warning more informative
%        Mar. 19, 2013 - doc update
%        Mar. 25, 2013 - use mt_change/mt_check
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 23:55 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check input & convert to grid style
error(mt_check(varargin{:}));
mt=mt_change('g',varargin{:});

% calculate moment
mo=nan(sz(3),1);
for i=1:sz(3)
    mo(i,1)=sqrt(trace(mt(:,:,i)*mt(:,:,i))/2);
end

end
