function [mo]=scalarmoment(mt)
%SCALARMOMENT    Returns the scalar moment of a moment tensor
%
%    Usage:    mo=scalarmoment(mt)
%
%    Description:
%     MO=SCALARMOMENT(MT) returns the scalar moment of the moment tensor
%     MT in dyne*cm (a measure of energy).  MT must be a Nx6 or 3x3xN array
%     where N is the number of moment tensors.  MO is a Nx1 column vector.
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
%     mo=scalarmoment(mt_s2v(cmts));
%     hist(cmts.scalarmoment-mo,100);
%
%     % How does the principal axes compare?
%     t=[cmts.eigval1 cmts.plunge1 cmts.azimuth1];
%     p=[cmts.eigval3 cmts.plunge3 cmts.azimuth3];
%     b=[cmts.eigval2 cmts.plunge2 cmts.azimuth2];
%     mo2=scalarmoment(tpb2mt(t,p,b));
%     hist(cmts.scalarmoment-mo2,100);
%
%    See also: FINDCMTS, FINDCMT, MOMENTMAG, MO2HD

%     Version History:
%        Mar. 11, 2011 - initial version
%        June  7, 2011 - made warning more informative
%        Mar. 19, 2013 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 19, 2013 at 23:55 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check inputs
sz=size(mt);
if(~isnumeric(mt) || ~isreal(mt))
    error('seizmo:scalarmoment:badInput',...
        'MT must be a real-valued numeric array!');
elseif(numel(sz)==2 && sz(2)==6) % (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)
    % convert to 3x3
    mt=reshape(mt(:,[1 4 5 4 2 6 5 6 3])',[3 3 sz(1)]);
    sz=[3 3 sz(1)];
elseif(isequal(sz(1:2),[3 3])) % [Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp]
    % leave alone
    sz=[sz(1:2) prod(sz(3:end))];
else
    error('seizmo:scalarmoment:badInput',...
        'MT is not in a valid format!');
end
mo=nan(sz(3),1);
for i=1:sz(3)
    mo(i,1)=sqrt(trace(mt(:,:,i)*mt(:,:,i))/2);
end

end
