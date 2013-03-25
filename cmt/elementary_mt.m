function [varargout]=elementary_mt(n)
%ELEMENTARY_MT    Returns one of six elementary moment tensors
%
%    Usage:    mt=elementary_mt(idx)
%              [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=elementary_mt(idx)
%
%    Description:
%     MT=ELEMENTARY_MT(IDX) returns elementary moment tensors commonly
%     used in moment tensor inversion.  There are 6 different elementary
%     moment tensors: 1 & 2 correspond to strike-slip on vertical faults,
%     3 & 4 are dip-slip on vertical faults, 5 is a 45deg dipping dip-slip
%     fault, & 6 is an isotropic source (explosion).  IDX should be a
%     scalar or vector of integers from 1 to 6 indicating the elementary
%     moment tensor(s) to return.  IDX is optional and the default is 1:6
%     which returns all of the tensors.  MT is a Nx6 array of moment
%     tensor components where N is the number of integers in IDX.  MT
%     components are in the Harvard format (Up, South, East axes) giving
%     a layout of [Mrr Mtt Mpp Mrt Mrp Mtp].
%
%     [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]=ELEMENTARY_MT(IDX) outputs the moment
%     tensor components separately.
%
%    Notes:
%     - References:
%        Kikuchi and Kanamori (1991)
%        Aki and Richards (2002)
%
%    Examples:
%     % Plot the 6 elementary moment tensors:
%     plotmt(1:6,zeros(1,6),elementary_mt(1:6));
%     axis tight equal off
%
%    See also: PLOTMT

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 22, 2010 - added docs
%        Feb. 11, 2011 - mass nargchk fix
%        June  1, 2011 - improve docs, use harvard format
%        Mar. 19, 2013 - doc update, cleanup output code
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 19, 2013 at 15:05 GMT

% todo:

% the six elementary moment tensors
m=[ 0  0  0  0  0 -1
    0  1 -1  0  0  0
    0  0  0  0 -1  0
    0  0  0  1  0  0
    1 -1  0  0  0  0
    1  1  1  0  0  0];

% check nargin
error(nargchk(1,1,nargin));

% check n
if(~isreal(n) || any(~ismember(n,1:6)))
    error('seizmo:elementary_mt:badInput',...
        'N must be 1 or more integers with value from 1 to 6!');
end
n=n(:);

% output as combined array or individual components
if(nargout<2)
    varargout={m(n,:)};
else
    varargout={m(n,1) m(n,2) m(n,3) m(n,4) m(n,5) m(n,6)};
end

end
