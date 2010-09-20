function [pol]=polysqrt(p,tol)
%POLYSQRT    Returns the square root of a polynomial if it exists
%
%    Usage:    pol=polysqrt(p,tol)
%
%    Description:
%     POL=POLYSQRT(P) returns a vector P, if it exists, such that
%     conv(POL,POL)=P. P is a vector whose elements are the coefficients 
%     of a polynomial in descending powers.
%
%    Notes:
%     - FEX #75171
%
%    Examples:
%     % Some simple examples:
%     polysqrt([1 2 1])
%     polysqrt([4 -12 9])
%     polysqrt(conv(1:5,1:5))
%
%    See also: POLYFIT, POLYVAL, POLY, ROOTS, POLYSTR

%     Version History:
%        Nov.  4, 2009 - initial version
%        Sep. 19, 2010 - doc update and some code refactoring
%
%     Written by Andre Fioravanti (fioravanti at gmail dot com)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 19, 2010 at 10:00 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default tol
if(nargin<2 || isempty(tol)); tol=1e-6*norm(p,1); end

% check inputs
if(~isreal(p) || ~isvector(p))
    error('seizmo:polysqrt:badInput',...
        'P must be a vector of real-valued polynomial coefficients!');
elseif(~isreal(tol) || ~isscalar(tol) || tol<0)
    error('seizmo:polysqrt:badInput',...
        'TOL must be a positive real!');
end

% get polynomial degree
deg_p=numel(p)-1;

% no square root if not even
if(mod(deg_p,2))
    pol=[];
    return;
end

% order of output poly
n_var=deg_p/2+1;
pol=zeros(1,n_var);

% first coefficient is just the root of the input poly's first coefficient
pol(1)=sqrt(p(1));

% determine the remaining coefficients
for pos=2:n_var
    indep_term=p(pos)-calc_delta(pol(2:pos-1));
    pol(pos)=indep_term/(2*pol(1));
end

% check if polynomial is the square root of the other
% and return nothing if it is not within the tolerance
test_pol=conv(pol,pol);
res=max(abs(p-test_pol));
if(res>tol); pol=[]; end

end
    
function res=calc_delta(t) 
k=numel(t);
res=0;
for i=1:k
    res=res+t(i)*t(k-i+1);
end
end
