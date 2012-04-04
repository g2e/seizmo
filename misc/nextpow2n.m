function [p]=nextpow2n(n)
%NEXTPOW2N   Returns the next higher power of 2 for all array elements
%
%    Usage:    p=nextpow2n(n)
%
%    Description:
%     P=NEXTPOW2N(N) returns the next higher power of 2 for all elements in
%     N such that 2.^P>=abs(N).  This differs from NEXTPOW2 in that vectors
%     and arrays are operated on element-wise. 
%
%    Notes:
%     - No checks.
%
%    Examples:
%     % Get the next power of 2 for a sequence:
%     nextpow2n(1:10)
%
%    See also: NEXTPOW2, LOG2, POW2

%     Version History:
%        June 25, 2009 - initial version
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 08:25 GMT

% todo:

% n=f.*2.^p
[f,p]=log2(abs(n));

% handle exact powers of 2
if(~isempty(f)); k=(f==0.5); if(any(k)); p(k)=p(k)-1; end; end

% handle infinities and NaNs
k=~isfinite(f); p(k)=f(k);

end
