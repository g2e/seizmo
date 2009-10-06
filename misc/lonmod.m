function [c]=lonmod(a,b)
%LONMOD    Returns a longitude modulus
%
%    Usage:    c=lonmod(a,b)
%
%    Description: LONMOD(A,B) is A-N.*X where N=round(A./B) if B~=0.  Thus
%     LONMOD(A,B) is always within the range +/-(B/2).  The inputs A and B
%     must be real arrays of the same size, or real scalars.
%
%    Notes:
%     By convention:
%      lonmod(A,0) returns A
%      lonmod(A,A) returns 0
%
%    Examples:
%     To get longitude values LON within the range of +/-180:
%      LON=lonmod(LON,360);
%
%    See also: latmod, mod, rem

%     Version History:
%        Mar. 24, 2009 - initial version
%        Apr. 23, 2009 - move usage up
%        Sep.  7, 2009 - minor doc update
%        Sep. 30, 2009 - changed name from CMOD to LONMOD
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 30, 2009 at 15:30 GMT

% todo:

c=a-round(a./b).*b;
if(isscalar(b))
    if(b==0)
        c=a;
    end
else
    d=b==0;
    c(d)=a(d);
end

end
