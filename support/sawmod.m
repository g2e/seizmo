function [c]=sawmod(a,b)
%SAWMOD    Returns a sawtooth modulus
%
%    Usage: [c]=sawmod(a,b)
%
%    Description: SAWMOD(A,B) is S.*(A-N.*X) where N=round(A./B) if B~=0
%     and S=1-2.*MOD(N,2).  Thus SAWMOD(A,B) is always within the range
%     +/-B and forms a continuous function (but discontinuous in the 1st
%     dirivative).  The function is mainly useful for unwrapping latitude
%     values back to a valid range.  The inputs A and B must be real arrays
%     of the same size, or real scalars.
%
%    Notes:
%     By convention:
%      sawmod(A,0) returns A
%      sawmod(A,A) returns A
%
%    Tested on: Matlab r2007b
%
%    Examples:
%     To get latitude values LAT within the range of +/-90:
%      LAT=sawmod(LAT,90);
%
%    See also: cmod, mod, rem

%     Version History:
%        Mar. 24, 2009 - initial version
%        Apr. 23, 2009 - move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 22:20 GMT

n=round(0.5.*a./b);
s=1-2*mod(n,2);
c=s.*(a-2.*n.*b);
if(isscalar(b))
    if(b==0)
        c=a;
    end
else
    d=b==0;
    c(d)=a(d);
end

end
