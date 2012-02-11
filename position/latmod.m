function [c,p]=latmod(a,b)
%LATMOD    Returns a latitude modulus (ie unwraps latitudes)
%
%    Usage:    lat=latmod(lat)
%              c=latmod(a,b)
%              [c,p]=latmod(a,b)
%
%    Description:
%     LAT=LATMOD(LAT) returns latitudes LAT to be within the range +/-90.
%     For example, if LAT=100 then the output is 80.  See the other usage
%     forms for more algorithm details.
%
%     LATMOD(A,B) is S.*(A-N.*B) where N=round(A./B) if B~=0 and
%     S=1-2.*MOD(N,2).  Thus LATMOD(A,B) is always within the range +/-B
%     and forms a continuous function (but discontinuous in the 1st
%     derivative) that looks like a "triangle wave".  The function is
%     primarily intended for unwrapping latitude values back to a valid
%     range.  The inputs A and B must be real arrays of the same size, or
%     real scalars.  B is optional and defaults to 90.
%
%     [C,P]=LATMOD(A,B) also returns the number of pole-crossings each
%     latitude in A made in P.  For example, if A=100 & B=90 then C=80 &
%     P=1.  This is useful for coupling LATMOD with LONMOD to preserve the
%     positions while reducing the values to reasonable ranges.  See the
%     Examples section below for instructions on how to use this output.
%
%    Notes:
%     By convention:
%      [C,P]=latmod(A,0) returns C=A, P=Inf
%      [C,P]=latmod(A,A) returns C=A, P=0
%
%    Examples:
%     % Modifying the latitude should also take into account the longitude
%     % shift necessary to preserve the actual position.  This may be done
%     % by utilizing the second output to shift the longitude by 180
%     % degrees if there are an odd number of pole-crossings:
%     [LAT,P]=latmod(LAT,90);
%     LON=lonmod(LON+mod(P,2)*180,360);
%
%     % or (for convenience)
%     [LAT,LON]=fixlatlon(LAT,LON);
%
%    See also: LONMOD, FIXLATLON, MOD, REM

%     Version History:
%        Mar. 24, 2009 - initial version
%        Apr. 23, 2009 - move usage up
%        Sep.  7, 2009 - minor doc update
%        Sep. 30, 2009 - changed name from SAWMOD to LATMOD
%        Feb. 16, 2010 - added P output, doc update
%        Feb.  9, 2012 - B=90 by default, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 01:00 GMT

% todo:

if(nargin<2 || isempty(b)); b=90; end
n=round(0.5.*a./b);
s=1-2*mod(n,2);
c=s.*(a-2.*n.*b);
p=ceil((abs(a)-b)./(2*b));
if(isscalar(b))
    if(b==0)
        c=a;
    end
else
    d=b==0;
    c(d)=a(d);
end

end
