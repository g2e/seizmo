function [mean]=nanmean(x,dim)
%NANMEAN    Return mean excluding NaNs
%
%    Description:  Returns the means along a dimension excluding NaN 
%     elements.  Operates exactly like mean.
%
%    Usage: mean=nanmean(x)
%           mean=nanmean(x,dim)
%
%    Examples:
%
%    See also: mean, nanstd

% IDENTIFY NaN ELEMENTS AND SET THEM TO ZERO (NO INFLUENCE ON SUM)
nans=isnan(x);
x(nans)=0;

% LET SUM DECIDE WHICH DIMENSION TO WORK ON
if(nargin==1)
    nne=sum(~nans);     % COUNT NUMBER OF NON-NaN ELEMENTS (FOR THE AVERAGING)
    nne(nne==0)=nan;    % IF EVERY ELEMENT IS NaN RETURN NaN AS THE MEAN
    mean=sum(x)./nne;   % MEAN EXCLUDING NaN ELEMENTS (UNLESS ALL NaN)
% USER DEFINED DIMENSION
else
    nne=sum(~nans,dim);     % COUNT NUMBER OF NON-NaN ELEMENTS (FOR THE AVERAGING)
    nne(nne==0)=nan;        % IF EVERY ELEMENT IS NaN RETURN NaN AS THE MEAN
    mean=sum(x,dim)./nne;   % MEAN EXCLUDING NaN ELEMENTS (UNLESS ALL NaN)
end

end
