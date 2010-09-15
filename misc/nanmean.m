function [mean]=nanmean(x,dim)
%NANMEAN    Return mean excluding NaNs
%
%    Usage: mean=nanmean(x)
%           mean=nanmean(x,dim)
%
%    Description:
%     MEAN=NANMEAN(X) returns the means along the first non-
%     singleton dimension of X excluding NaN elements.  MEAN is equal in
%     size to X except for the dimension that the mean is computed for,
%     which has a size of 1.
%
%     MEAN=NANMEAN(X,DIM) specifies the dimension DIM that the mean is
%     taken across.  Use [] to get the default behavior.
%
%    Notes:
%     - Equivalent to NANMEAN of the Statistics Toolbox.
%
%    Examples:
%     % Mean down the 3rd dimension of a
%     % random valued matrix with some NaNs:
%     x=rand(10,10,10);
%     x(x>0.4 & x<0.6)=nan;
%     mean=nanmean(x,3)
%
%    See also: MEAN, NANVARIANCE

%     Version History:
%        Sep. 13, 2010 - brought back, better documentation, nargchk
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 13, 2010 at 14:45 GMT

% CHECK NARGIN
error(nargchk(1,2,nargin));

% IDENTIFY NaN ELEMENTS AND SET THEM TO ZERO (NO INFLUENCE ON SUM)
nans=isnan(x);
x(nans)=0;

% LET SUM DECIDE WHICH DIMENSION TO WORK ON
if(nargin==1 || isempty(dim))
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
