function [mean]=nanmean(x,dim)
%NANMEAN    Return mean excluding NaN elements
%
%    Usage:    mean=nanmean(x)
%              mean=nanmean(x,dim)
%
%    Description: MEAN=NANMEAN(X) calculates the mean of X ignoring NaN
%     elements.  The mean is taken across the first non-singleton dimension
%     and the returned mean is equal in size to X except for the dimension
%     that the mean is found for, which has a size of 1.
%
%     MEAN=NANMEAN(X,DIM) specifies the dimension DIM that the mean is
%     taken across.  Use [] to get the default behavior.
%
%    Notes:
%
%    Examples:
%     Mean down the 3rd dimension of a random valued matrix with some NaNs:
%      x=rand(10,10,10);
%      x(x>0.4 & x<0.6)=nan;
%      mean=nanmean(x,3);
%
%    See also: MEAN, NANVARIANCE

%     Version History:
%        Oct. 14, 2009 - brought back, added checks, better documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2009 at 02:55 GMT

% todo:

% CHECK NARGIN
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% CHECKS
if(~isnumeric(x) && ~islogical(x))
    error('seizmo:nanmean:badInput','X must be numeric or logical!');
elseif(nargin>1 && ~isempty(dim) ...
        && (~isscalar(dim) || dim~=fix(dim) || dim<=0))
    error('seizmo:nanmean:badInput','DIM must be a valid dimension!');
end

% REMOVE INFLUENCE OF NaNs BY SETTING THEM TO 0
nans=isnan(x);
x(nans)=0;

% WHICH DIMENSION?
if(nargin>1 && ~isempty(dim))
    % USER DEFINED DIMENSION
    nne=sum(~nans,dim);    % GET DIVISOR
    nne(nne==0)=nan;       % IF NO NON-NaNs, RETURN NaN AS THE MEAN
    mean=sum(x,dim)./nne;  % MEAN EXCLUDING NaN ELEMENTS
else
    % FIRST NON-SINGLETON DIMENSION
    nne=sum(~nans);    % GET DIVISOR
    nne(nne==0)=nan;   % IF NO NON-NaNs, RETURN NaN AS THE MEAN
    mean=sum(x)./nne;  % MEAN EXCLUDING NaN ELEMENTS
end

end
