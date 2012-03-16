function [m,covm] = wlinem(dd,tt,power,covd,We)
%WLINEM    Linear fit to traveltime data with weighted least squares
%
%    Usage:  [m,covm]=wlinem(dd,tt)
%            [m,covm]=wlinem(dd,tt,power)
%            [m,covm]=wlinem(dd,tt,power,covd)
%            [m,covm]=wlinem(dd,tt,power,covd,We)
%
%    Description:
%     [M,COVM]=WLINEM(DD,TT) uses the degree distances DD and the travel
%     times TT to solve for the best straight-line fit to the travel time
%     curve using a least squares inversion.  The first output is the
%     y-intercept and slope organized in M as : [y-intercept slope].  The
%     second output is the corresponding model covariance matrix COVM
%     assuming that the data are uncorrelated and all have equal variances
%     of 1 (ie the identity matrix).  The weighting matrix is also set
%     equal to the identity matrix for this case.
%
%     [M,COVM]=WLINEM(DD,TT,POWER) fits the travel times using an nth
%     degree polynomial, where n is POWER.  The default value of POWER is
%     1 (straight line fit).  Setting POWER to 0, finds the mean.  Setting
%     POWER to 2 gives a parabolic fit.  M is arranged so that the
%     polynomial coefficients are in ascending order (opposite of POLYFIT).
%
%     [M,COVM]=WLINEM(DD,TT,POWER,COVD) includes the data covariance matrix
%     COVD in the computation of the model covariance matrix COVM.  The
%     weighting matrix is set equal to the identity matrix for this case.
%
%     [M,COVM]=WLINEM(DD,TT,POWER,COVD,WE) utilizes the weighting matrix WE
%     in the inversion.
%
%    Notes:
%     - Make sure that your data covariance matrix is filled with variances
%       and covariances, NOT standard deviations!
%     - The weighting matrix usually only has non-zero elements on the
%       diagonal.  Make sure you know what your doing if your weighting
%       matrix has off diagonal elements!
%
%    Examples:
%     % Fit a line through some predicted P arrivals:
%     P=findpicks(arrivals2picks(data),'P',1);
%     gcarc=getheader(data,'gcarc');
%     [m,covm]=wlinem(gcarc,P);
%
%    See also: POLYFIT, POLYVAL

%     Version History:
%        Aug. 25, 2009 - initial version
%        Sep.  9, 2009 - minor doc fix, fix weighting matrix checks
%        Oct. 13, 2009 - minor doc update, add power option
%        Mar.  2, 2010 - minor doc update
%        Sep. 13, 2010 - nargchk fix
%        Feb. 11, 2011 - drop inv call
%        Mar.  5, 2012 - minor doc update
%        Mar. 15, 2012 - fix example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 15, 2012 at 18:00 GMT

% todo:

% check nargin
error(nargchk(2,5,nargin));

% default power, covd, We
ndd=numel(dd);
if(nargin<3 || isempty(power)); power=1; end
if(nargin<4 || isempty(covd))
    covd=sparse(1:ndd,1:ndd,ones(ndd,1),ndd,ndd);
end
if(nargin<5 || isempty(We))
    We=sparse(1:ndd,1:ndd,ones(ndd,1),ndd,ndd);
end

% make sure dd,tt are equal sized column vectors
if(ndd~=numel(tt))
    error('seizmo:wlinem:badInput',...
        ['Distance and Traveltime arrays must\n'...
        'have the same number of elements!']);
end
dd=dd(:);
tt=tt(:);

% make sure power is scalar, whole and >=0
if(~isscalar(power) || power~=fix(power) || power<0)
    error('seizmo:wlinem:badInput',...
        'POWER must be a scalar integer >=0!');
end

% make sure covd is square and correct size
sz=size(covd);
if(numel(sz)>2 || sz(1)~=sz(2) || sz(1)~=ndd)
    error('seizmo:wlinem:badInput',...
        ['Data covariance matrix must be square and\n'...
        'match with Distance and Traveltime arrays!']);
end

% make sure We is square and correct size
sz=size(We);
if(numel(sz)>2 || sz(1)~=sz(2) || sz(1)~=ndd)
    error('seizmo:wlinem:badInput',...
        ['Weighting matrix must be square and\n'...
        'match with Distance and Traveltime arrays!']);
end

% kernel matrix
p=1:power;
G=[ones(ndd,1) dd(:,ones(power,1)).^p(ones(ndd,1),:)];

% generalized inverse of kernel matrix, weighted least squares
Gg=(G.'*We*G)\G.'*We;

% covariance of model parameters
covm=Gg*covd*Gg.';

% least squares solution
m=Gg*tt;

end
