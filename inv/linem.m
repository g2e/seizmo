function [m,covm] = linem(dd,tt,covd)
%LINEM    Straight line fit to traveltime data with simple least squares
%
%    Usage:  [m,covm]=linem(dd,tt)
%            [m,covm]=linem(dd,tt,covd)
%
%    Description: [M,COVM]=LINEM(DD,TT) uses the degree distances DD and
%     the travel times TT to solve for the best straight-line fit to the
%     travel time curve using a least squares inversion.  The first output
%     is the y-intercept and slope organized in M as : [y-intercept slope].
%     The second output is the corresponding model covariance matrix COVM
%     assuming that the data are uncorrelated and all have equal variances
%     of 1.
%
%     [M,COVM]=LINEM(DD,TT,COVD) includes the supplied data covariance
%     matrix COVD in the computation of the model covariance matrix COVM.
%
%    Notes:
%     - Make sure that your data covariance matrix is filled with variances
%       and covariances, NOT standard deviations!
%
%    Examples:
%     Fit a line through some predicted P arrivals:
%      data=addarrivals(data);
%      P=getarrival(data,'P');
%      gcarc=getheader(data,'gcarc');
%      [m,covm]=linem(gcarc,P);
%
%    See also: correlate, wlinem

%     Version History:
%        Aug. 25, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 25, 2009 at 16:15 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end

% default covd
ndd=numel(dd);
if(nargin==2); covd=diag(ones(ndd,1)); end

% make sure dd,tt are equal sized column vectors
if(ndd~=numel(tt))
    error('seizmo:linem:badInput',...
        ['Distance and Traveltime arrays must\n'...
        'have the same number of elements!']);
end
dd=dd(:);
tt=tt(:);

% make sure covd is square and correct size
sz=size(covd);
if(numel(sz)>2 || sz(1)~=sz(2) || sz(1)~=ndd)
    error('seizmo:linem:badInput',...
        ['Data covariance matrix must be square and\n'...
        'match with Distance and Traveltime arrays!']);
end

% kernel matrix
G=[ones(ndd,1) dd];

% generalized inverse of kernel matrix, simple least squares
M=inv(G.'*G)*G.';

% covariance of model parameters
covm=M*covd*M.';

% least squares solution
m=M*tt;

end
