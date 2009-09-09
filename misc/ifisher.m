function [R]=ifisher(Z)
%IFISHER    Converts Z statistics to correlation coefficients
%
%    Usage:    r=ifisher(z)
%
%    Description:  R=IFISHER(Z) uses the inverse Fisher transform to
%     convert z statistics Z to correlation coefficients R.  Correlation
%     coefficients will be bounded by -1 and 1.  NaNs will be preserved as
%     NaNs.
%
%    Notes:
%
%    Examples:
%     Convert correlations to z statistic get the mean and std dev and
%      convert them back to correlation values:
%       z=fisher(r)
%       z_mean=mean(z)
%       z_std=sqrt(var(z))
%       [r_lower,r_mean,r_upper]=ifisher(z_mean-z_std,z_mean,z_mean+z_std)
%
%    See also: fisher

%     Version History:
%        Sep.  9, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  9, 2009 at 01:50 GMT

% todo:

% INVERSE FISHER TRANSFORM -> Z STATISTIC SPACE TO CORRELATION SPACE
R=(exp(2*Z)-1)./(exp(2*Z)+1);
R(Z==inf)=1; % otherwise inf returns nan

end
