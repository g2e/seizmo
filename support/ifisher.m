function [R]=ifisher(Z)
%IFISHER    Converts Z statistics to correlation coefficients
%
%    Description:  This uses the inverse Fisher transform to convert z 
%     statistics to correlation coefficients.  Correlation coefficients 
%     will be bounded by -1 and 1.  NaNs will be preserved as NaNs.
%
%    Usage: r=ifisher(z)
%
%    Examples:
%     Convert correlations to z statistic get the mean and std dev and
%      convert them back to correlation values
%       z=fisher(r)
%       z_mean=mean(z)
%       z_std=sqrt(var(z))
%       [r_lower,r_mean,r_upper]=ifisher(z_mean-z_std,z_mean,z_mean+z_std)
%
%    See also: fisher

% INVERSE FISHER TRANSFORM -> Z STATISTIC SPACE TO CORRELATION SPACE
R=(exp(2*Z)-1)./(exp(2*Z)+1);
R(Z==inf)=1; % otherwise inf returns nan

end
