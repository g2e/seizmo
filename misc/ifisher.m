function [R]=ifisher(Z)
%IFISHER    Converts Z statistics to correlation coefficients
%
%    Usage:    r=ifisher(z)
%
%    Description:
%     R=IFISHER(Z) uses the inverse Fisher transform to convert z
%     statistics Z to correlation coefficients R.  Correlation coefficients
%     will be bounded by -1 and 1.  NaNs will be preserved as NaNs.
%
%    Notes:
%
%    Examples:
%     % Convert correlations to z statistic get the mean and
%     % std dev and convert them back to correlation values:
%     z=fisher(r);
%     z_mean=mean(z);
%     z_std=std(z);
%     r_lower=ifisher(z_mean-z_std);
%     r_mean=ifisher(z_mean);
%     r_upper=ifisher(z_mean+z_std);
%
%    See also: FISHER

%     Version History:
%        Sep.  9, 2009 - minor doc update
%        Mar. 11, 2010 - fixed example
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 10:35 GMT

% todo:

% INVERSE FISHER TRANSFORM -> Z STATISTIC SPACE TO CORRELATION SPACE
R=(exp(2*Z)-1)./(exp(2*Z)+1);
R(Z==inf)=1; % otherwise inf returns nan

end
