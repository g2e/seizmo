function []=ttsolve(xc,mri,pd,snr,info)
%TTSOLVE    Solves relative arrival times & polarities

% additional parameters
% - max number of iterations
% - when to force polarity consistency
%   - can we delay this?
%     - 0 == no delay
%     - 1 == force consistency after first refinement
%     - >mri == never force consistency
% - verbosity
verbose=seizmoverbose;

% check inputs

% convert correlation values to z-statistic
% - more appropriate for weighting
% - more appropriate for error estimation
xc.zg=fisher(xc.cg);
xc.zg(abs(xc.cg)>1-eps)=nan; % ignore perfect correlations

% initial alignment & polarity estimate
% - 2 choices:
%   - based on best correlating peaks (weighted by snr & z)
%   - times & polarities given as input
if(isempty(info))
    if(verbose); disp('Inverting for Initial Alignment & Polarity'); end
    [dt]=ttsolve(lag,lagw);
else
    if(verbose); disp('Initial Alignment & Polarity Given as Input'); end
    dt=info(:,1);
    pol=info(:,2);
end

% refine alignment estimate
% - selects lag with minimum weighted misfit to initial alignment among the
%   peaks with the same polarity as the initial estimate
% - reinvert using with these more consistent peaks
% - repeat until max number of iterations or converged

% now get weighted errors based on consistency of lags to estimate


end
