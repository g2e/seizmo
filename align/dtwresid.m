function [S]=dtwresid(DT,CG,LG,xcpow)
%DTWRESID    Gets weighted standard deviations of inverted arrival times 
%            from the misfit to xc lags
%
% INPUTS REQUIRED:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DT - calculated relative arrival times
% CG - grid of xc correlation values (used for weighting)
% LG - grid of xc lag times
% xcpow - power(s) to apply to xc correlation values for weighting (0=no xc value weight)

% NUMBER OF RECORDS
nr=size(DT,1);
if(any(size(CG)~=nr) || any(size(LG)~=nr)); error('size mismatch'); end

S=sqrt(nanwvar(LG+DT(:,ones(nr,1))-DT(:,ones(nr,1)).',0,1,CG.^xcpow));

end
