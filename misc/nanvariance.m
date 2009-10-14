function [var]=nanvariance(x,biased,dim,w)
%NANVARIANCE    Return variance excluding NaNs
%
%    Usage: var=nanvariance(x)
%           var=nanvariance(x,biased)
%           var=nanvariance(x,biased,dim)
%           var=nanvariance(x,biased,dim,w)
%
%    Description: VAR=NANVARIANCE(X) calculates the variance of X ignoring
%     NaN elements.  The variance(s) are taken across the first non-
%     singleton dimension and VAR is equal in size to X except for the
%     dimension that the variance is found for, which has a size of 1.
%
%     VAR=NANVARIANCE(X,BIASED) sets whether the variance is biased or
%     unbiased.  If BIASED is TRUE, the variance is a biased variance (no
%     degree of freedom removed for using the mean).  BIASED set FALSE (the
%     default) accounts for using the mean in the variance calculation.
%
%     VAR=NANVARIANCE(X,BIASED,DIM) specifies the dimension DIM that the
%     variance is taken across.  Use [] to get the default behavior.
%
%     VAR=NANVARIANCE(X,BIASED,DIM,W) uses weights W to calculate a
%     weighted variance.  W must be equal in size to X.  The formula for
%     the unbiased weighted variance is as follows:
%
%               SUM(ABS(W(NN)*(X(NN)-NANMEAN(X))).^2,DIM)
%        VAR = ________________________________________
%                SUM(W(NN),DIM)-(SUM(W(NN),DIM)/NNN)
%
%     where NN are the non-NaN element indices and NNN is the number of
%     non-NaN elements.  Note the ABS which assures a real-valued variance.
%     Also note that elements with a weight on NaN are ignored.
%
%    Notes:
%     - there is no special case for the first input of []
%     - nanvariance is NOT compatible nanvar!
%
%    Examples:
%     Mean and unbiased variance down the 3rd dimension of a random valued
%     matrix with some NaNs:
%      x=rand(10,10,10);
%      x(x>0.4 & x<0.6)=nan;
%      mean=nanmean(x,3);
%      var=nanvariance(x,[],3);
%
%    See also: VAR, NANMEAN, NANVAR

%     Version History:
%        Oct. 14, 2009 - brought back, added checks, better documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 14, 2009 at 05:25 GMT

% todo:

% CHECK NARGIN
msg=nargchk(1,4,nargin);
if(~isempty(msg)); error(msg); end

% CHECK X
if(~isnumeric(x) && ~islogical(x))
    error('seizmo:nanvariance:badInput','X must be numeric or logical!');
end

% CHECK BIASED
if(nargin<2 || isempty(biased)); biased=false; end
if(~isscalar(biased) || (~isnumeric(biased) && ~islogical(biased)))
    error('seizmo:nanvariance:badInput','BIASED must be TRUE or FALSE!');
end

% CHECK DIM
if(nargin<3 || isempty(dim))
    % SPECIAL BEHAVIOR TO MATCH NANVAR
    if(isequal(x,[])); var=nan(class(x)); return; end
    
    % DIMENSION TO WORK ON
    dim=find(size(x)~=1,1);
    if(isempty(dim)); dim=1; end
elseif(~isscalar(dim) || dim~=fix(dim) || dim<=0)
    error('seizmo:nanvariance:badInput','DIM must be a valid dimension!');
end

% INDEXING TO EXPAND MEAN TO SIZE OF X
idx(1:max(ndims(x),dim))={':'};
idx{dim}=ones(size(x,dim),1);

% UNWEIGHTED VARIANCE
if(nargin<4 || isempty(w))
    % GET RESIDUALS
    % --> NOTE THAT NANMEAN IS USED TO GET THE MEAN
    mean=nanmean(x,dim);
    resid=x-mean(idx{:});
    
    % NON-NaN ELEMENTS
    nans=isnan(resid);
    nne=sum(~nans,dim);
    
    % REMOVE INFLUENCE OF NaNs BY SETTING THEM TO 0
    resid(nans)=0;
    
    % UNBIASED VARIANCE
    if(~biased); nne=nne-1; end
    
    % RETURN NaN WHEN DOING UNBIASED VARIANCE WITH ONLY ONE ELEMENT
    % --> NOTE THAT NANMEAN RETURNS NaN IF THERE ARE ZERO ELEMENTS, SO
    %     THE VARIANCE WILL ALREADY BE NaN IN THESE CASES.
    nne(nne<1)=nan;
    
    % GET VARIANCE
    % --> NOTE THAT ABS FORCES REAL-VALUED OUTPUT
    var=sum(abs(resid).^2,dim)./nne;
% WEIGHTED VARIANCE
else
    % CHECK W
    if(~isequal(size(x),size(w)))
        error('seizmo:nanvariance:badInput','W must equal size of X!');
    end
    
    % GET WEIGHTED RESIDUALS
    % --> NOTE THAT NANMEAN IS USED TO GET THE MEAN
    mean=nanmean(x,dim);
    wresid=(x-mean(idx{:})).*w;
    
    % NON-NaN ELEMENTS
    nans=isnan(wresid);
    nne=sum(~nans,dim);
    
    % REMOVE INFLUENCE OF NaNs BY SETTING THEM TO 0
    wresid(nans)=0;
    
    % GET WEIGHTED DEGREES OF FREEDOM (IGNORING NaNs)
    w(nans)=0;
    wnne=sum(w,dim);
    
    % UNBIASED WEIGHTED VARIANCE
    if(~biased); wnne=wnne-wnne./nne; end
    
    % RETURN NaN WHEN DOING UNBIASED VARIANCE WITH ONLY ONE ELEMENT
    % --> NOTE THAT NANMEAN RETURNS NaN IF THERE ARE ZERO ELEMENTS, SO
    %     THE VARIANCE WILL ALREADY BE NaN IN THESE CASES.
    wnne(wnne==0)=nan;
    
    % GET VARIANCE
    % --> NOTE THAT ABS FORCES REAL-VALUED OUTPUT
    var=sum(abs(wresid).^2,dim)./wnne;
end

end
