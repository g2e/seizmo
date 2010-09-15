function [var]=nanvariance(x,biased,dim,w)
%NANVARIANCE    Return variance excluding NaNs
%
%    Usage: var=nanvariance(x)
%           var=nanvariance(x,biased)
%           var=nanvariance(x,biased,dim)
%           var=nanvariance(x,biased,dim,w)
%
%    Description:
%     VAR=NANVARIANCE(X) calculates the variance of X ignoring
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
%               SUM(W(NN),DIM)*SUM(W(NN).*ABS((X(NN)-X0(NN))).^2,DIM)
%        VAR = _______________________________________________________
%                        SUM(W(NN),DIM)^2-SUM(W(NN)^2,DIM)
%
%     where X0 is the weighted mean and NN are the non-NaN element indices.
%     Note the usage of ABS, which assures a real-valued variance.  Also
%     note that elements with a weight on NaN are ignored.
%
%    Notes:
%     - nanvar incompatibilities:
%       - weights are always the 4th argument in nanvariance while they are
%         the second input for nanvar (weighted variance in nanvar is
%         always biased)
%       - variances of 1 sample is always biased in nanvar (returns 0),
%         while nanvariance does an unbiased variance (returns nan) unless
%         the biased flag is set
%       - the weighting matrix, w, must be equal sized with x in
%         nanvariance, while nanvar requires a vector equal in length to
%         the dimension of x that the variance is being calculated on
%         (nanvar replicates the weights for array cases)
%
%    Examples:
%     % Unbiased variance down the 3rd dimension
%     % of a random valued matrix with some NaNs:
%     x=rand(10,10,10);
%     x(x>0.4 & x<0.6)=nan;
%     var=nanvariance(x,[],3)
%
%    See also: VAR, NANMEAN, NANVAR

%     Version History:
%        Oct. 14, 2009 - brought back, added checks, better documentation
%        Oct. 15, 2009 - working version (checks out with nanvar), dropped
%                        nanmean calls
%        Apr. 28, 2010 - appropriate formula for unbiased weighted variance
%        June  4, 2010 - fixed bug in unbiased weights
%        Sep. 13, 2010 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 13, 2010 at 14:45 GMT

% todo:

% CHECK NARGIN
error(nargchk(1,4,nargin));

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
    % NaNs
    nans=isnan(x);
    nne=sum(~nans,dim);
    
    % WEIGHTED MEAN
    x(nans)=0;
    x0=sum(x,dim)./nne;
    
    % GET RESIDUALS
    % --> NOTE THAT NANMEAN IS USED TO GET THE MEAN
    resid=x-x0(idx{:});
    
    % REMOVE INFLUENCE OF NaNs BY SETTING THEM TO 0
    resid(nans)=0;
    
    % UNBIASED VARIANCE
    if(~biased); nne=nne-1; end
    
    % RETURN NaN WHEN DOING UNBIASED VARIANCE WITH ONLY ONE ELEMENT
    % --> NOTE THAT THE MEAN RETURNS NaN IF THERE ARE ZERO ELEMENTS, SO
    %     THE VARIANCE WILL ALREADY BE NaN IN THESE CASES.
    % --> THIS IS ONLY NECESSARY TO AVOID DIVIDE BY ZERO WARNINGS
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
    
    % NaNs
    nans=isnan(x) | isnan(w);
    
    % GET WEIGHTED DEGREES OF FREEDOM (IGNORING NaNs)
    w(nans)=0;
    wnne=sum(w,dim);
    w2nne=sum(w.^2,dim);
    
    % WEIGHTED MEAN
    x(nans)=0;
    x0=sum(w.*x,dim)./wnne;
    
    % GET RESIDUALS
    resid=x-x0(idx{:});
    
    % REMOVE INFLUENCE OF NaNs BY SETTING THEM TO 0
    resid(nans)=0;
    
    % UNBIASED WEIGHTED VARIANCE
    if(~biased); wnne=(wnne.^2-w2nne)./wnne; end
    
    % RETURN NaN WHEN DOING UNBIASED VARIANCE WITH ONLY ONE ELEMENT
    % --> NOTE THAT THE MEAN RETURNS NaN IF THERE ARE ZERO ELEMENTS, SO
    %     THE VARIANCE WILL ALREADY BE NaN IN THESE CASES.
    % --> THIS IS ONLY NECESSARY TO AVOID DIVIDE BY ZERO WARNINGS
    wnne(wnne==0)=nan;
    
    % GET VARIANCE
    % --> NOTE THAT ABS FORCES REAL-VALUED OUTPUT
    var=sum(w.*abs(resid).^2,dim)./wnne;
end

end
