function [var]=nanvar(x,opt,dim,w)
%NANMEAN    Return variance excluding NaNs
%
%    Description:  Returns the variances along a dimension excluding NaN 
%     elements.
%     
%     Other differences from var: 
%       - weights for a weighted variance come after the dimension 
%         parameter so weighted variance can also be done unbiased style
%       - weight matrix must be the same size as the first input so every
%         element has its own weight
%       - there is no special case for the first input of []
%
%     Elements with a weight of NaN are ignored.  Weights with a 
%     corresponding element of NaN are ignored.
%
%    Usage: var=nanvar(x)
%           var=nanvar(x,opt)
%           var=nanvar(x,opt,dim)
%           var=nanvar(x,opt,dim,w)
%
%    Examples:
%
%    See also: var, nanmean, nanstd

% ASSURE INPUT IS FLOAT
if(isinteger(x)); x=double(x); end

% CHECK OPT (UNBIASED=0, BIASED=1)
if(nargin<2 || isempty(opt)); opt=0; end
if(~isscalar(opt) || ~any(opt==[0 1]))
    error('Second argument (bias option) must be 0 or 1'); 
end

% DIMENSION TO OPERATE ON
if(nargin<3 || isempty(dim))
    dim=find(size(x)~=1,1);
    if(isempty(dim)); dim=1; end
end

% TILING FOR MEAN (INPUT INTO REPMAT)
tile=ones(1,max(ndims(x),dim)); tile(dim)=size(x,dim);

% UNWEIGHTED VARIANCE
if(nargin<4 || isempty(w))
    % GET RESIDUALS
    % --> NOTE THAT NANMEAN IS USED TO GET THE MEAN
    resid=x-repmat(nanmean(x,dim),tile);
    
    % IDENTIFY AND COUNT NUMBER OF NON-NaN ELEMENTS
    nans=isnan(resid);
    nne=sum(~nans,dim);
    
    % SET NaN ELEMENTS TO ZERO (NO INFLUENCE ON SUM)
    resid(nans)=0;
    
    % UNBIASED VARIANCE
    if(~opt); nne=nne-1; end
    
    % RETURN NaN WHEN DOING UNBIASED VARIANCE WITH ONLY ONE ELEMENT
    % --> NOTE THAT NANMEAN RETURNS NaN IF THERE ARE ZERO ELEMENTS, SO
    %     THE VARIANCE WILL ALREADY BE NaN IN THESE CASES.
    nne(nne<1)=nan;
    
    % GET VARIANCE
    % --> NOTE THAT ABS RETURNS THE MAGNITUDE GUARANTEEING A REAL OUTPUT
    var=sum(abs(resid).^2,dim)./nne;
% WEIGHTED VARIANCE
else
    % CHECK WEIGHTING MATRIX
    if(~isequal(size(x),size(w)))
        error('Weight matrix must equal size of data matrix')
    end
    
    % ASSURE INPUT IS FLOAT
    if(isinteger(w)); w=double(w); end
    
    % GET WEIGHTED RESIDUALS
    % --> NOTE THAT NANMEAN IS USED TO GET THE MEAN
    wresid=(x-repmat(nanmean(x,dim),tile)).*w;
    
    % IDENTIFY AND COUNT NUMBER OF NON-NaN ELEMENTS
    nans=isnan(wresid);
    nne=sum(~nans,dim);
    
    % SET NaN ELEMENTS TO ZERO (NO INFLUENCE ON SUM)
    wresid(nans)=0;
    
    % GET WEIGHTED DEGREES OF FREEDOM
    w(nans)=0; wnne=sum(w,dim);
    
    % UNBIASED WEIGHTED VARIANCE
    if(~opt); wnne=wnne-wnne./nne; end
    
    % RETURN NaN WHEN DOING UNBIASED VARIANCE WITH ONLY ONE ELEMENT
    % --> NOTE THAT NANMEAN RETURNS NaN IF THERE ARE ZERO ELEMENTS, SO
    %     THE VARIANCE WILL ALREADY BE NaN IN THESE CASES.
    wnne(wnne==0)=nan;
    
    % GET VARIANCE
    % --> NOTE THAT ABS RETURNS THE MAGNITUDE GUARANTEEING A REAL OUTPUT
    var=sum(abs(wresid).^2,dim)./wnne;
end

end
