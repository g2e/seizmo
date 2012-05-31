function [med]=nanmedian(x,dim)
%NANMEDIAN    Return median excluding NaNs
%
%    Usage: med=nanmedian(x)
%           med=nanmedian(x,dim)
%
%    Description:
%     MED=NANMEDIAN(X) returns the medians along the first non-
%     singleton dimension of X excluding NaN elements.  MEAN is equal in
%     size to X except for the dimension that the mean is computed for,
%     which has a size of 1.
%
%     MED=NANMEDIAN(X,DIM) specifies the dimension DIM that the median is
%     taken across.  Use [] to get the default behavior.
%
%    Notes:
%     - Equivalent to NANMEDIAN of the Statistics Toolbox.
%
%    Examples:
%     % Median down the 3rd dimension of a
%     % random valued matrix with some NaNs:
%     x=rand(3,4,5);
%     x(x>0.4 & x<0.6)=nan;
%     med=nanmedian(x,3)
%
%    See also: MEDIAN, NANMEAN, NANVARIANCE

%     Version History:
%        May   4, 2012 - initial version
%        May  15, 2012 - fixed reshaping bug, further optimization
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  15, 2012 at 14:45 GMT

% CHECK NARGIN
error(nargchk(1,2,nargin));

% IDENTIFY NaN ELEMENTS
nans=isnan(x);

% PASS ON TO MEDIAN IF NO NANS
if(~any(nans(:)))
    if(nargin==1); med=median(x);
    else med=median(x,dim);
    end
else % NANS EXIST SO GO COLUMN BY COLUMN
    sx=size(x);
    
    % GET DIMENSION TO WORK ON
    if(nargin==1 || isempty(dim))
        dim=find(sx~=1,1);
        if(isempty(dim)); dim=1; end
    end
    
    % PERMUTE DIMENSION TO FRONT
    x=permute(x,[dim:max(numel(sx),dim) 1:dim-1]);
    if(dim>numel(sx)); sx=[sx ones(1,dim-numel(sx))]; end
    sx=sx([dim:max(numel(sx),dim) 1:dim-1]);
    
    % MAKE A MATRIX
    nr=sx(1); nc=prod(sx)/nr;
    x=reshape(x,nr,nc);
    
    % SETUP
    x=sort(x,1);
    nn=~isnan(x);
    nnc=all(nn);
    med=nan(1,nc,class(x));
    
    % PROCESS NON-NAN COLUMNS SEPARATELY
    med(nnc)=median(x(:,nnc));
    for i=find(~nnc)
        last=find(nn(:,i),1,'last');
        if(last)
            if(mod(last,2)); med(i)=x((last+1)/2,i);
            else med(i)=(x(last/2,i)+x(last/2+1,i))/2;
            end
        end
    end
    
    % RESIZE
    sx(1)=1;
    med=reshape(med,sx);
    
    % UNDO PERMUTE
    med=ipermute(med,[dim:max(numel(sx),dim) 1:dim-1]);
end

end
