function [d]=gcdn(x,dim)
%GCDN    Returns greatest common devisor for a vector of numbers
%
%    Usage:    d=gcdn(x)
%              d=gcdn(x,dim)
%
%    Description:
%     D=GCDN(X) returns the greatest common devisor of the elements of X if
%     X is a vector.  If X is an array GCDN computes the greatest common
%     devisor across the first non-singleton dimension.  For a 2D array
%     this would be across the rows, resulting in a separate value for each
%     column.
%
%     D=GCDN(X,DIM) indicates the dimension that the greatest common
%     devisor is computed across.  The default is the first non-singleton
%     dimension.
%
%    Notes:
%
%    Examples:
%     % GCD of 20, 30, & 45
%     gcd([20 30 45])
%
%     % GCD down the pages of a 6x6x6 array of linear indices:
%     gcdn(lind(6,6,6),3)
%
%    See also: GCD

%     Version History:
%        Aug.  8, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  8, 2013 at 17:45 GMT

% todo:

% check x
if(~isnumeric(x) || ~isreal(x))
    error('seizmo:gcdn:badX','X must be a real-valued numeric array!');
elseif(numel(x)<2)
    error('seizmo:gcdn:badX','X must have 2 or more elements!');
elseif(any(x(:)~=fix(x(:))) || any(isinf(x(:))))
    error('seizmo:gcdn:badX','X must contain only integers!');
end

% class check x
cx=class(x);
switch cx
    case 'double'
        bad=2^53-1;
    case 'single'
        bad=2^24-1;
end
if(any(x(:)>bad))
    warning('seizmo:gcdn:badPrecision',...
        ['Values in X exceed floating point integer precision level!\n' ...
        'Output is likely incorrect for these values!']);
end

% default dimension
sx=size(x);
ndx=numel(sx);
if(nargin==1 || isempty(dim))
    dim=find(sx~=1,1);
end

% check dimension
if(~isscalar(dim) || dim~=fix(dim))
    error('seizmo:gcdn:badDim','DIM must be a scalar integer!');
elseif(ndx<dim || sx(dim)==1)
    error('seizmo:gcdn:badDim','X dimension %d must have size>=2!',dim);
end

% abs & permute
x=permute(abs(x),[dim 1:dim-1 dim+1:ndx]);
sx1=sx(dim);
sx(dim)=[];

% loop over each "column"
d=nan([1 sx],cx);
d(:)=x(1,:); % fill with 1st "row"
for i=1:prod(sx)
    % loop down the vector
    for j=2:sx1
        while(x(j,i))
            tmp=d(1,i)-x(j,i)*floor(d(1,i)/x(j,i));
            d(1,i)=x(j,i);
            x(j,i)=tmp;
        end
    end
end

% ipermute
d=ipermute(d,[dim 1:dim-1 dim+1:ndx]);

end
