function [x]=ndsquareform(x,method,flag)
%NDSQUAREFORM    Reshapes between an n-d distance matrix and "vector"
%
%    Usage:    m=ndsquareform(v)
%              v=ndsquareform(m)
%              m=ndsquareform(v,'tomatrix')
%              v=ndsquareform(m,'tovector')
%              m=ndsquareform(v,method,flag)
%              v=ndsquareform(m,method,flag)
%
%    Description: V=NDSQUAREFORM(M) & M=NDSQUAREFORM(V) reshapes M & V
%     between an n-dimensional, symmetric, square matrix M and the
%     corresponding n-d "vector" V (quoted because it isn't quite a vector,
%     as it may have >1 non-singletone dimensions - but at least one of the
%     first two must be singleton).  The vector form is just the upper
%     triangle section of the matrix form and thus consumes less space.
%     The disadvantage to the vector form is that one does not immediately
%     know where a specific element comes from in the matrix form.  Vectors
%     are always returned as column-vectors.  Matrices are always returned
%     with the diagonal elements set to zero.
%
%     M=NDSQUAREFORM(V,'TOMATRIX') requires NDSQUAREFORM to interpret the
%     input as a vector to be transformed into a matrix.  This option
%     should never be necessary except for readibility.  Leaving this field
%     empty sets NDSQUAREFORM to auto-detect (the default).
%
%     V=NDSQUAREFORM(M,'TOVECTOR') requires NDSQUAREFORM to interpret the
%     input as a matrix to be transformed into a vector.  This is necessary
%     if the matrix is 1x1x..., as this will be auto-detected as a vector.  
%     Leaving this field empty sets NDSQUAREFORM to auto-detect (the
%     default).
%
%     M=NDSQUAREFORM(V,METHOD,FLAG) allows switching between creation of a
%     symmetric and an anti-symmetric matrix.  Setting FLAG to TRUE (the
%     default) converts V to an n-d symmetric, square matrix.  With FLAG
%     set to FALSE, V is converted to an n-d anti-symmetric, square matrix
%     (the lower triangle elements have opposite sign to the upper triangle
%     elements).  Remember that V always corresponds to the upper triangle
%     values in M!
%
%     V=NDSQUAREFORM(M,METHOD,FLAG) allows switching between extraction of
%     the upper and lower triangle portions of M.  When FLAG is set to TRUE
%     (the default), the upper triangle portion of M is returned in vector
%     form as V.  Setting FLAG to FALSE returns the lower triangle portion.
%
%    Note:
%     - 1x1x... is always assumed a vector unless METHOD is set to
%       'ismatrix'!
%     - diagonal is NOT required to be zeros
%
%    Example:
%      v=repmat(1:10,[1 1 2 2])
%      m=ndsquareform(v)
%      m(:,:,1,1)=[0  1  2  3  4
%                  1  0  5  6  7
%                  2  5  0  8  9
%                  3  6  8  0 10
%                  4  7  9 10  0];
%      m(:,:,2,1)=[0  1  2  3  4
%                  1  0  5  6  7
%                  2  5  0  8  9
%                  3  6  8  0 10
%                  4  7  9 10  0];
%      m(:,:,1,2)=[0  1  2  3  4
%                  1  0  5  6  7
%                  2  5  0  8  9
%                  3  6  8  0 10
%                  4  7  9 10  0];
%      m(:,:,2,2)=[0  1  2  3  4
%                  1  0  5  6  7
%                  2  5  0  8  9
%                  3  6  8  0 10
%                  4  7  9 10  0];
%      v=ndsquareform(m)
%
%    See also: squareform, tril, triu, permute

%     Version History:
%        Sep.  8, 2009 - rewrite and added documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  8, 2009 at 19:30 GMT

% todo:

% check nargin
msg=nargchk(1,3,nargin);
if(~isempty(msg)); error(msg); end

% size up input
sz=size(x);
ndims=numel(sz);
pages=max([1 prod(sz(3:end))]);

% default method
if(nargin<2 || isempty(method))
    % auto-detect
    method='tovector';
    if(any(sz(1:2)==1)); method='tomatrix'; end
end

% default flag
if(nargin<3 || isempty(flag)); flag=true; end

% proceed by method
switch lower(method)
    case 'tovector'
        % check squareness and symmetry (only first page)
        if(sz(1)~=sz(2) || (~isequal(x(:,:,1),x(:,:,1).') && ...
                ~isequal(x(:,:,1),-x(:,:,1).')))
            error('seizmo:ndsquareform:badInput',...
                'M must be a square, (anti-)symmetric matrix!');
        end
        
        if(flag)
            % swap upper triangle with lower to make extraction easy
            x=permute(x,[2 1 3:ndims]);
        end
        
        % logical indexing
        lti=tril(true(sz(1)),-1);
        
        % extract lower triangle, replicating logical indices to match
        % the number of pages (this preserves the n-dimensionality too)
        x=x(lti(:,:,ones(1,pages)));
    case 'tomatrix'
        % get size of output matrix
        len=prod(sz(1:2));
        n=ceil(sqrt(2*len));
        
        % check vectorness and proper length
        if(~any(sz(1:2)==1) || len~=(n^2-n)/2)
            error('seizmo:ndsquareform:badInput',...
                'V is not a ''vector'' or is an improper length!');
        end
        
        % preallocate output (preserving class)
        if(isnumeric(x))
            m=zeros([n n sz(3:end)],class(x));
        elseif(islogical(x))
            m=false([n n sz(3:end)]);
        else
            error('seizmo:ndsquareform:badInput',...
                'V must be logical or numeric!');
        end
        
        % logical indices of lower triangle
        lti=tril(true(n),-1);
        
        % fill lower triangle
        m(lti(:,:,ones(1,pages)))=x(:);
        
        % permute to upper and fill lower
        if(flag)
            % symmetric matrix
            x=m+permute(m,[2 1 3:ndims]);
        else
            % anti-symmetric matrix
            x=permute(m,[2 1 3:ndims])-m;
        end
    otherwise
        error('seizmo:ndsquareform:badMethod',...
            'METHOD must be either ''tovector'' or ''tomatrix''!');
end

end
