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
%    Description:
%     V=NDSQUAREFORM(M) & M=NDSQUAREFORM(V) reshapes M & V between an
%     n-dimensional, symmetric, square matrix M and the corresponding n-d
%     "vector" V (quoted because it isn't quite a vector, as it may have >1
%     non-singletone dimensions - think of it as an array of vectors).  The
%     vector form is the lower triangle section of the matrix form and thus
%     consumes less space.  The primary disadvantage to the vector form is
%     that one does not immediately know what a specific element
%     corresponds to in the matrix form.  Vectors are always returned as
%     row vectors.  Matrices are always returned with the diagonal elements
%     set to zero.
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
%     default).  V is a multi-page row vector.
%
%     M=NDSQUAREFORM(V,'TOMATRIX',FLAG) switches between creation of a
%     symmetric and an anti-symmetric matrix.  Setting FLAG to TRUE (the
%     default) converts V to an n-d symmetric, square matrix.  With FLAG
%     set to FALSE, V is converted to an n-d anti-symmetric, square matrix
%     (the lower triangle elements have opposite sign to the upper triangle
%     elements).  If V is a logical vector then the upper and lower
%     triangles have opposite logical values in the anti-symmetric case.
%     Remember that V always corresponds to the lower triangle values in M!
%
%     V=NDSQUAREFORM(M,'TOVECTOR',FLAG) switches between extraction of
%     the lower and upper triangle portions of M.  When FLAG is set to TRUE
%     (the default), the lower triangle portion of M is returned in vector
%     form as V.  Setting FLAG to FALSE returns the upper triangle section.
%
%    Note:
%     - 1x1x... is always assumed a vector unless METHOD is set to
%       'tovector'!
%     - The diagonal is NOT required to be zeros.
%
%    Example:
%     % Create an n-d vector and convert to matrix form:
%     v=repmat(1:10,[1 1 2 2])
%     m=ndsquareform(v)
%     % gives:
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
%     % and going back to vector form:
%     v=ndsquareform(m)
%     % gives:
%      v(:,:,1,1)=[1 2 3 4 5 6 7 8 9 10];
%      v(:,:,2,1)=[1 2 3 4 5 6 7 8 9 10];
%      v(:,:,1,2)=[1 2 3 4 5 6 7 8 9 10];
%      v(:,:,2,2)=[1 2 3 4 5 6 7 8 9 10];
%
%    See also: SQUAREFORM, TRIL, TRIU, PERMUTE

%     Version History:
%        Sep.  8, 2009 - rewrite and added documentation
%        Oct. 14, 2009 - (fully?) compatible with SQUAREFORM, note switch
%                        to lower triangle, lots of fixes
%        Mar.  2, 2010 - minor doc update
%        Mar. 12, 2010 - doc update
%        Mar. 22, 2010 - improved checking of anti-symmetric matrices
%        Feb. 11, 2011 - mass nargchk fix
%        Apr.  3, 2012 - minor doc update
%        Sep. 28, 2012 - fixed nasty eps usage bug
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 28, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

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
        xd=x((1:sz(1))+((1:sz(1))-1)*sz(1));
        if(sz(1)~=sz(2) || (~isnumeric(x) && ~islogical(x)))
            error('seizmo:ndsquareform:badInput',...
                'M must be a square numeric or logical array!');
        elseif(isnumeric(x) && (~isequal(x(:,:,1),x(:,:,1).') && ...
                any(any(((x(:,:,1)-diag(xd))-(diag(xd)-x(:,:,1).'))...
                >eps(max(max(x(:,:,1))))))))
            % the above uses diag to ignore the diagonal
            error('seizmo:ndsquareform:badInput',...
                'M must be a(n) (anti-)symmetric matrix!');
        elseif(islogical(x) && (~isequal(x(:,:,1),x(:,:,1).') && ...
                ~isequal(x(:,:,1) | diag(true(sz(1),1)), ...
                ~x(:,:,1).' | diag(true(sz(1),1)))))
            % the above uses diag to ignore the diagonal
            error('seizmo:ndsquareform:badInput',...
                'M must be a(n) (anti-)symmetric matrix!');
        end
        
        if(~flag)
            % swap upper triangle with lower
            x=permute(x,[2 1 3:ndims]);
        end
        
        % logical indexing
        lti=tril(true(sz(1)),-1);
        
        % extract lower triangle, replicating logical
        % indices to match the number of pages
        % - we also force a row vector here to match SQUAREFORM & PDIST
        x=reshape(x(lti(:,:,ones(1,pages))),...
            [1 (sz(1)^2-sz(1))/2 sz(3:end)]);
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
            lgc=false;
            m=zeros([n n sz(3:end)],class(x));
        elseif(islogical(x))
            lgc=true;
            m=false([n n sz(3:end)]);
        else
            error('seizmo:ndsquareform:badInput',...
                'V must be logical or numeric!');
        end
        
        % logical indices of lower triangle
        lti=tril(true(n),-1);
        
        % fill lower triangle
        m(lti(:,:,ones(1,pages)))=x(:);
        
        % fill upper
        if(flag)
            % symmetric matrix
            if(lgc)
                x=m | permute(m,[2 1 3:ndims]);
            else
                x=m+permute(m,[2 1 3:ndims]);
            end
        else
            % anti-symmetric matrix
            if(lgc)
                % preserve the diagonal & lower triangle
                m=permute(m,[2 1 3:ndims]);
                m(lti(:,:,ones(1,pages)))=~x(:);
                x=permute(m,[2 1 3:ndims]);
            else
                x=m-permute(m,[2 1 3:ndims]);
            end
        end
    otherwise
        error('seizmo:ndsquareform:badMethod',...
            'METHOD must be either ''tovector'' or ''tomatrix''!');
end

end
