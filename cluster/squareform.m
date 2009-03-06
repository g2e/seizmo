function Z = squareform(Y,dir)
%SQUAREFORM Reformat a distance matrix between upper triangular and square form.
%   Z = SQUAREFORM(Y), if Y is a vector as created by the PDIST function,
%   converts Y into a symmetric, square format, so that Z(i,j) denotes the
%   distance between the i and j objects in the original data.
%
%   Y = SQUAREFORM(Z), if Z is a symmetric, square matrix with zeros along
%   the diagonal, creates a vector Y containing the Z elements below the
%   diagonal.  Y has the same format as the output from the PDIST function.
%
%   Z = SQUAREFORM(Y,'tovector') forces SQUAREFORM to treat Y as a vector.
%   Y = SQUAREFORM(Z,'tomatrix') forces SQUAREFORM to treat Z as a matrix.
%   These formats are useful if the input has a single element, so it is
%   ambiguous as to whether it is a vector or square matrix.
%
%   Example:  If Y = (1:6) and X = [0  1  2  3
%                                   1  0  4  5
%                                   2  4  0  6
%                                   3  5  6  0],
%             then squareform(Y) is X, and squareform(X) is Y.
%
%   See also PDIST.

%   Copyright 1993-2007 The MathWorks, Inc. 
%   $Revision: 1.7.2.8 $

if ~(isnumeric(Y) || islogical(Y)) || ndims(Y) > 2
   error('stats:squareform:BadInput','Y must be a numeric or logical vector or matrix.');
end

[m, n] = size(Y);
if nargin<2 || isempty(dir)
   if isvector(Y)
      dir = 'tomatrix';
   else
      dir = 'tovector';
   end
end
okdirs = {'tovector' 'tomatrix'};
j = strmatch(dir,okdirs);
if ~isscalar(j)
   error('stats:squareform:BadDirection',...
         'Direction argument must be ''tomatrix'' or ''tovector''.');
end
dir = okdirs{j};

switch(dir)
 case 'tomatrix'
   if ~isvector(Y)
      error('stats:squareform:BadInput',...
            'Y must be a vector to convert it to a distance matrix');
   end
   if m~=1
      Y = Y';
      n = m;
   end

   m = ceil(sqrt(2*n)); % (1 + sqrt(1+8*n))/2, but works for large n
   if m*(m-1)/2 ~= n
      error('stats:squareform:BadInput',...
            'The size of the vector Y is not correct');
   end

   if islogical(Y)
      Z = false(m);
      if m>1
         Z(tril(true(m),-1)) = Y;
         Z = Z | Z';
      end
   else % isnumeric(Y)
      Z = zeros(m,class(Y));
      if m>1
         Z(tril(true(m),-1)) = Y;
         Z = Z + Z';
      end
   end

 case 'tovector'
   if m~=n || ~all(diag(Y)==0)
      error('stats:squareform:BadInput',...
            'The distance matrix Z must be square with 0 along the diagonal.');
   end

   Z = Y(tril(true(n),-1));
   Z = Z(:)';                 % force to a row vector, even if empty
end
