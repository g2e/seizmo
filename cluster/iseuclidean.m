function t = iseuclidean(D)
%ISEUCLIDEAN Is a distance matrix Euclidean?
%   T = ISEUCLIDEAN(D) returns a logical indicating whether or not the
%   dissimilarity matrix D is a Euclidean distance matrix, i.e., whether
%   there exist n points in p-dimensional space (for some p < n) such that
%   their Euclidean distances are given by D.  D may be specified as either
%   a full (square, symmetric) dissimilarity matrix, or as the lower
%   triangle (e.g., output by PDIST).
%
%   This algorithm is essentially classical multidimensional scaling.
%
%   See also CMDSCALE, PDIST, LINKAGE.
%
%   References:
%     [1] Seber, G.A.F., Multivariate Observations, Wiley, 1984

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.4.2.4 $ $Date: 2006/11/11 22:57:32 $

[n m] = size(D);
del = 10*eps(class(D));

% lower triangle form for D
if n == 1
    % make sure it's a valid dissimilarity matrix
    n = ceil(sqrt(2*m)); % (1+sqrt(1+8*m))/2, but works for large m
    if n*(n-1)/2 == m & all(D >= 0)
        D = squareform(D);
    else
        warning('stats:iseuclidean:NotDistanceMatrix',...
                'Not a valid dissimilarity or distance matrix.')
        t = logical(0);
        return
    end
    
% full matrix form, make sure it's valid dissimilarity matrix
elseif n ~= m | any(any(D < 0 | abs(D - D') > del*max(max(D)))) | any(diag(D) > del)
    warning('stats:iseuclidean:NotDistanceMatrix',...
            'Not a valid dissimilarity or distance matrix.')
    t = logical(0);
    return
end

P = eye(n) - repmat(1/n,n,n);
B = P * (-.5 * D .* D) * P;
g = eig((B+B')./2); % guard against spurious complex e-vals from roundoff
t = all(-eps(class(g))^(3/4) * max(abs(g)) <= g); % all non-negative eigenvals (within roundoff)?
