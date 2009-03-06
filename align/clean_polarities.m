function [A,D,OK]=clean_polarities(B,method)
%
%
%    Description: [A,D,OK]=CLEAN_POLARITIES(B,METHOD) returns an orthogonal
%    square sign matrix A that satisfies A=v'*v where v is a row vector of
%    signs (1 or -1) under the stipulation that A is a minimally cleaned
%    version of B (an orthogonal square sign matrix with defects) with the
%    relationship A=B.*D.  METHOD is either 1 or 2.  OK is either true or
%    false and indicates if the solution is valid.
%    

sz=size(B);
if(~isnumeric(B) || numel(sz)~=2 || sz(1)~=sz(2) || ~isequal(B,B') || ~all(B(:)==1 | B(:)==-1))
    error('fuck:you','Input matrix must be a numeric orthogonal square sign matrix');
end

switch method
    % iterative average
    % get the matrix for each row
    % average those to get a new matrix (sign it)
    % iterate until convergence
    case 1
        A=zeros(sz(1));
        C=B;
        while(~isequal(A,C))
            A=C;
            C=sum(repmat(permute(A,[2 3 1]),[1 sz(1) 1]).*repmat(permute(A,[3 2 1]),[sz(1) 1 1]),3);
            C=sign(C);
        end
    case 2
        % iterative difference
        % get the matrix for each row
        % sum difference with B to get defect matrix (sign it)
        % iterate until convergence
        A=zeros(sz(1));
        C=B;
        while(~isequal(A,C))
            A=C;
            C=sum(repmat(permute(A,[2 3 1]),[1 sz(1) 1]).*repmat(permute(A,[3 2 1]),[sz(1) 1 1]).*repmat(B,[1 1 sz(1)]),3);
            C=B.*sign(C);
        end
    %case 3
        % mode
    otherwise
        error('what:the:hell','Unknown method');
end

% it ain't perfect
OK=1;
if(rank(A)~=1 || any(A(:)==0))
    OK=0;
    warning('damn:it','Could not solve for polarities!');
end

D=B.*A;

end