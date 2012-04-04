function [pol,def,ok]=ttpolar(pg)
%TTPOLAR    Solves for polarities given relative polarities
%
%    Usage:  [pol,defects,ok]=ttpolar(pg)
%
%    Description:
%     [POL,DEFECTS,OK]=TTPOLAR(PG) returns estimated polarities POL that
%     produce the input matrix of relative polarities PG. PG & POL are
%     related by:
%         PG=(POL*POL').*DEFECTS
%     where PG is an NxN matrix, POL is an Nx1 column vector, and DEFECTS
%     is an NxN matrix.  All contain only 1s and -1s.  DEFECTS indicates
%     the relative polarities in PG that are not matched by POL.  The
%     solution is the polarities that minimize the number of defects.  In
%     some situations the relative polarities will balance exactly equal so
%     that there is no solution for some polarities (they will be set as
%     0).  In such cases a warning is issued and OK is set FALSE.
%
%    Notes:
%     - The problem is posed as an eigenvalue/eigenvector case.  We
%       seek a vector POL that produces a rank 1 matrix POL*POL' that is
%       minimally different from the relative polarity matrix PG.  The
%       eigenvector of PG with the highest eigenvalue is the most
%       significant component of PG and it's signs may be used as a
%       solution for POL.
%
%    Examples:
%     % Correlate records, solve for alignment & polarities, plot results:
%     plot0(data);
%     xc=correlate(data,'npeaks',1);
%     dt=ttalign(xc.lg);
%     data=timeshift(data,dt);
%     pol=ttpolar(xc.pg);
%     data=multiply(data,pol);
%     plot0(data);
%
%    See also: TTALIGN, TTSTDERR, TTREFINE, TTSOLVE

%     Version History:
%        Mar. 11, 2010 - initial version (from clean_polarities)
%        Sep. 13, 2010 - doc update, highest abs eigenvalue, nargchk fix
%        Jan. 23, 2011 - flag/warning for edge case
%        Apr.  2, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 11:10 GMT

% todo:
% - it would really be nice if this was a weighted inversion, then it could
%   done in conjunction with the arrival determination

% check nargin
error(nargchk(1,1,nargin));

% check input
pg=v2m(pg,'PG');

% get signs of eigenvector with largest eigenvalue
[v,d]=eig(pg);
[d,idx]=max(abs(diag(d)));

% this is inexact but works well
pol=sign(v(:,idx));

% check/fix polarity of zero
ok=true;
if(any(~pol))
    warning('seizmo:ttpolar:noSolution',...
        ['Polarity can not be solved for records:\n' ...
        sprintf('%g ',find(~pol))]);
    ok=false;
end

% return defect if necessary
if(nargout>1)
    def=(pol*pol').*pg;
end

end

function [x,nr]=v2m(x,str)
%V2M    Checks array and converts to matrix form if necessary

% allow either vector or matrix form
if(isvector(x))
    len=numel(x);                         % NUMBER OF LAGS
    nr=ceil(sqrt(2*len));                 % NUMBER OF RECORDS (FAST)
    %nr=round(max(roots([1 -1 -2*len]))); % NUMBER OF RECORDS (OLD)
    
    % assure length is ok
    if(~isreal(x) || (nr^2-nr)/2~=len || any(abs(x(:))~=1))
        error('seizmo:ttpolar:badInput',...
            '%s is not a properly lengthed vector of 1s & -1s!',str);
    end
    
    % convert to matrix form (with 1s down the middle)
    x=ndsquareform(x)+eye(nr);
else % matrix/grid form
    % number of records
    nr=size(x,1);

    % check input
    if(~isreal(x) || ~isequal(size(x),[nr nr]) ...
            || any(abs(x(:))~=1) || ~isequal(x,x'))
        error('seizmo:ttpolar:badInput',...
            '%s must be a square symmetric matrix of 1s & -1s!',str);
    end
end

end
