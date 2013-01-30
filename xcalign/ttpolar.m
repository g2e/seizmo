function [pol,def,ok]=ttpolar(pg)
%TTPOLAR    Solves for polarities given relative polarities
%
%    Usage:    [pol,defects,ok]=ttpolar(pg)
%
%    Description:
%     [POL,DEFECTS,OK]=TTPOLAR(PG) returns estimated polarities POL that
%     produce the input matrix of relative polarities PG. PG & POL are
%     related by:
%         PG=(POL*POL').*DEFECTS
%     where PG is an NxN matrix, POL is an Nx1 column vector, and DEFECTS
%     is an NxN matrix.  All contain only 1s and -1s.  DEFECTS indicates
%     the relative polarities in PG that are not matched by POL (usually
%     due to cycle skipping).  The solution is the polarities that minimize
%     the number of defects.  In some situations there is not a unique
%     solution some entries in POL because the relative polarities are in
%     balance.  In such cases a warning is issued, OK is set FALSE, & the
%     corresponding POL entries are set to 0.
%
%    Notes:
%     - The solution is found by posing the problem as an eigenvalue/
%       eigenvector decomposition.  We seek a vector POL that produces a
%       rank 1 matrix POL*POL' that is minimally different from the
%       relative polarity matrix PG.  The eigenvector of PG with the
%       highest eigenvalue is the most significant component of PG and the
%       signs of the elements may be used as a solution for POL.
%
%    Examples:
%     % Correlate records, solve for alignment & polarities, plot results:
%     plot0(data);
%     xc=correlate(data,'mcxc','noa','a','p',{'npeaks',1});
%     dt=ttalign(xc.lg);
%     data1=timeshift(data,-dt);
%     pol=ttpolar(xc.pg);
%     data1=multiply(data1,pol);
%     plot0(data1);
%
%    See also: TTALIGN, TTSTDERR, TTREFINE, TTSOLVE

%     Version History:
%        Mar. 11, 2010 - initial version (from clean_polarities)
%        Sep. 13, 2010 - doc update, highest abs eigenvalue, nargchk fix
%        Jan. 23, 2011 - flag/warning for edge case
%        Apr.  2, 2012 - minor doc update
%        Jan. 29, 2013 - doc update for new correlate functions
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 29, 2013 at 11:10 GMT

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
