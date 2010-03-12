function [pol,def]=ttpolar(pg)
%TTPOLAR    Solves for polarities given relative polarities
%
%    Usage:  [pol,defects]=ttpolar(pg)
%
%    Description: [POL,DEFECTS]=TTPOLAR(PG) returns the likely polarities
%     POL that produced the matrix of relative polarities PG.  PG & POL are
%     related by:
%                           PG=(POL*POL').*DEFECTS
%     where PG is an NxN matrix, POL is an Nx1 column vector, and DEFECTS
%     is an NxN matrix.  All contain only 1s and -1s.  DEFECTS indicates
%     the relative polarities in PG that are not matched by POL.
%
%    Notes:
%     - solved as an eigenvalue/eigenvector problem
%
%    Examples:
%     Correlate records, solve for alignment & polarities, plot results:
%      plot0(data);
%      xc=correlate(data,'npeaks',1);
%      dt=ttalign(xc.lg);
%      data=timeshift(data,dt);
%      pol=ttpolar(xc.pg);
%      data=multiply(data,pol);
%      plot0(data);
%
%    See also: TTALIGN, TTSTDERR, TTREFINE, TTSOLVE

%     Version History:
%        Mar. 11, 2010 - initial version (from clean_polarities)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 11, 2010 at 11:10 GMT

% todo:
% - it would really be nice if this was a weighted inversion

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check input
pg=v2m(pg,'PG');

% get signs of eigenvector with largest eignevalue
[v,d]=eig(pg);
[d,idx]=max(diag(d));
pol=sign(v(:,idx));

% return defect if necessary
if(nargout==2)
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
