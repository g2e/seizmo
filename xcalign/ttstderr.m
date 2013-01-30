function [std]=ttstderr(m,lag,lagw)
%TTSTDERR    Weighted standard error of arrival times
%
%    Usage:    std=ttstderr(m,lag)
%              std=ttstderr(m,lag,lagw)
%
%    Description:
%     STD=TTSTDERR(M,LAG) returns the standard error of the arrival times
%     in M to the lags in LAG.  M is expected to be formatted like the
%     output of TTALIGN.  Note that this error measurement is more of a
%     measure of the inconsistency of the lags rather than an accurate
%     assessment of the true arrival time error.  LAG must be a NxN matrix
%     where N is the number of signals being compared or a N*(N-1)/2
%     element vector.
%
%     STD=TTSTDERR(M,LAG,LAGW) allows weighting the standard error solution
%     using LAGW.  LAGW and LAG should be equal size.
%
%    Notes:
%     - SNR2MAXPHASEERROR also provides a measure of the phase error
%       although even this has assumptions (narrowbandedness and that the
%       waveforms being correlated have strong similarity).
%     - Diagonal elements of LAG are ignored because they are assumed to
%       correspond to autocorrelations.  These always have 0 lag and 0
%       residual and are not useful error measures.
%
%    Examples:
%     % Get arrival time solution and standard error:
%     m=ttalign(lag);
%     mstd=ttstderr(m,lag);
%
%    See also: TTALIGN, TTPOLAR, TTREFINE, TTSOLVE, SNR2MAXPHASEERROR

%     Version History:
%        Mar.  2, 2010 - initial version (from dtwresid)
%        Mar. 11, 2010 - doc update, column vector output
%        Mar. 12, 2010 - better way to get column vector output, more
%                        input checking
%        Mar. 22, 2010 - account for ttalign fix, ignore diagonal elements
%                        in error calculation (makes error slightly higher)
%        Sep. 13, 2010 - nargchk fix
%        Apr.  2, 2012 - minor doc update
%        Jan. 30, 2013 - doc update, better checking, require 2D lag input
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 30, 2013 at 01:05 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% check inputs
if(~isreal(m) || ~isvector(m))
    error('seizmo:ttstderr:badInput',...
        'M must be a real-valued vector!');
elseif(~isreal(lag))
    error('seizmo:ttstderr:badInput',...
        'LAG must be a real-valued array!');
elseif(ndims(lag)>2)
    error('seizmo:ttstderr:badInput',...
        'LAG must be a 2D matrix or vector!');
end

% vector to grid
if(isvector(lag))
    lag=ndsquareform(lag,'tomatrix',false);
else
    if(~isequal(lag,-lag.'))
        error('seizmo:ttstderr:badInput',...
            'LAG must be a anti-symmetric matrix of real values!');
    end
end

% number of records
nr=size(lag,1);
if(nr~=numel(m))
    error('seizmo:ttstderr:badInput',...
        'M & LAG are inconsistent in size!');
end

% handle weights
if(nargin==2 || isempty(lagw)); lagw=ones(size(lag)); end
if(isscalar(lagw)); lagw=lagw(ones(size(lag))); end
if(~isreal(lagw))
    error('seizmo:ttstderr:badInput',...
        'LAGW must be a real-valued array!');
elseif(ndims(lagw)>2)
    error('seizmo:ttstderr:badInput',...
        'LAGW must be a 2D matrix or vector!');
end

% vector to grid
if(isvector(lagw))
    lagw=ndsquareform(lagw,'tomatrix');
else
    if(~isequal(lagw,lagw.'))
        error('seizmo:ttstderr:badInput',...
            'LAGW must be a symmetric matrix of real values!');
    end
end

% number of records
nr=size(lagw,1);
if(nr~=numel(m))
    error('seizmo:ttstderr:badInput',...
        'LAG & LAGW are inconsistent in size!');
end

% force m to be column vector
m=m(:);

% ignore diagonal
lag=lag+diag(nan(nr,1));

% get standard error of arrivals
std=sqrt(nanvariance(lag-m(:,ones(nr,1))+m(:,ones(nr,1)).',0,2,lagw));

end
