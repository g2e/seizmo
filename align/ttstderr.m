function [std]=ttstderr(m,lag,lagw)
%TTSTDERR    Weighted standard error of arrival times
%
%    Usage:    stderr=ttstderr(m,lag)
%              stderr=ttstderr(m,lag,lagw)
%
%    Description: STDERR=TTSTDERR(M,LAG) returns the standard error of the
%     arrival times in M to the lags in LAG.  M is expected to be formatted
%     like the output of TTALIGN.  Note that this error measurement is more
%     of a measure of the inconsistency of the lags rather than an accurate
%     assessment of the true arrival time error.
%
%     STDERR=TTSTDERR(M,LAG,LAGW) allows weighting the standard error
%     solution using LAGW.  LAGW and LAG should correspond to one another.
%
%    Notes:
%     - For narrow band situations, SNR2MAXPHASEERROR is likely a better
%       measure of the error although even this has assumptions (that the
%       waveforms being correlated have a strong similarity).
%     - Diagonal elements of LAG are ignored because they are assumed to
%       correspond to autocorrelations.  These always have 0 lag and 0
%       residual and are not useful error measures.
%
%    Examples:
%     Get arrival time solution and standard error:
%      m=ttalign(lag);
%      mstd=ttstderr(m,lag);
%
%    See also: TTALIGN, TTPOLAR, TTREFINE, TTSOLVE, SNR2MAXPHASEERROR

%     Version History:
%        Mar.  2, 2010 - initial version (from dtwresid)
%        Mar. 11, 2010 - doc update, column vector output
%        Mar. 12, 2010 - better way to get column vector output, more
%                        input checking
%        Mar. 22, 2010 - account for ttalign fix, ignore diagonal elements
%                        in error calculation (makes error slightly higher)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 22, 2010 at 23:35 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end

% check inputs
if(~isreal(m) || ~isvector(m))
    error('seizmo:ttwstd:badInput',...
        'M must be a real-valued vector!');
elseif(~isreal(lag))
    error('seizmo:ttwstd:badInput',...
        'LAG must be a real-valued array!');
end

% vector to grid
if(isvector(lag))
    lag=ndsquareform(lag,'tomatrix',false);
else
    if(~isequal(lag,permute(-lag,[2 1 3])))
        error('seizmo:ttstderr:badInput',...
            'LAG must be a anti-symmetric matrix of real values!');
    end
end

% number of records
nr=size(lag,1);
if(nr~=numel(m))
    error('seizmo:ttwstd:badInput',...
        'M & LAG are inconsistent in size!');
end

% handle weights
if(nargin==2 || isempty(lagw)); lagw=ones(size(lag)); end
if(isscalar(lagw)); lagw=lagw(ones(size(lag))); end
if(~isreal(lagw))
    error('seizmo:ttwstd:badInput',...
        'LAGW must be a real-valued array!');
end

% vector to grid
if(isvector(lagw))
    lagw=ndsquareform(lagw,'tomatrix');
else
    if(~isequal(lagw,permute(lagw,[2 1 3])))
        error('seizmo:ttstderr:badInput',...
            'LAGW must be a symmetric matrix of real values!');
    end
end

% number of records
nr=size(lagw,1);
if(nr~=numel(m))
    error('seizmo:ttwstd:badInput',...
        'LAG & LAGW are inconsistent in size!');
end

% force m to be column vector
m=m(:);

% ignore diagonal
lag=lag+diag(nan(nr,1));

% get standard error of arrivals
std=sqrt(nanvariance(lag-m(:,ones(nr,1))+m(:,ones(nr,1)).',0,2,lagw));

end
