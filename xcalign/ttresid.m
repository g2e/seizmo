function [resid]=ttresid(lag,dt)
%TTRESID    Returns the time residuals of lags from the solution
%
%    Usage:    resid=ttresid(lag,m)
%
%    Description:
%     RESID=TTRESID(LAG,M) returns the time residual of each lag in LAG
%     from the lags estimated by relative arrivals in M.
%
%    Notes:
%     - LAG is required to be 2D or 3D where the 1st 2 dimensions are NxN
%       or N*(N-1)/2x1 where N is the number of signals being compared.
%
%    Examples:
%     % Plot of record number vs the lag residuals
%     % to the estimated relative arrival times:
%     nrecs=numel(m);
%     figure; plot(repmat(1:nrecs,nrecs,1),ttresid(lag,m),'x');
%     xlabel('Record Number')
%     ylabel('Time Residual of Lags to Estimated Arrivals')
%
%    See also: TTSTDERR, TTALIGN, CORRELATE

%     Version History:
%        Mar. 22, 2010 - initial version
%        Sep. 13, 2010 - nargchk fix
%        Apr.  2, 2012 - minor doc update
%        Jan. 29, 2013 - change lg to lag for readibility, fix warnings
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 29, 2013 at 01:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check inputs
[nr,nxc,np,vector]=check_xc_info(lag);
check_xc_solutions(nr,dt);

% find number of deviations
if(vector)
    % i=row, j=col in matrix form
    [s,m]=ind2sub([nr nr],find(tril(true(nr),-1)));
    
    % misfit
    resid=lag-dt(s,1,ones(1,np))+dt(m,1,ones(1,np));
else % matrix
    % misfit
    resid=lag-dt(:,ones(1,nr),ones(1,np))...
        +permute(dt(:,ones(1,nr),ones(1,np)),[2 1 3]);
end

end


function [nr,nxc,np,vector]=check_xc_info(lag)
% size up
sz=size(lag);

% require only 3 dimensions
if(numel(sz)==2)
    sz(3)=1;
elseif(numel(sz)>3)
    error('seizmo:ttresid:badInput',...
        'LAG has too many dimensions!');
end

% allow either vector or matrix form
if(sz(2)==1)
    vector=true;
    
    % get number of records/peaks
    nr=ceil(sqrt(2*sz(1)));
    nxc=sz(1);
    np=sz(3);
    
    % assure length is ok
    if((nr^2-nr)/2~=sz(1))
        error('seizmo:ttresid:badInput',...
            'LAG is not a properly lengthed vector!');
    end
    
    % check values are in range
    if(~isreal(lag))
        error('seizmo:ttresid:badInput',...
            'LAG must be a vector of real values!');
    end
else % matrix form
    vector=false;
    
    % get number of records/peaks
    nr=sz(1);
    nxc=nr^2;
    np=sz(3);
    
    % check grids are square
    if(sz(1)~=sz(2))
        error('seizmo:ttresid:badInput',...
            'LAG are not square matrices!');
    end
    
    % check grids are symmetric & values are in range
    if(~isequal(lag,permute(-lag,[2 1 3])))
        error('seizmo:ttresid:badInput',...
            'LAG must be a anti-symmetric matrix of real values!');
    end
end

end


function []=check_xc_solutions(nr,dt)
% checking that all are correctly sized and valued
if(~isreal(dt) || ~isequal(size(dt),[nr 1]))
    error('seizmo:ttresid:badInput',...
        'DT is not a properly sized real-valued column vector!');
end

end
