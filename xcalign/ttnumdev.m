function [numdev]=ttnumdev(lg,dt,std)
%TTNUMDEV    Returns number of standard deviations of lags from solution
%
%    Usage:    numdev=ttnumdev(lg,m,stderr)
%
%    Description: NUMDEV=TTNUMDEV(LG,M,STDERR) returns the number of
%     standard deviations, NUMDEV, each of the lags in LG are from the
%     lags estimated by relative arrivals in M.  This is useful for
%     identifying outliers and for converting lags to a measure that is
%     more effective in weighting schemes.  STDERR is the standard error
%     for each of the relative arrivals in M and must be equal sized with
%     M.
%
%    Notes:
%
%    Examples:
%     Plot the record number vs number of standard deviations of the
%     corresponding lags:
%      nrecs=numel(m);
%      figure; plot(repmat(1:nrecs,nrecs,1),ttnumdev(lag,m,stderr),'x');
%      xlabel('Record Number')
%      ylabel('Number of Standard Deviations')
%
%    See also: TTSTDERR, TTALIGN, CORRELATE

%     Version History:
%        Mar. 22, 2010 - initial version
%        Sep. 13, 2010 - nargchk fix, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 13, 2010 at 01:05 GMT

% todo:

% check nargin
error(nargchk(3,3,nargin));

% check inputs
[nr,nxc,np,vector]=check_xc_info(lg);
check_xc_solutions(nr,dt,std);

% square standard deviations to be variances so we can add them
std=std.^2;

% find number of deviations
if(vector)
    % i=row, j=col in matrix form
    [i,j]=ind2sub([nr nr],find(tril(true(nr),-1)));
    
    % misfit
    numdev=abs((lg-dt(i,1,ones(1,np))+dt(j,1,ones(1,np)))...
        ./sqrt(std(i,1,ones(1,np))+std(j,1,ones(1,np))));
else % matrix
    % misfit
    numdev=abs((lg-dt(:,ones(1,nr),ones(1,np))...
        +permute(dt(:,ones(1,nr),ones(1,np)),[2 1 3]))...
        ./sqrt(std(:,ones(1,nr),ones(1,np))...
        +permute(std(:,ones(1,nr),ones(1,np)),[2 1 3])));
end

end


function [nr,nxc,np,vector]=check_xc_info(lg)
% size up
sz=size(lg);

% require only 3 dimensions
if(numel(sz)==2)
    sz(3)=1;
elseif(numel(sz)>3)
    error('seizmo:ttnumdev:badInput',...
        'LG has too many dimensions!');
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
        error('seizmo:ttnumdev:badInput',...
            'LG is not a properly lengthed vector!');
    end
    
    % check values are in range
    if(~isreal(lg))
        error('seizmo:ttnumdev:badInput',...
            'LG must be a vector of real values!');
    end
else % matrix form
    vector=false;
    
    % get number of records/peaks
    nr=sz(1);
    nxc=nr^2;
    np=sz(3);
    
    % check grids are square
    if(sz(1)~=sz(2))
        error('seizmo:ttnumdev:badInput',...
            'LG are not square matrices!');
    end
    
    % check grids are symmetric & values are in range
    if(~isequal(lg,permute(-lg,[2 1 3])))
        error('seizmo:ttnumdev:badInput',...
            'LG must be a anti-symmetric matrix of real values!');
    end
end

end


function []=check_xc_solutions(nr,dt,std)
% checking that all are correctly sized and valued
if(~isreal(dt) || ~isequal(size(dt),[nr 1]))
    error('seizmo:ttnumdev:badInput',...
        'DT is not a properly sized real-valued column vector!');
elseif(~isreal(std) || ~isequal(size(std),[nr 1]))
    error('seizmo:ttnumdev:badInput',...
        'STD is not a properly sized real-valued column vector!');
end

end
