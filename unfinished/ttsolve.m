function [dt,std,pol,zmean,zstd,nc]=ttsolve(xc,mtri,mpri,snr,dt,std,pol)
%TTSOLVE    Solves relative arrival times & polarities
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also: TTALIGN, TTREFINE, TTPOLAR, TTSTDERR

%     Version History:
%        Mar. 12, 2010 - initial version (derived from groupeval)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 12, 2010 at 14:05 GMT

% todo:
% - everything but xc needs to be an option given by
%    ttsolve(xc,'option1',value1,'option2',value2,...)

% check nargin
msg=nargchk(1,7,nargin);
if(~isempty(msg)); error(msg); end

% verbosity
verbose=seizmoverbose;

% additional arguments not allowed as input
minstd=0.1;    % lowest lag standard error in misfit calc
nonconv=1;     % non-convergence method
maxcorr=0.999; % max correlation value for weights/stats
wp=2;          % weight power

% check xc struct
if(~isstruct(xc) || any(~ismember({'cg' 'lg' 'pg'},fieldnames(xc))))
    error('seizmo:ttsolve:badInput',...
        'XC must be a struct with fields .cg .lg .pg !');
end
[nr,nxc,np,vector]=check_xc_info(xc.cg,xc.lg,xc.pg);

% get defaults
if(nargin<2 || isempty(mtri)); mtri=inf; end
if(nargin<3 || isempty(mpri)); mpri=inf; end
if(nargin<4 || isempty(snr)); snr=ones(nr,1); end
if(nargin<5); dt=[]; end
if(nargin<6); std=[]; end
if(nargin<7); pol=[]; end

% check optional arguments
if(~isscalar(mtri) || mtri~=round(mtri'))
    error('seizmo:ttsolve:badInput',...
        'MTRI must be a scalar integer!');
elseif(~isscalar(mpri) || mpri~=round(mpri'))
    error('seizmo:ttsolve:badInput',...
        'MPRI must be a scalar integer!');
end
check_xc_estimates(nr,snr,dt,std,pol);

% force correlations to be under some maximum
% - this stabilizes the inversion for cases
%   when there are identical records involved
xc.cg(xc.cg>maxcorr)=maxcorr;

% convert correlation values to z-statistic
% - more appropriate for weighting
% - more appropriate for error estimation
xc.zg=fisher(xc.cg);

% now get weights (z-statistic * root of product of rescaled snr)
snr=rescale_snr(snr);
xc.wg=sqrt(snr*snr');
xc.wg=xc.wg(:,:,ones(np,1));
if(vector)
    xc.wg=permute(ndsquareform(xc.wg),[2 1 3]).*(xc.zg);
else % matrix
    xc.wg=xc.wg.*xc.zg;
end
xc.wg=xc.wg.^wp;

% initial alignment/error/polarity estimates
if(isempty(dt))
    if(verbose); disp('Inverting for Initial Relative Arrivals'); end
    dt=ttalign(xc.lg(:,:,1),xc.wg(:,:,1));
else
    if(verbose); disp('Initial Relative Arrivals Given as Input'); end
end
if(isempty(std))
    if(verbose); disp('Determining Standard Errors'); end
    std=ttstderr(dt,xc.lg(:,:,1),xc.wg(:,:,1));
else
    if(verbose); disp('Standard Errors Given as Input'); end
end
if(isempty(pol))
    if(verbose); disp('Inverting for Initial Relative Polarities'); end
    pol=ttpolar(xc.pg(:,:,1));
else
    if(verbose); disp('Initial Relative Polarities Given as Input'); end
end

% refine alignment estimate
% - selects lag with minimum weighted misfit to initial alignment among the
%   peaks with the same polarity as the initial estimate
% - reinvert using with these more consistent peaks
% - repeat until max number of iterations or converged
iter=1; nc=nan(0,1);
while(iter<=mtri)
    % are we refining polarities?
    if(iter<=mpri)
        % yes
        if(verbose)
            disp('Refining Relative Arrivals & Errors & Polarities');
        end
        pflag=false;
    else
        % no
        if(verbose); disp('Refining Relative Arrivals & Errors'); end
        pflag=true;
    end
    
    % refine
    [li,best]=ttrefine(...
        xc.cg,xc.lg,xc.pg,dt,std,pol,xc.wg,minstd,pflag,true);
    [xc.cg,xc.lg,xc.pg,nc(iter,1)]=ttrefine(...
        xc.cg,xc.lg,xc.pg,dt,std,pol,xc.wg,minstd,pflag);
    
    % check for change
    if(nc(iter,1)>0)
        if(verbose)
            disp(['Reordered ' num2str(nc(iter,1)) ' Peaks']);
        end
    else
        % drop this iteration
        nc(iter)=[];
        break;
    end
    
    % reorder z-stats & weights
    tmp=xc.zg(best);
    xc.zg(best)=xc.zg(li);
    xc.zg(li)=tmp;
    tmp=xc.wg(best);
    xc.wg(best)=xc.wg(li);
    xc.wg(li)=tmp;
    
    % reinvert
    if(verbose); disp('Reinverting Reordered Peak Info'); end
    dt=ttalign(xc.lg(:,:,1),xc.wg(:,:,1));
    std=ttstderr(dt,xc.lg(:,:,1),xc.wg(:,:,1));
    if(pflag)
        % check that polarity solution did not change
        % - this is just to make sure I coded this right
        if(~isequal(pol,ttpolar(xc.pg(:,:,1))))
            error('seizmo:ttsolve:brokenPolarities',...
                'Polarity solution changed when it should not have!');
        end
    else
        % new polarity solution
        pol=ttpolar(xc.pg(:,:,1));
    end
    
    % detect non-convergence (after at least 3 iterations)
    % - number peaks reordered stays constant over 3 iterations
    %   ie 2, 2, 2
    % - number changed alternates between 2 numbers twice
    %   ie 5, 4, 5, 4
    if((iter>=3 && isscalar(unique(nc(iter-2:iter,1)))) ...
            || (iter>=4 && isequal(nc(iter-3:iter-2,1),nc(iter-1:iter,1))))
        % non-convergence!
        if(verbose); disp('Non-Convergence Detected'); end
        
        % handle non-convergence
        switch nonconv
            case 1
                % just break
                if(verbose); disp('Immediately Breaking Refinement'); end
                break;
            case 2
                % downweight troubled peaks, re-solve, break
                if(verbose); disp('Downweighting Troubled Peaks'); end
                xc.wg(li)=0;
                if(verbose); disp('Reinverting 1 More Time'); end
                dt=ttalign(xc.lg(:,:,1),xc.wg(:,:,1));
                std=ttstderr(dt,xc.lg(:,:,1),xc.wg(:,:,1));
                if(verbose); disp('Breaking Refinement'); end
                break;
            case 3
                % downweight troubled peaks, re-solve, continue on
                if(verbose); disp('Downweighting Troubled Peaks'); end
                xc.wg(li)=0;
                if(verbose); disp('Reinverting Downweighted Peaks'); end
                dt=ttalign(xc.lg(:,:,1),xc.wg(:,:,1));
                std=ttstderr(dt,xc.lg(:,:,1),xc.wg(:,:,1));
            otherwise
                error('seizmo:ttsolve:improperMethod',...
                    'NONCONV must be 1, 2, or 3!');
        end
    end
    
    % increment counter
    iter=iter+1;
end

% get zmean, zstd for reordered matrices
if(vector)
    zmean=nanmean(ndsquareform(xc.zg(:,:,1)),2);
    zstd=sqrt(nanvariance(ndsquareform(xc.zg(:,:,1)),...
        0,2,ndsquareform(xc.wg(:,:,1))));
else
    zmean=nanmean(xc.zg(:,:,1),2);
    zstd=sqrt(nanvariance(xc.zg(:,:,1),0,2,xc.wg(:,:,1)));
end

end


function [nr,nxc,np,vector]=check_xc_info(cg,lg,pg)
% size up
sz=size(cg);

% all must be the same size
if(~isequal(sz,size(lg),size(pg)))
    error('seizmo:ttsolve:badInput',...
        'CG/LG/PG must be equal sized arrays!');
end

% require only 3 dimensions
if(numel(sz)==2)
    sz(3)=1;
elseif(numel(sz)>3)
    error('seizmo:ttsolve:badInput',...
        'CG/LG/PG have too many dimensions!');
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
        error('seizmo:ttsolve:badInput',...
            'CG/LG/PG are not properly lengthed vectors!');
    end
    
    % check values are in range
    if(~isreal(cg) || any(abs(cg(:)>1)))
        error('seizmo:ttsolve:badInput',...
            'CG must be a vector of real values within -1 & 1!');
    elseif(~isreal(lg))
        error('seizmo:ttsolve:badInput',...
            'LG must be a vector of real values!');
    elseif(any(abs(pg(:))~=1))
        error('seizmo:ttsolve:badInput',...
            'PG must be a vector of 1s & -1s!');
    end
else % matrix form
    vector=false;
    
    % get number of records/peaks
    nr=sz(1);
    nxc=nr^2;
    np=sz(3);
    
    % check grids are square
    if(sz(1)~=sz(2))
        error('seizmo:ttsolve:badInput',...
            'CG/LG/PG are not square matrices!');
    end
    
    
    % check grids are symmetric & values are in range
    if(~isequal(cg,permute(cg,[2 1 3])) || any(abs(cg(:))>1))
        error('seizmo:ttsolve:badInput',...
            'CG must be a symmetric matrix of real values within -1 & 1!');
    elseif(~isequal(lg,permute(-lg,[2 1 3])))
        error('seizmo:ttsolve:badInput',...
            'LG must be a anti-symmetric matrix of real values!');
    elseif(~isequal(pg,permute(pg,[2 1 3])) || any(abs(pg(:))~=1))
        error('seizmo:ttsolve:badInput',...
            'PG must be a symmetric matrix of 1s & -1s!');
    end
end

end


function []=check_xc_estimates(nr,snr,dt,std,pol)
% checking that all are correctly sized and valued
if(~isreal(snr) || ~isequal(size(snr),[nr 1]) || any(snr<0))
    error('seizmo:ttsolve:badInput',...
        'SNR is not a properly sized column vector of positive reals!');
elseif(~isempty(dt) && (~isreal(dt) || ~isequal(size(dt),[nr 1])))
    error('seizmo:ttsolve:badInput',...
        'DT is not a properly sized real-valued column vector!');
elseif(~isempty(std) && (~isreal(std) || ~isequal(size(std),[nr 1])))
    error('seizmo:ttsolve:badInput',...
        'STD is not a properly sized real-valued column vector!');
elseif(~isempty(pol) && (~isreal(pol) || any(abs(pol)~=1) ...
        || ~isequal(size(pol),[nr 1])))
    error('seizmo:ttsolve:badInput',...
        'POL is not a properly sized column vector of 1s & -1s!');
end

end

