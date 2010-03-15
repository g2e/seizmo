function [dt,std,pol,zmean,zstd,nc]=ttsolve(xc,varargin)
%TTSOLVE    Solves relative arrival times & polarities
%
%    Usage:    [arr,err,pol,zmean,zstd,nc]=ttsolve(xc)
%
%            weighting options:
%              [...]=ttsolve(...,'snr',snr,...)
%              [...]=ttsolve(...,'wgtpow',value,...)
%              [...]=ttsolve(...,'inrange',[low high],...)
%              [...]=ttsolve(...,'outrange',[low high],...)
%              [...]=ttsolve(...,'minlag',stddev,...)
%              [...]=ttsolve(...,'maxcor',rval,...)
%
%            starting solution options:
%              [...]=ttsolve(...,'estarr',relarr,...)
%              [...]=ttsolve(...,'esterr',stderr,...)
%              [...]=ttsolve(...,'estpol',relpol,...)
%
%            style/depth options:
%              [...]=ttsolve(...,'mtri',niter,...)
%              [...]=ttsolve(...,'mpri',niter,...)
%              [...]=ttsolve(...,'noncnv',method,...)
%
%    Description: [ARR,ERR,POL,ZMEAN,ZSTD,NC]=TTSOLVE(XC) takes the
%     cross-correlation struct from a CORRELATE call that was passed with
%     option 'NPEAKS' set to 1 or above and solves for the relative arrival
%     times and polarities between the correlated signals.  The inversion
%     utilizes an iterative, weighted least-squares approach that makes use
%     of multiple peaks in a correlogram.  This approach is particularly
%     suited for aligning noisy, narrow-band signals.  There are a
%     multitude of options that allow for adjusting the weighting in the
%     inversion, the starting conditions of the inversion, and the style
%     and depth of the inversion.
%
%     %%%%%%%%%%%%%%%%%%%%%%%%
%     SOLVER WEIGHTING OPTIONS
%     %%%%%%%%%%%%%%%%%%%%%%%%
%
%     [...]=TTSOLVE(...,'SNR',SNR,...) sets the signal to noise ratio to
%     something other than one (the default).  This is for weighting when
%     inverting for relative arrival times and errors.  Please note that
%     the values given in SNR will be rescaled utilizing RESCALE_SNR with
%     its default ranges.  Use options 'INRANGE' & 'OUTRANGE' to adjust the
%     rescaling.  SNR must be a column vector of positive values.
%
%     [...]=TTSOLVE(...,'WGTPOW',POWER,...) applies raises the weights used
%     in the inversion by POWER.  The default is 1 (no change in weight).
%     Setting POWER to 0 will cause the inversion to be completely driven
%     by lag misfit and polarities (ie SNR and cross-correlation
%     coefficient will have no influence).  Setting POWER to 2+ will cause
%     the SNR and cross-correlation coefficients to have a stronger
%     influence in the inversion.
%
%     [...]=TTSOLVE(...,'INRANGE',[LOW HIGH],...) sets the in rescale range
%     to [LOW HIGH] for rescaling SNR values to weights (see option 'SNR'
%     to set the SNR values).  See function RESCALE_SNR for more details.
%
%     [...]=TTSOLVE(...,'OUTRANGE',[LOW HIGH],...) sets the out rescale
%     range to [LOW HIGH] for rescaling SNR values to weights (see option
%     'SNR' to set the SNR values).  See function RESCALE_SNR for more
%     details.
%
%     [...]=TTSOLVE(...,'MINLAG',STDDEV,...) resets any lag's standard
%     deviation that is lower than STDDEV to STDDEV.  This allows for
%     limiting the influence of lags near the estimated arrival time's lag
%     prediction.  The default is 0.1 and seems to do pretty well.  This
%     means that a lag off by 1 standard deviation has 10 times the
%     unweighted misfit as that of a lag right at the predicted time.
%     Setting this value higher usually slows down the convergence to a
%     solution but favors highly correlated signals.  Setting this number
%     too low could cause the solution to align on a poorly correlated
%     signal (although this may be useful to align on a particular signal
%     even when there are other strong signals degrading the correlation).
%
%     [...]=TTSOLVE(...,'MAXCOR',RVAL,...) resets any cross-correlation
%     coefficient above RVAL to RVAL.  This limits the influence of highly
%     correlated signals in the inversion and in the ZMEAN/ZSTD output.
%     The default is 0.999 and mainly eliminates the problem of identical
%     signals that destabilize the inversion.
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     SOLVER STARTING SOLUTION OPTIONS
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     [...]=TTSOLVE(...,'ESTARR',RELARR,...) sets the initial estimate of
%     the relative arrival times to RELARR.  The inversion will guide the
%     refinement of relative arrival times to be near these estimates.  The
%     default is no estimate, which causes TTSOLVE to get an estimate using
%     the lags given in XC.lg(:,:,1).  This usually corresponds to the
%     strongest peaks in the correlograms which generally give a close
%     estimate of the actual relative arrival times.
%
%     [...]=TTSOLVE(...,'ESTERR',STDERR,...) sets the initial estimate of
%     error in the relative arrival times to STDERR.  This can be useful
%     for driving the first iteration of the refinement.  The default is no
%     estimate, which cause TTSOLVE to estimate the error using the
%     consistency of the lags in XC to the estimated relative arrival times
%     (see option 'ESTARR' above).
%
%     [...]=TTSOLVE(...,'ESTPOL',RELPOL,...) sets the initial estimate of
%     the relative polarity to RELPOL.  This will guide the refinement of
%     relative arrival times (see option 'MPRI' for more information on
%     solving for relative polarities).  The default is no estimate, which
%     causes TTSOLVE to get an estimate using the relative polarities given
%     in XC.
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%     SOLVER STYLE/DEPTH OPTIONS
%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     [...]=TTSOLVE(...,'MTRI',NITER,...) sets the maximum number of
%     relative arrival time refinement iterations to NITER.  If NITER
%     iterations are performed, the solver will exit with the most refined
%     (the last iteration's) solution. The default NITER is infinity.
%     Setting NITER to 0 will perform no refinement while setting it to 10
%     will allow up to 10 iterations.  TTSOLVE may exit before NITER
%     iterations if no more peaks are to be reordered or non-convergence is
%     detected (see option 'NONCNV').
%
%     [...]=TTSOLVE(...,'MPRI',NITER,...) sets the maximum number of
%     relative polarity refinement iterations to NITER.  The default is 0
%     which means the polarity is not refined from the initial estimate
%     (see option 'ESTPOL').  In fact, the solution is forced to be
%     consistent with the polarity estimate after NITER iterations.  This
%     allows for strong consistency between the polarity and arrival
%     solutions.  It is extremely rare for the polarity to be refined so
%     this is best left alone.  It is also a strong driving force for
%     solution convergence.
%
%     [...]=TTSOLVE(...,'NONCNV',METHOD,...) alters how non-convergence is
%     handled.  The default METHOD is 'BREAK' which just exits with the
%     last solution when the non-convergence was detected.  Other options
%     are 'IGNORE' & 'REWEIGHT'.  The 'IGNORE' method do nothing when a
%     non-convergence is detected (not setting option 'MTRI' to a finite
%     number when using the 'IGNORE' method may lead to an infinite loop).
%     The 'REWEIGHT' method will downweight the troubled peaks (those that
%     keep changing) temporarily to get a solution not influenced by those
%     peaks and then continue the refinement with that solution.  This may
%     lead to a better solution in difficult cases.
%
%    Notes:
%
%    Examples:
%     Get the SNR of some signals roughly aligned near 0:
%      snr=quicksnr(data,[-100 -10],[-10 100]);
%     Correlate the records:
%      xc=correlate(data,'npeaks',3,'spacing',10);
%     Solve using signal-to-noise info in the weights:
%      [arr,err,pol,zmean,zstd,nc]=ttsolve(xc,'snr',snr);
%     Now plot the alignment:
%      recordsection(multiply(timeshift(normalize(data),dt),pol));
%
%    See also: TTALIGN, TTREFINE, TTPOLAR, TTSTDERR

%     Version History:
%        Mar. 12, 2010 - initial version (derived from groupeval)
%        Mar. 14, 2010 - options worked out, documentation added
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 14, 2010 at 04:05 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% verbosity
verbose=seizmoverbose;

% check xc struct
if(~isstruct(xc) || any(~ismember({'cg' 'lg' 'pg'},fieldnames(xc))))
    error('seizmo:ttsolve:badInput',...
        'XC must be a struct with fields .cg .lg .pg !');
end
[nr,nxc,np,vector]=check_xc_info(xc.cg,xc.lg,xc.pg);

% get options
opt=parse_ttsolve_args(nr,varargin{:});

% force correlations to be under some maximum
% - this stabilizes the inversion for cases
%   when there are identical records involved
xc.cg(xc.cg>opt.MAXCOR)=opt.MAXCOR;

% convert correlation values to z-statistic
% - more appropriate for weighting
% - more appropriate for error estimation
xc.zg=fisher(xc.cg);

% now get weights (z-statistic * root of product of rescaled snr)
opt.SNR=rescale_snr(opt.SNR,opt.INRANGE,opt.OUTRANGE);
xc.wg=sqrt(opt.SNR*opt.SNR');
xc.wg=xc.wg(:,:,ones(np,1));
if(vector)
    xc.wg=permute(ndsquareform(xc.wg),[2 1 3]).*(xc.zg);
else % matrix
    xc.wg=xc.wg.*xc.zg;
end
xc.wg=xc.wg.^opt.WGTPOW;

% initial alignment/error/polarity estimates
if(isempty(opt.ESTARR))
    if(verbose); disp('Inverting for Initial Relative Arrivals'); end
    dt=ttalign(xc.lg(:,:,1),xc.wg(:,:,1));
else
    if(verbose); disp('Initial Relative Arrival Estimate Given'); end
    dt=opt.ESTARR;
end
if(isempty(opt.ESTERR))
    if(verbose); disp('Determining Standard Errors'); end
    std=ttstderr(dt,xc.lg(:,:,1),xc.wg(:,:,1));
else
    if(verbose); disp('Initial Error Estimate Given'); end
    std=opt.ESTERR;
end
if(isempty(opt.ESTPOL))
    if(verbose); disp('Inverting for Initial Relative Polarities'); end
    pol=ttpolar(xc.pg(:,:,1));
else
    if(verbose); disp('Initial Relative Polarity Estimate Given'); end
    pol=opt.ESTPOL;
end

% refine alignment estimate
% - selects lag with minimum weighted misfit to initial alignment among the
%   peaks with the same polarity as the initial estimate
% - reinvert using with these more consistent peaks
% - repeat until max number of iterations or converged
iter=1; nc=nan(0,1);
while(iter<=opt.MTRI)
    % are we refining polarities?
    if(iter<=opt.MPRI)
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
        xc.cg,xc.lg,xc.pg,dt,std,pol,xc.wg,opt.MINLAG,pflag,true);
    [xc.cg,xc.lg,xc.pg,nc(iter,1)]=ttrefine(...
        xc.cg,xc.lg,xc.pg,dt,std,pol,xc.wg,opt.MINLAG,pflag);
    
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
        % handle non-convergence
        switch lower(opt.NONCNV)
            case 'break'
                % just break
                if(verbose); disp('Non-Convergence Detected'); end
                if(verbose); disp('Immediately Breaking Refinement'); end
                break;
            case 'reweight'
                % downweight troubled peaks, re-solve, continue on
                if(verbose); disp('Non-Convergence Detected'); end
                if(verbose); disp('Downweighting Troubled Peaks'); end
                tmp=xc.wg(li);
                xc.wg(li)=0;
                if(verbose); disp('Reinverting Downweighted Version'); end
                dt=ttalign(xc.lg(:,:,1),xc.wg(:,:,1));
                std=ttstderr(dt,xc.lg(:,:,1),xc.wg(:,:,1));
                xc.wg(li)=tmp;
            case 'ignore'
                % free as can be
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


function [opt]=parse_ttsolve_args(nr,varargin)

% check nargin
if(mod(nargin-1,2))
    error('seizmo:ttsolve:badNumInputs',...
        'Unpaired OPTION/VALUE!');
end

% valid string options
valid.NONCNV={'break' 'reweight' 'ignore'};

% option defaults
opt.MTRI=inf;       % max relative arrival refinement iterations
opt.MPRI=0;         % max relative polarity refinement iterations
opt.MINLAG=0.100;   % min standard error of lags (keeps misfit sane)
opt.MAXCOR=0.999;   % max correlation value (keeps misfit/stats sane)
opt.NONCNV='break'; % break/reweight/ignore
opt.WGTPOW=1;       % power to apply to misfit weights
opt.ESTARR=[];      % initial relative arrival estimate
opt.ESTERR=[];      % initial relative arrival error estimate
opt.ESTPOL=[];      % initial relative polarity estimate
opt.SNR=ones(nr,1); % snr (for misfit weights)
opt.INRANGE=[];     % snr range over which rescaling varies
opt.OUTRANGE=[];    % range of rescaled snr values

% option must be specified by a string
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:ttsolve:badInput',...
        'OPTION must be specified by a string!');
end

% loop over option/value pairs
for i=1:2:nargin-1
    % skip if empty (even unknown)
    if(isempty(varargin{i+1}))
        continue;
    end
    
    % check and add
    switch lower(varargin{i})
        case {'mt' 'mtri' 'maxtimeiter'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}~=fix(varargin{i+1}))
                error('seizmo:ttsolve:badInput',...
                    'MTRI option must be a scalar integer!');
            end
            opt.MTRI=varargin{i+1};
        case {'mp' 'mpri' 'maxpoliter'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}~=fix(varargin{i+1}))
                error('seizmo:ttsolve:badInput',...
                    'MPRI option must be a scalar integer!');
            end
            opt.MPRI=varargin{i+1};
        case {'ml' 'minlag' 'minlagerr'}
            if(~isscalar(varargin{i+1}) || varargin{i+1}<0)
                error('seizmo:ttsolve:badInput',...
                    'MINLAG option must be a scalar >=0!');
            end
            opt.MINLAG=varargin{i+1};
        case {'mc' 'maxcor'}
            if(~isscalar(varargin{i+1}) ...
                    || varargin{i+1}<0 || varargin{i+1}>1)
                error('seizmo:ttsolve:badInput',...
                    'MAXCOR option must be a scalar in range 0-1!');
            end
            opt.MAXCOR=varargin{i+1};
        case {'nc' 'noncon' 'nonconv' 'noncnv'}
            if(~ischar(varargin{i+1}) || ...
                    ~ismember(varargin{i+1},valid.NONCNV))
                error('seizmo:ttsolve:badInput',...
                    ['NONCNV option must be one of the following:\n' ...
                    sprintf('''%s'' ',valid.NONCNV)]);
            end
            opt.NONCNV=varargin{i+1};
        case {'wp' 'wgtpow'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('seizmo:ttsolve:badInput',...
                    'WGTPOW option must be a real-valued scalar!');
            end
            opt.WGTPOW=varargin{i+1};
        case {'ra' 'relarr' 'dt' 'estarr' 'arr'}
            if(~isreal(varargin{i+1}) ...
                    || ~isequal(size(varargin{i+1}),[nr 1]))
                error('seizmo:ttsolve:badInput',...
                    ['ESTARR option must be a properly ' ...
                    'sized real-valued column vector!']);
            end
            opt.ESTARR=varargin{i+1};
        case {'se' 'std' 'stderr' 'esterr' 'err'}
            if(~isreal(varargin{i+1}) || any(varargin{i+1}(:)<0) ...
                    || ~isequal(size(varargin{i+1}),[nr 1]))
                error('seizmo:ttsolve:badInput',...
                    ['ESTERR option must be a properly ' ...
                    'sized column vector of positive reals!']);
            end
            opt.ESTERR=varargin{i+1};
        case {'rp' 'pol' 'relpol' 'estpol' 'pol'}
            if(~isreal(varargin{i+1}) || any(abs(varargin{i+1}(:))~=1) ...
                    || ~isequal(size(varargin{i+1}),[nr 1]))
                error('seizmo:ttsolve:badInput',...
                    ['ESTPOL option must be a properly ' ...
                    'sized column vector of 1s & -1s!']);
            end
            opt.ESTPOL=varargin{i+1};
        case {'snr'}
            if(~isreal(varargin{i+1}) || any(varargin{i+1}(:)<0) ...
                    || ~isequal(size(varargin{i+1}),[nr 1]))
                error('seizmo:ttsolve:badInput',...
                    ['SNR option must be a properly ' ...
                    'sized column vector of positive reals!']);
            end
            opt.SNR=varargin{i+1};
        case {'in' 'inrange' 'rescale_inrange'}
            if(~isequal(size(varargin{i+1}),[1 2]) ...
                    || ~isreal(varargin{i+1}) || any(varargin{i+1}<0))
                error('seizmo:ttsolve:badInput',...
                    'RESCALE_INRANGE must be a 1x2 positive real array!');
            end
            opt.INRANGE=varargin{i+1};
        case {'out' 'outrange' 'rescale_outrange'}
            if(~isequal(size(varargin{i+1}),[1 2]) ...
                    || ~isreal(varargin{i+1}) || any(varargin{i+1}<0))
                error('seizmo:ttsolve:badInput',...
                    'RESCALE_OUTRANGE must be a 1x2 positive real array!');
            end
            opt.OUTRANGE=varargin{i+1};
        otherwise
            error('seizmo:ttsolve:badInput',...
                'Unknown Option: %s',varargin{i});
    end
end

end
