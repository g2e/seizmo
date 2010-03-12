function [cg,lg,pg,NCHANGED]=ttrefine(...
    cg,lg,pg,dt,std,pol,wg,MINSTD,FORCEPOLAR,SKIP)
%TTREFINE   Rearranges cross-correlation peak info
%
%    Usage:    [cg,lg,pg,nc]=ttrefine(cg,lg,pg,dt,std,pol)
%              [cg,lg,pg,nc]=ttrefine(cg,lg,pg,dt,std,pol,wg)
%              [cg,lg,pg,nc]=ttrefine(cg,lg,pg,dt,std,pol,wg,minstd)
%              [cg,lg,pg,nc]=ttrefine(cg,lg,pg,dt,std,pol,wg,minstd,pflag)
%
%    Description: [CG,LG,PG,NC]=TTREFINE(CG,LG,PG,DT,STD,POL) takes the
%     output matrices from a multi-peak run of CORRELATE (CG, LG, PG) and
%     the solutions from TTALIGN (DT), TTSTDERR (STD), & TTPOLAR (POL) to
%     reorder the CG, LG, & PG so the first page of those matrices contain
%     the peaks with the smallest timing difference from DT and a polarity
%     matching that of POL.  If the peak on the first page is within 2
%     standard deviations as given by STD then that peak will not be
%     reordered.  See the third usage form with MINSTD to adjust this.  The
%     output NC is the number of peaks changed.
%
%     [CG,LG,PG,NC]=TTREFINE(CG,LG,PG,DT,STD,POL,WG) weights the misfits
%     using WG.  WG should be empty/scalar (no weighting) or a matrix of
%     size equal to CG/LG/PG.  Weights must be >=0.  Higher weights give
%     lower misfits.  The default is no weighting.
%
%     [CG,LG,PG,NC]=TTREFINE(CG,LG,PG,DT,STD,POL,WG,MINSTD) adjusts the
%     limit at which timing influences the misfits to MINSTD.  By default
%     MINSTD is 2.  Setting MINSTD higher will tend to adjust less peaks
%     based on timing, preferring to adjust them based on their weights
%     given in WG.  The default does well when the standard deviations in
%     STD are from TTSTDERR (ie weighted).
%
%     [CG,LG,PG,NC]=TTREFINE(CG,LG,PG,DT,STD,POL,WG,MINSTD,PFLAG) sets if
%     peaks with the wrong polarity are allowed or not.  By default PFLAG
%     is TRUE and this indicates the peaks must match the polarities in
%     POL.  Setting PFLAG to FALSE will allow the polarity to be refined
%     too.
%
%    Notes:
%     - Why reorder peak info?  Re-solving the reordered peak info gives a
%       more consistent/precise and generally more accurate alignment of
%       correlated signals when there is noise in the data.  Often in this
%       case the highest peak in a correlogram does not correspond to the
%       correct shift between 2 signals as the noise may raise another peak
%       while lowering the correct peak.  But noise being chaotic, this is
%       not always the case and generally the highest peaks scatter in time
%       about the true solution (which means that TTSOLVE gets close the
%       first time).  TTREFINE chooses peaks that best match the solution
%       from TTSOLVE (& TTPOLAR).  Essentially the accuracy of the
%       reordered solution heavily depends on the accuracy of DT & POL.
%
%    Examples:
%     Correlate, solve, plot, refine, solve again, plot:
%      plot0(data);
%      xc=correlate(data,'npeaks',3,'spacing',10);
%      dt=ttalign(xc.lg(:,:,1));
%      data1=timeshift(data,dt);
%      pol=ttpolar(xc.pg(:,:,1));
%      data1=multiply(data1,pol);
%      std=ttstderr(dt,xc.lg(:,:,1));
%      plot0(data1);
%      [xc2.cg,xc2.lg,xc2.pg,nc]=...
%          ttrefine(xc.cg,xc.lg,xc.pg,dt,std,pol,[],0.1,0);
%      dt=ttalign(xc2.lg(:,:,1));
%      data1=timeshift(data,dt);
%      pol=ttpolar(xc2.pg(:,:,1));
%      data1=multiply(data1,pol);
%      std=ttstderr(dt,xc2.lg(:,:,1));
%      plot0(data1);
%
%    See also: TTSOLVE, TTALIGN, TTSTDERR, TTPOLAR

%     Version History:
%        Mar. 11, 2010 - initial version (derived from dtwrefine)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 11, 2010 at 17:45 GMT

% todo:

% check nargin
msg=nargchk(7,10,nargin);
if(~isempty(msg)); error(msg); end

% defaults
if(nargin<7 || isempty(wg)); wg=1; end
if(nargin<8 || isempty(MINSTD)); MINSTD=2; end
if(nargin<9 || isempty(FORCEPOLAR)); FORCEPOLAR=true; end
if(nargin<10 || isempty(SKIP)); SKIP=false; end

% check inputs
% - cg, lg, pg, wg should all be the same
%   - wg can be scalar (expanded)
% - dt, std, pol should all be the same
[nr,nxc,np,vector]=check_xc_info(cg,lg,pg);
check_xc_solutions(nr,dt,std,pol);
wg=check_xc_weights(cg,wg);
if(~isscalar(MINSTD) || ~isreal(MINSTD))
    error('seizmo:ttrefine:badInput',...
        'MINSTD must be a real-valued scalar!');
elseif(~isscalar(FORCEPOLAR) ...
        || (~islogical(FORCEPOLAR) && ~isreal(FORCEPOLAR)))
    error('seizmo:ttrefine:badInput',...
        'FORCEPOLAR must be a logical scalar!');
elseif(~isscalar(SKIP) || (~islogical(SKIP) && ~isreal(SKIP)))
    error('seizmo:ttrefine:badInput',...
        'SKIP must be a logical scalar!');
end

% find misfit
% - essentially the number of standard deviations the measured lag is from
%    the expected lag
% - lags under MINSTD standard deviations are all assigned MINSTD standard
%    deviations (limits the effect of lags on the misfit)
if(vector)
    % i=row, j=col in matrix form
    [i,j]=ind2sub([nr nr],find(tril(true(nr),-1)));
    
    % misfit
    M=1./wg.*max(MINSTD,...
        sqrt(((lg+dt(i,1,ones(1,np))-dt(j,1,ones(1,np)))...
        ./(std(i,1,ones(1,np))+std(j,1,ones(1,np)))).^2));
    
    % to get number adjusted
    factor=1;
else
    % misfit
    M=1./wg.*max(MINSTD,...
        sqrt(((lg+dt(:,ones(1,nr),ones(1,np))...
        -permute(dt(:,ones(1,nr),ones(1,np)),[2 1 3]))...
        ./(std(:,ones(1,nr),ones(1,np))...
        +permute(std(:,ones(1,nr),ones(1,np)),[2 1 3]))).^2));
    
    % to get number adjusted
    factor=2;
end

% IGNORE ANY PEAKS NOT OF THE CORRECT POLARITY
if(FORCEPOLAR)
    pol=pol*pol';
    pol=pol(:,:,ones(npeaks,1));
    M(pg~=pol)=nan;
end

% FIND BEST PEAKS NOT ON PAGE 1 (COMPARES DOWN PAGES = 3RD DIM)
[v,mi]=min(M,[],3);
li=find(mi>1);
NCHANGED=size(li,1)/factor;

% QUICK EXIT
if(~SKIP && NCHANGED==0); return; end

% LINEAR INDICES OF BEST PEAKS
best=li+(mi(li)-1)*nxc;

% SPECIAL EXIT (FOR DEBUGGING & EDGE CASES)
if(SKIP)
    % output li, best, M, NCHANGED
    cg=li;
    lg=best;
    pg=M;
    return;
end

% MOVING BEST PEAKS NOT ON PAGE 1 TO TEMPORARY VECTOR
ctmp=cg(best);
ltmp=lg(best);
ptmp=pg(best);

% MOVING PAGE 1 PEAKS THAT ARE NOT THE BEST
cg(best)=cg(li);
lg(best)=lg(li);
pg(best)=pg(li);

% MOVING BEST PEAKS NOT ON PAGE 1 TO PAGE 1
cg(li)=ctmp;
lg(li)=ltmp;
pg(li)=ptmp;

end

function [nr,nxc,np,vector]=check_xc_info(cg,lg,pg)
% size up
sz=size(cg);

% all must be the same size
if(~isequal(sz,size(lg),size(pg)))
    error('seizmo:ttrefine:badInput',...
        'CG/LG/PG must be equal sized arrays!');
end

% require only 3 dimensions
if(numel(sz)==2)
    sz(3)=1;
elseif(numel(sz)>3)
    error('seizmo:ttrefine:badInput',...
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
        error('seizmo:ttrefine:badInput',...
            'CG/LG/PG are not properly lengthed vectors!');
    end
    
    % check values are in range
    if(~isreal(cg) || any(abs(cg(:)>1)))
        error('seizmo:ttrefine:badInput',...
            'CG must be a vector of real values within -1 & 1!');
    elseif(~isreal(lg))
        error('seizmo:ttrefine:badInput',...
            'LG must be a vector of real values!');
    elseif(any(abs(pg(:))~=1))
        error('seizmo:ttrefine:badInput',...
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
        error('seizmo:ttrefine:badInput',...
            'CG/LG/PG are not square matrices!');
    end
    
    % check grids are symmetric & values are in range
    if(~isequal(cg,permute(cg,[2 1 3])) || any(abs(cg(:))>1))
        error('seizmo:ttrefine:badInput',...
            'CG must be a symmetric matrix of real values within -1 & 1!');
    elseif(~isequal(lg,-lg'))
        error('seizmo:ttrefine:badInput',...
            'LG must be a anti-symmetric matrix of real values!');
    elseif(~isequal(pg,pg') || any(abs(pg(:))~=1))
        error('seizmo:ttrefine:badInput',...
            'PG must be a symmetric matrix of 1s & -1s!');
    end
end

end

function [wg]=check_xc_weights(cg,wg)
% check values
if(~isreal(wg) || any(wg(:)<0))
    error('seizmo:ttrefine:badInput',...
        'WG must be an array of positive reals!');
end
% check size
if(isscalar(wg))
    % expand scalar
    wg=wg(ones(size(cg)));
elseif(~isequal(size(cg),size(wg)))
    error('seizmo:ttrefine:badInput',...
        'Non-scalar WG must match size of CG/LG/PG!');
end

end

function []=check_xc_solutions(nr,dt,std,pol)
% checking that all are correctly sized and valued
if(~isreal(dt) || ~isequal(size(dt),[nr 1]))
    error('seizmo:ttrefine:badInput',...
        'DT is not a properly sized real-valued column vector!');
elseif(~isreal(std) || ~isequal(size(std),[nr 1]))
    error('seizmo:ttrefine:badInput',...
        'STD is not a properly sized real-valued column vector!');
elseif(~isreal(pol) || any(abs(pol)~=1) || ~isequal(size(pol),[nr 1]))
    error('seizmo:ttrefine:badInput',...
        'POL is not a properly sized column vector of 1s & -1s!');
end

end
