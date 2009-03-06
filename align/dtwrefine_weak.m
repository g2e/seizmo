function [cv,lv,pv,NCHANGED,li,M]=dtwrefine_weak(dt,cv,lv,pv,xcpow,xcw,frac,SKIP)
%DTWREFINE_WEAK    Rearranges peak info based on misfit to inverted arrival times
%
% NOTE: Refines lag sets that have misfits exceeding some fraction of the worst misfit.
%
% INPUTS REQUIRED:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dt    - calculated relative arrival times
% cv/cg - vector/grid of xc correlation values (used for weighting)
% lv/lg - vector/grid of corresponding xc lag times
% pv/pg - vector/grid of corresponding xc polarities
% xcpow - power(s) to apply to xc correlation values for weighting (0=no xc value weight)
% xcw   - weighting coefficient(s) to apply in addition to xc correlation values
% frac  - refine only misfits above frac times the max misfit
% SKIP  - logical flag -> TRUE = don't rearrange
%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cv/cg, lv/lg, pv/pg - Rearranged vector/grid unless SKIP = TRUE
% NCHANGED            - Number of correlations changed (*1/2 to account for symmetry if grid)
% li                  - page 1 linear indices of peaks to be exchanged
% M                   - misfit matrix (generally for debugging)

% CHECK ARRIVAL TIME VECTOR, FORCE TO COLUMN VECTOR, GET NUMBER OF RECORDS
if(~isvector(dt) || size(dt,2)~=1); error('dt must be column vector'); end
dt=dt(:); nr=length(dt);

% CHECK THAT CORRELATION MATRICES MATCH
[nrows,ncols,npeaks]=size(cv);
if(~isequal([nrows,ncols,npeaks],size(lv),size(pv))); error('correlation matrices not same size'); end
if(~isscalar(xcw) && ~isequal([nrows,ncols,npeaks],size(xcw))); error('xcw bad size'); end
if(~isscalar(xcpow) && ~isequal([nrows,ncols,npeaks],size(xcpow))); error('xcpow bad size'); end

% GENERATE MISFIT FUNCTION (HIGH VALUES ARE WORSE)
%  - MISFIT OF LAGS TO INVERTED RELATIVE ARRIVALS
%  - WEIGHTED BY A COEFFICIENT AND INVERSILY BY THE CORRELATION TO SOME POWER
if(any([nrows ncols]==1)) 
    % CORRELATION MATRICES ARE VECTORS
    if(length(cv)~=(nr^2-nr)/2); error('correlation vectors have bad length'); end
    [j,i]=ind2sub([nr nr],find(tril(true(nr),-1)));  % MATRIX SUBSCRIPTS FOR VECTORS
    M=abs((lv+dt(i,1,ones(1,npeaks))-dt(j,1,ones(1,npeaks))).*xcw./cv.^xcpow);
    factor=1;
else
    % CORRELATION MATRICES ARE GRIDS
    if(~isequal(nrows,ncols,nr)); error('correlation grids must be square'); end
    M=abs((lv+dt(:,ones(nr,1),ones(npeaks,1))-permute(dt(:,ones(nr,1),ones(npeaks,1)),[2 1 3])).*xcw./cv.^xcpow);
    factor=2;
end

% ONLY CHANGE THOSE THAT EXCEED THRESHOLD
thresh=max(max(M(:,:,1)))*frac;
M(M(:,:,ones(1,npeaks))<thresh)=0;

% FIND BEST PEAKS NOT ON PAGE 1 (COMPARES DOWN PAGES = 3RD DIM)
[v,mi]=min(M,[],3);
li=find(mi>1);
NCHANGED=size(li,1)/factor;

% QUICK EXIT
if (NCHANGED==0 || SKIP); return; end

% INFORM
disp(['CHANGING ' num2str(NCHANGED) ' CORRELATIONS'])

% MOVING BEST PEAKS NOT ON PAGE 1 TO TEMPORARY VECTOR
ct=cv(li+(mi(li)-1)*nr^2);
lt=lv(li+(mi(li)-1)*nr^2);
pt=pv(li+(mi(li)-1)*nr^2);

% MOVING PAGE 1 PEAKS THAT ARE NOT THE BEST
cv(li+(mi(li)-1)*nr^2)=cv(li);
lv(li+(mi(li)-1)*nr^2)=lv(li);
pv(li+(mi(li)-1)*nr^2)=pv(li);

% MOVING BEST PEAKS NOT ON PAGE 1 TO PAGE 1
cv(li)=ct;
lv(li)=lt;
pv(li)=pt;

end
