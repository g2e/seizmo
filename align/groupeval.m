function [A,CG,LG,PG,li,M]=groupeval(T,cv,lv,pv,dt,nsd,mri,xcpow,xcw,frac)
%GROUPEVAL    Solves relative arrival times for each waveform group

% ADDITIONAL INPUTS
%nsd=2;   % number of standard deviations for confidence intervals
%mri=20;  % max refinement iterations
%xcpow=2; % weighting power of xc values
%xcw=1;   % weighting coeff. of xc values

% NUMBER OF RECORDS
nrecs=numel(T);

% CHECKS
if(~isequal(size(cv),size(lv),size(pv)))
    error('correlation matrices not same size'); 
end

% MAKE INTO GRIDS
sz=size(cv); np=prod(sz(3:end));
if(isvector(cv(:,:,1)))
    cv=ndsquareform(cv);
    lv=ndsquareform(lv,false); % ANTI-SYMMETRIC
    pv=ndsquareform(pv);
end

% CHECKS
if(~isequal(nrecs,size(cv,1)))
    error('correlation matrices dont match group vector'); 
end

% CORRELATION VALUES TO Z-STATISTIC (FOR NORMAL DISTRIBUTION)
cv(abs(cv)>1-eps)=nan; % IGNORE PERFECT CORRELATIONS
zv=fisher(cv);

% ALLOCATE RESULTS MATRIX
%
% COLUMNS:
% 1 - GROUP
% 2 - ARRIVAL
% 3 - ARRIVAL ERROR
% 4 - POLARITY
% 5 - CORR LOWER BOUND
% 6 - CORR MEAN
% 7 - CORR UPPER BOUND
A=nan(nrecs,7);
A(:,1)=T;

% GROUP INFO
ngrps=max(T);               % NUMBER OF GROUPS
gi=cell(ngrps,1);           % ALLOCATION
gpop=zeros(ngrps,1);        % ALLOCATION
for i=1:ngrps
    gi{i}=find(T==i);       % INDICES OF GROUP i MEMBERS
    gpop(i)=length(gi{i});  % POPULATION OF GROUP i
end

% BIG GROUPS ONLY (2+ MEMBERS)
big=find(gpop>1);
nbg=length(big);

% SEPARATE GROUPS
CG=cell(nbg,1);
LG=CG; PG=CG;
for i=1:nbg
    % MAKES NDIMS<=3
    CG{i}=zv(gi{big(i)},gi{big(i)},:)+repmat(eye(gpop(big(i))),[1 1 np]);
    LG{i}=lv(gi{big(i)},gi{big(i)},:);
    PG{i}=pv(gi{big(i)},gi{big(i)},:)+repmat(eye(gpop(big(i))),[1 1 np]);
end

% how are we gonna do this?
%  1. invert for relative arrival solution (need to add reweight option)
%  2. gently converge on a more consistent solution (refinement limited to outliers)
%  3. strongly converge on a more consistent solution
%  4. get statistics
%  5. get polarities (note the solution may not be consistent in
%  regard to polarities - this needs to be worked on)

% INITIAL INVERSION
% - DO WE ALIGN WEIGHTED TO THE ESTIMATED ARRIVAL TIMES
% - OR DO WE ALIGN BASED ON BEST CORRELATING PEAKS
disp('INVERTING FOR INTRAGROUP RELATIVE ARRIVALS')
if(isempty(dt))
    for i=1:nbg
        A(gi{big(i)},2)=dtwalign(CG{i}(:,:,1),LG{i}(:,:,1),xcpow,xcw);
    end
else
    % MAKE SURE SMALL GROUPS STAY AS NAN
    for i=1:nbg
        A(gi{big(i)},2)=dt(gi{big(i)});
    end
end

% LIGHTLY REFINE ARRIVAL TIMES USING THE MULTI-CANDIDATE METHOD
%
%  - CHOOSES CANDIDATE LAGS MOST CONSISTENT WITH INVERTED ARRIVAL TIMES
%  - LIMITS REFINEMENT TO OUTLIERS, LIMITING THEIR INFLUENCE
%  - HELPS AUTOMATICALLY DEAL WITH CYCLE SKIPPING
%
nri=zeros(nbg,2);
change=zeros(nbg,2*mri);
for i=1:nbg
    disp(['WEAK REFINING OF GROUP ' num2str(big(i)) '''S RELATIVE ARRIVALS'])
    for j=1:mri
        [CG{i},LG{i},PG{i},change(i,j),li,M]=...
            dtwrefine_weak(A(gi{big(i)},2),CG{i},LG{i},PG{i},xcpow,xcw,frac,0);
        
        % NO CHANGE - BREAK EARLY
        if (change(i,j)==0)
            nri(i)=j-1;
            break;
        end
        
        % CHANGED - REINVERT
        disp('REINVERTING')
        A(gi{big(i)},2)=dtwalign(CG{i}(:,:,1),LG{i}(:,:,1),xcpow,xcw);
        nri(i)=j; % INCREMENT COUNTER
        
        % DETECT NON-CONVERGENCE (AFTER AT LEAST 3 ITERATIONS)
        %
        % DETECTS TWO CASES: 
        %  1. NUMBER CHANGED STAYS CONSTANT OVER THREE ITERATIONS
        %  2. NUMBER CHANGED ALTERNATES BETWEEN 2 NUMBERS TWICE
        %
        if ((j>=3 && change(i,j)==change(i,j-1)...
                && change(i,j)==change(i,j-2))...
                || (j>=4 && change(i,j)==change(i,j-2)...
                && change(i,j-1)==change(i,j-3)))
            
            disp('NON-CONVERGENCE DETECTED')
            
            % CAN NOT RELY ON NON-CONVERGING TERMS
            % WEIGHT THEM DOWN AND SOLVE ONCE MORE
            disp('DOWN WEIGHTING TROUBLED CORRELATIONS AND REINVERTING')
            [CG{i},LG{i},PG{i},change(i,j),li,M]=...
                dtwrefine_weak(A(gi{big(i)},2),CG{i},LG{i},PG{i},xcpow,xcw,frac,1);
            CG{i}(li)=0; % NO WEIGHT
            A(gi{big(i)},2)=dtwalign(CG{i}(:,:,1),LG{i}(:,:,1),xcpow,xcw);
            break;
        end
    end
end

% REFINE ARRIVAL TIMES USING THE MULTI-CANDIDATE METHOD
%
%  - CHOOSES CANDIDATE LAGS MOST CONSISTENT WITH INVERTED ARRIVAL TIMES
%  - HELPS AUTOMATICALLY DEAL WITH CYCLE SKIPPING
%
nri=zeros(nbg,1);
change=zeros(nbg,mri);
for i=1:nbg
    disp(['STRONG REFINING OF GROUP ' num2str(big(i)) '''S RELATIVE ARRIVALS'])
    for j=1:mri
        [CG{i},LG{i},PG{i},change(i,j),li,M]=...
            dtwrefine_strong(A(gi{big(i)},2),CG{i},LG{i},PG{i},xcpow,xcw,0);
        
        % NO CHANGE - BREAK EARLY
        if (change(i,j)==0)
            nri(i)=j-1;
            break;
        end
        
        % CHANGED - REINVERT
        disp('REINVERTING')
        A(gi{big(i)},2)=dtwalign(CG{i}(:,:,1),LG{i}(:,:,1),2,1);
        nri(i)=j; % INCREMENT COUNTER
        
        % DETECT NON-CONVERGENCE (AFTER AT LEAST 3 ITERATIONS)
        %
        % DETECTS TWO CASES: 
        %  1. NUMBER CHANGED STAYS CONSTANT OVER THREE ITERATIONS
        %  2. NUMBER CHANGED ALTERNATES BETWEEN 2 NUMBERS TWICE
        %
        if ((j>=3 && change(i,j)==change(i,j-1)...
                && change(i,j)==change(i,j-2))...
                || (j>=4 && change(i,j)==change(i,j-2)...
                && change(i,j-1)==change(i,j-3)))
            
            disp('NON-CONVERGENCE DETECTED')
            
            % CAN NOT RELY ON NON-CONVERGING TERMS
            % WEIGHT THEM DOWN AND SOLVE ONCE MORE
            disp('DOWN WEIGHTING TROUBLED CORRELATIONS AND REINVERTING')
            [CG{i},LG{i},PG{i},change(i,j),li,M]=...
                dtwrefine_strong(A(gi{big(i)},2),CG{i},LG{i},PG{i},xcpow,xcw,1);
            CG{i}(li)=0; % NO WEIGHT
            A(gi{big(i)},2)=dtwalign(CG{i}(:,:,1),LG{i}(:,:,1),xcpow,xcw);
            break;
        end
    end
end

% GETTING ERROR BOUNDS ON RELATIVE ARRIVAL TIMES
disp('FINDING RELATIVE ARRIVAL ERROR BOUNDS')
for i=1:nbg
    A(gi{big(i)},3)=...
        dtwresid(A(gi{big(i)},2),CG{i}(:,:,1),LG{i}(:,:,1),2)*nsd;
end

% GETTING AVERAGE CORRELATIONS AND STANDARD DEVIATIONS
disp('EVALUATING CORRELATION STATISTICS')
for i=1:nbg
    Zmn=nanmean(CG{i}(:,:,1));
    Zstd=sqrt(nanvar(CG{i}(:,:,1)));
    A(gi{big(i)},6)=ifisher(Zmn);
    A(gi{big(i)},5)=ifisher(Zmn-Zstd*nsd);
    A(gi{big(i)},7)=ifisher(Zmn+Zstd*nsd);
end

% FIXING POLARITIES
%  - CURRENTLY JUST USES A MASTER TRACE TO DO THIS
%
%  - SHOULD MAKE A MORE ADVANCED METHOD THAT TIES IN WITH DTWREFINE TO
%    ASSURE POLARITY CONSISTENCY.  WOULD BREAK DOWN TO AN ALGORITHM
%    FOR CONVERTING PG TO A MATRIX OF ONES BY MULTIPLYING CERTAIN
%    ROW/COLUMN PAIRS BY -1.  THOSE PAIRS GIVE THE RECORDS WITH 
%    NEGATIVE POLARITY.  THIS IS NON-LINEAR THOUGH...
%
disp('SOLVING FOR BEST-FIT POLARITIES')
for i=1:nbg
    % USE RECORD WITH BEST AVG CORR IN GROUP
    [cm,seed]=max(A(gi{big(i)},6));
    A(gi{big(i)},4)=PG{i}(seed,:,1);
end


end
