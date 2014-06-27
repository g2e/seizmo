function [w]=pairweights(st,ev,method,stddev)
%PAIRWEIGHTS    Returns weights for a coarray frequency slowness spectra
%
%    Usage:    w=pairweights(st)
%              w=pairweights(st,matrix)
%              w=pairweights(st,ev)
%              w=pairweights(st,ev|matrix,method)
%              w=pairweights(st,ev|matrix,'gaussian',stddev)
%
%    Description:
%     W=PAIRWEIGHTS(ST) returns weights useful for improving the resolution
%     of frequency slowness spectra estimates.  This calculates weights for
%     the N*(N-1)/2 pairings between stations given by location matrix ST.
%     ST must be a Nx2 matrix of [LAT LON] in degrees.  W is a N*(N-1)/2
%     element column vector.  This will only work with the 'coarray' method
%     of FSS or a similar function.
%
%     W=PAIRWEIGHTS(ST,MATRIX) selects for which matrix type the automatic
%     pairing is done.  Either 'full' or 'coarray' is allowed.  The default
%     is 'coarray'.  W is N^2 element column vector if MATRIX='full'.
%
%     W=PAIRWEIGHTS(ST,EV) uses the pairings explicitly defined by ST & EV.
%     Both are Nx2 matrices of [LAT LON].  W is an Nx1 column vector.  This
%     is useful for more complicated sets of pairs (missing pairs, etc).
%
%     W=PAIRWEIGHTS(ST,EV|MATRIX,METHOD) sets the pair weight method.  The
%     available choices are:
%      'cartesian' - weight is minimum (non-zero) cartesian distance
%      'polar'     - weight is based on minimum (non-zero) radius
%                    difference multiplied by the normalized cross product
%      'uniform'   - weights approximate uniformly distributed array
%      'gaussian'  - weights approximate normally distributed array
%     The default is METHOD='gaussian'.
%
%     W=PAIRWEIGHTS(ST,EV|MATRIX,'GAUSSIAN',STDDEV) alters the number of
%     standard deviations of the gaussian distribution that the array
%     samples.  The default is STDDEV=2.5.  High values give lower slowness
%     resolution with some tradeoff for lower noise levels.
%
%    Notes:
%     - Pair weights are determined by the uniqueness of the individual
%       spatial difference vectors (which form the coarray).  Explicitly
%       the weight for a specific pair is based on the minimum distance
%       between the corresponding spatial difference vector and the
%       remainder of the coarray.
%
%    Examples:
%     % The 0.1Hz array response of LASA:
%     plotarf(arf(lasa,50,201,0,0,1/10,'m','coarray'));
%     % and now compare that to a weighted version:
%     plotarf(arf(lasa,50,201,0,0,1/10,...
%         'm','coarray','w',pairweights(lasa)));
%
%    See also: FSS, FSSXC, FSSHORZ, FSSHORZXC, ARF, ARFHORZ

%     Version History:
%        Sep. 27, 2012 - initial version
%        Mar. 26, 2014 - added uniform/normal array weights
%        Mar. 27, 2014 - sdv clustering, lower sdv grid resolution,
%                        gaussian is now default, stddev option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 27, 2014 at 15:05 GMT

% todo:

% check number of inputs
error(nargchk(1,4,nargin));

% defaults
if(nargin<2 || isempty(ev)); ev='coarray'; end
if(nargin<3 || isempty(method)); method='gaussian'; end
if(nargin<4 || isempty(stddev)); stddev=2.5; end

% checks
if(~isnumeric(st) || ~isreal(st) || size(st,2)<2 || ndims(st)>2)
    error('seizmo:pairweights:badInput',...
        'ST must be given as a Nx2 array of [LAT LON] in degrees!');
elseif(~ischar(ev) && ~isnumeric(ev))
    error('seizmo:pairweights:badInput',...
        'Second argument must be METHOD or EV!');
elseif(~ischar(method))
    error('seizmo:pairweights:badInput',...
        'METHOD must be a string!');
elseif(~isnumeric(stddev) || ~isreal(stddev) ...
        || ~isscalar(stddev) || stddev<=0)
    error('seizmo:pairweights:badInput',...
        'STDDEV must be a real scalar >0!');
end

% determine pairs automatically or given explicitly
if(ischar(ev)) % auto pairing (for fss, fsshorz, arf, arfhorz)
    % get indices
    nrecs=size(st,1);
    switch lower(ev)
        case {'coarray' 'co' 'ca' 'c'}
            [master,slave]=find(triu(true(nrecs),1));
            [master2,slave2]=find(tril(true(nrecs),-1));
            master=[master; master2];
            slave=[slave; slave2];
        case {'full' 'f' 'capon'}
            [master,slave]=find(true(nrecs));
        otherwise
            error('seizmo:pairweights:badInput',...
                'Unknown method: %s !',ev);
    end
    npairs=numel(master);
    
    % projected locations
    [clat,clon]=arraycenter(st(:,1),st(:,2));
    [e,n]=geographic2enu(st(:,1),st(:,2),0,clat,clon,0);
    radius=arrayradius(st(:,1),st(:,2));
    
    % spatial difference vectors
    e=e(slave)-e(master);
    n=n(slave)-n(master);
else % explicit pairings (for fssxc, fsshorzxc)
    % check & expand location inputs
    if(~isnumeric(ev) || ~isreal(ev) || size(ev,2)<2 || ndims(ev)>2)
        error('seizmo:pairweights:badInput',...
            'EV must be given as a Nx2 array of [LAT LON] in degrees!');
    end
    if(size(st,1)==1); st=st(ones(size(ev,1),1),:); end % expand
    if(size(ev,1)==1); ev=ev(ones(size(st,1),1),:); end % scalars
    if(size(st,1)~=size(ev,1))
        error('seizmo:pairweights:badInput',...
            'ST & EV must be equally sized arrays!');
    end
    
    % get indices
    npairs=size(st,1);
    [st0,idx1,idx2]=unique([st;ev],'rows');
    idx2=reshape(idx2,[],2);
    slave=idx2(:,1);
    master=idx2(:,2);
    
    % projected locations
    [clat,clon]=arraycenter(st0(:,1),st0(:,2));
    [e,n]=geographic2enu(st0(:,1),st0(:,2),0,clat,clon,0);
    radius=arrayradius(st0(:,1),st0(:,2));
    
    % spatial difference vectors
    e=e(slave)-e(master);
    n=n(slave)-n(master);
end

% find weights for sdv
switch lower(method)
    case {'cartesian' 'cart' 'c'}
        % calc all spatial difference vector differences
        d=sqrt((e(:,ones(npairs,1))-e(:,ones(npairs,1))').^2 ...
            +(n(:,ones(npairs,1))-n(:,ones(npairs,1))').^2);
        
        % account for 0s
        nz=sum(d==0);
        d(d==0)=inf;
        
        % weights are minimum non-zero spatial
        % difference vector difference distance
        % - divided by number of zero "sdvdd" (at least 1 for auto sdvdd)
        w=min(d)./nz;
        
        % weighting statistics to suppress extreme weighting
        medw=median(log10(w));
        meddev=median(abs(log10(w)-medw));
        w(log10(w)-medw>3*meddev)=10^(medw+3*meddev);
    case {'polar' 'pol' 'p'}
        % calc all spatial difference vector differences
        r=sqrt(e.^2+n.^2);
        e=e./r; n=n./r;
        e(r==0)=0; n(r==0)=1;
        r=r(:,ones(npairs,1));
        d=abs((e*n'-n*e').*(r-r'));
        
        % account for 0s
        nz=sum(d==0);
        d(d==0)=inf;
        
        % weights are minimum non-zero spatial
        % difference vector difference distance
        % - divided by number of zero "sdvdd" (at least 1 for auto sdvdd)
        w=min(d)./nz;
        
        % weighting statistics to suppress extreme weighting
        medw=median(log10(w));
        meddev=median(abs(log10(w)-medw));
        w(log10(w)-medw>3*meddev)=10^(medw+3*meddev);
    case {'uniform' 'uni' 'u'}
        % calc all spatial difference vector differences
        d=sqrt((e(:,ones(npairs,1))-e(:,ones(npairs,1))').^2 ...
            +(n(:,ones(npairs,1))-n(:,ones(npairs,1))').^2);
        
        % nonzero minimum of sdvdd to get sdv grid interval
        maxres=100;
        d(d<max(hypot(e,n))/maxres)=inf;
        di=min(d(:));
        
        % get clusters
        [idx1,idx2]=find(isinf(d));
        gidx=1:npairs; gidxnew=gidx;
        done=false(size(idx1));
        gmem=cell(npairs,1);
        for i=1:numel(idx1)
            if(done(i)); continue; end
            in=idx1==idx1(i) | idx2==idx1(i);
            g=unique([idx1(in); idx2(in)]);
            gidxnew(g)=idx1(i);
            while(~isequal(gidx,gidxnew))
                gidx=gidxnew;
                in=ismember(idx1,g) | ismember(idx2,g);
                g=unique([idx1(in); idx2(in)]);
                gidxnew(g)=idx1(i);
            end
            done(in)=true;
            [gmem{g}]=deal(g);
        end
        
        % create sampling grid
        gi=di;
        maxg=gi*ceil(max(hypot(e,n))/gi);
        [ge,gn]=meshgrid(-maxg:gi:maxg,-maxg:gi:maxg);
        
        % sample cone with grid
        f=1-hypot(ge,gn)/maxg;
        f(f<0)=0;
        
        % determine which voronoi cell each grid point is in
        dmin=inf(size(ge));
        di=zeros(size(ge));
        for i=1:npairs
            d=sqrt((ge-e(i)).^2+(gn-n(i)).^2);
            di(d<dmin)=i;
            dmin(d<dmin)=d(d<dmin);
        end
        
        % compute weights
        w=zeros(npairs,1);
        for i=1:npairs
            w(i)=sum(f(ismember(di,gmem{i})))/numel(gmem{i});
        end
    case {'gaussian' 'gauss' 'g' 'normal' 'norm' 'n'}
        % calc all spatial difference vector differences
        d=sqrt((e(:,ones(npairs,1))-e(:,ones(npairs,1))').^2 ...
            +(n(:,ones(npairs,1))-n(:,ones(npairs,1))').^2);
        
        % nonzero minimum of sdvdd to get sdv grid interval
        maxres=100;
        d(d<max(hypot(e,n))/maxres)=inf;
        di=min(d(:));
        
        % get clusters
        [idx1,idx2]=find(isinf(d));
        gidx=1:npairs; gidxnew=gidx;
        done=false(size(idx1));
        gmem=cell(npairs,1);
        for i=1:numel(idx1)
            if(done(i)); continue; end
            in=idx1==idx1(i) | idx2==idx1(i);
            g=unique([idx1(in); idx2(in)]);
            gidxnew(g)=idx1(i);
            while(~isequal(gidx,gidxnew))
                gidx=gidxnew;
                in=ismember(idx1,g) | ismember(idx2,g);
                g=unique([idx1(in); idx2(in)]);
                gidxnew(g)=idx1(i);
            end
            done(in)=true;
            [gmem{g}]=deal(g);
        end
        
        % create sampling grid
        gi=di;
        maxg=gi*ceil(max(hypot(e,n))/gi);
        [ge,gn]=meshgrid(-maxg:gi:maxg,-maxg:gi:maxg);
        
        % sample gaussian with grid
        f=gaussian(hypot(ge,gn),0,sqrt(2)*radius/stddev);
        f(hypot(ge,gn)>maxg)=0;
        
        % determine which voronoi cell each grid point is in
        dmin=inf(size(ge));
        di=zeros(size(ge));
        for i=1:npairs
            d=sqrt((ge-e(i)).^2+(gn-n(i)).^2);
            di(d<dmin)=i;
            dmin(d<dmin)=d(d<dmin);
        end
        
        % compute weights
        w=zeros(npairs,1);
        for i=1:npairs
            w(i)=sum(f(ismember(di,gmem{i})))/numel(gmem{i});
        end
    otherwise
        error('seizmo:pairweights:badInput',...
            'Unknown pair weighting method: %s !',method);
end

% remove redundant coarray weights
if(ischar(ev))
    switch lower(ev)
        case {'coarray' 'co' 'ca' 'c'}
            w=w(1:(npairs/2));
    end
end

end
