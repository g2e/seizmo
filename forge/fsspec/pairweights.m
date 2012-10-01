function [w]=pairweights(st,ev,distance)
%PAIRWEIGHTS    Returns weights for a coarray frequency slowness spectra
%
%    Usage:    w=pairweights(st)
%              w=pairweights(st,method)
%              w=pairweights(st,ev)
%              w=pairweights(st,ev|method,distance)
%
%    Description:
%     W=PAIRWEIGHTS(ST) returns weights useful for improving the resolution
%     of frequency slowness spectra estimates.  This calculates weights for
%     the N*(N-1)/2 pairings between stations given by location matrix ST.
%     ST must be a Nx2 matrix of [LAT LON] in degrees.  W is a N*(N-1)/2
%     element column vector.  This will only work with the 'coarray' method
%     of FSS or a similar function.
%
%     W=PAIRWEIGHTS(ST,METHOD) selects for which method the automatic
%     pairing is done.  Either 'full', 'coarray' or 'capon' is allowed.
%     The default is 'coarray'.
%
%     W=PAIRWEIGHTS(ST,EV) uses the pairings explicitly defined by ST & EV.
%     Both are Nx2 matrices of [LAT LON].  W is an Nx1 column vector.  This
%     is useful for more complicated sets of pairs (missing pairs, etc).
%
%     W=PAIRWEIGHTS(ST,EV|METHOD,DISTANCE) sets the pair weight distance
%     metric.  The default is 'cartesian'.  An alternative is 'polar' which
%     uses the radius and a normalized cross product as a distance measure.
%
%    Notes:
%     - Pair weights are determined by the uniqueness of the individual
%       spatial difference vectors (which form the coarray).  Explicitly
%       the weight for a specific pair is the minimum distance between the
%       corresponding spatial difference vector and the remainder of the
%       coarray.  This does an excellent job of narrowing the central peak
%       while reducing nearby grating lobes but does enhance the background
%       level and grating lobes further away from the central peak.
%
%    Examples:
%     % The 0.1Hz array response of LASA:
%     plotarf(arf(lasa,50,201,0,0,1/10));
%     % and now compare that to a weighted cross-spectra version:
%     plotarf(arf(lasa,50,201,0,0,1/10,...
%         'm','coarray','w',pairweights(lasa)));
%
%    See also: FSS, FSSXC, FSSHORZ, FSSHORZXC, ARF, ARFHORZ

%     Version History:
%        Sep. 27, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 27, 2012 at 15:05 GMT

% todo:

% check number of inputs
error(nargchk(1,3,nargin));

% defaults
if(nargin<2 || isempty(ev)); ev='coarray'; end
if(nargin<3 || isempty(distance)); distance='cartesian'; end

% checks
if(~isnumeric(st) || ~isreal(st) || size(st,2)<2 || ndims(st)>2)
    error('seizmo:pairweights:badInput',...
        'ST must be given as a Nx2 array of [LAT LON] in degrees!');
elseif(~ischar(ev) && ~isnumeric(ev))
    error('seizmo:pairweights:badInput',...
        'Second argument must be METHOD or EV!');
elseif(~ischar(distance))
    error('seizmo:pairweights:badInput',...
        'DISTANCE must be a string!');
end

% determine pairs automatically or given explicitly
factor=0; % how many repeats of the same pairs
if(ischar(ev)) % auto pairing (for fss, fsshorz, arf, arfhorz)
    % get indices
    nrecs=size(st,1);
    switch ev
        case 'coarray'
            [master,slave]=find(triu(true(nrecs),1));
        case {'full' 'capon'}
            [master,slave]=find(true(nrecs));
            factor=1; % 1 repeat
        otherwise
            error('seizmo:pairweights:badInput',...
                'Unknown method: %s !',ev);
    end
    npairs=numel(master);
    
    % projected locations
    [clat,clon]=arraycenter(st(:,1),st(:,2));
    [e,n]=geographic2enu(st(:,1),st(:,2),0,clat,clon,0);
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
    e=e(slave)-e(master);
    n=n(slave)-n(master);
end

% flip sign of sdv with negative north
% - this chooses one direction between each station pair
%   so differences don't favor one direction over another
e(n<0)=-e(n<0);
n(n<0)=-n(n<0);

% calc all spatial difference vector differences
switch distance
    case 'cartesian'
        d=sqrt((e(:,ones(npairs,1))-e(:,ones(npairs,1))').^2 ...
            +(n(:,ones(npairs,1))-n(:,ones(npairs,1))').^2);
    case 'polar'
        r=sqrt(e.^2+n.^2);
        e=e./r; n=n./r;
        e(r==0)=0; n(r==0)=1;
        %r=log10(r);
        r=r(:,ones(npairs,1));
        d=abs((e*n'-n*e').*(r-r'));
    otherwise
        error('seizmo:pairweights:badInput',...
            'Unknown pair weighting method: %s !',method);
end


% account for 0s (& -inf)
nz=sum(d==0 | d==-inf)-factor;
d(d==0 | d==-inf)=inf;

% weights are minimum non-zero spatial
% difference vector difference distance
% - divided by number of zero "sdvdd" (at least 1 for auto sdvdd)
w=min(d)./nz;

% weighting statistics to suppress extreme weighting
medw=median(log10(w));
meddev=median(abs(log10(w)-medw));
w(log10(w)-medw>3*meddev)=10^(medw+3*meddev);

end
