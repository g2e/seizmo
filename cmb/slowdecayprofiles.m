function [pf]=slowdecayprofiles(results,azrng,gcrng)
%SLOWDECAYPROFILES    Returns multi-station profile measurements
%
%    Usage:    pf=slowdecayprofiles(results,azrng,gcrng)
%
%    Description:
%     PF=SLOWDECAYPROFILES(RESULTS,AZRNG,GCRNG) takes the relative arrival
%     time and amplitude measurements contained in RESULTS produced by
%     CMB_1ST_PASS, CMB_2ND_PASS, or CMB_OUTLIERS and calculates the
%     slowness and decay rate for a profile of stations within the
%     criteria set by azimuthal range AZRNG and distance range GCRNG.  Note
%     that AZRNG & GCRNG are absolute ranges, meaning an AZRNG of [0 360]
%     will not exclude any stations by azimuth.  AZRNG & GCRNG must both be
%     2-element vectors of [AZMIN AZMAX] & [GCMIN GCMAX].  They are by
%     default [0 360] & [0 180] (ie no exclusion) and are optional.
%
%    Notes:
%
%    Examples:
%     % Return station profiles with an azimuth
%     % of 200-220deg and a distance of 90-160deg:
%     pf=slowdecayprofiles(results,[200 220],[90 160])
%
%    See also: PREP_CMB_DATA, CMB_1ST_PASS, CMB_OUTLIERS, CMB_2ND_PASS,
%              SLOWDECAYPAIRS, CMB_CLUSTERING

%     Version History:
%        Dec. 12, 2010 - initial version
%        Jan. 18, 2011 - update for results struct standardization, added
%                        corrections & correlation coefficients to output,
%                        time is now a string, require common event
%        Jan. 23, 2011 - fix indexing bug
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2011 at 13:35 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check results struct
error(check_cmb_results(results));

% default azrng & gcrng
if(nargin<2 || isempty(azrng)); azrng=[0 360]; end
if(nargin<3 || isempty(gcrng)); gcrng=[0 180]; end

% check azrng & gcrng
if(~isreal(azrng) || numel(azrng)~=2)
    error('seizmo:slowdecayprofiles:badInput',...
        'AZRNG must be [AZMIN AZMAX]!');
elseif(any(abs(azrng)>540))
    error('seizmo:slowdecayprofiles:badInput',...
        'Keep AZRNG within +/-540deg!');
elseif(~isreal(gcrng) || numel(gcrng)~=2)
    error('seizmo:slowdecayprofiles:badInput',...
        'GCRNG must be [GCMIN GCMAX]!');
end

% verbosity
%verbose=seizmoverbose;

% loop over every result
total=0;
for a=1:numel(results)
    % skip if results.useralign is empty
    if(isempty(results(a).useralign)); continue; end
    
    % number of records
    %nrecs=numel(results(a).useralign.data);
    
    % extract header details
    [st,ev,delaz,kname]=getheader(results(a).useralign.data,...
        'st','ev','delaz','kname');
    
    % check event info matches
    ev=unique(ev,'rows');
    if(size(ev,1)>1)
        error('seizmo:slowdecaypairs:badInput',...
            'EVENT location varies between records!');
    end
    
    % corrected relative arrival times and amplitudes
    rtime=results(a).useralign.solution.arr;
    crtime=results(a).useralign.solution.arr...
        -results(a).corrections.ellcor...
        -results(a).corrections.crucor.prem...
        -results(a).corrections.mancor.hmsl06p.upswing;
    rtimeerr=results(a).useralign.solution.arrerr;
    rampl=results(a).useralign.solution.amp;
    crampl=results(a).useralign.solution.amp...
        ./results(a).corrections.geomsprcor;
    ramplerr=results(a).useralign.solution.amperr;
    
    % get cluster indexing
    cidx=results(a).usercluster.T;
    good=results(a).usercluster.good';
    
    % get outliers
    outliers=results(a).outliers.bad;
    
    % loop over "good" clusters
    for b=find(good)
        % indices of members
        m=cidx==b & ~outliers;
        
        % get stations in range (need the indices)
        delaz(:,2)=delaz(:,2)-360*ceil((delaz(:,2)-azrng(2))/360);
        idx=find(m & delaz(:,1)>=gcrng(1) & delaz(:,1)<=gcrng(2) ...
            & delaz(:,2)>=azrng(1) & delaz(:,2)<=azrng(2));
        nsta=numel(idx);
        
        % skip if none/one
        if(nsta<2)
            continue;
        else
            total=total+1;
        end
        
        % initialize struct
        pf(total)=struct('gcdist',[],'azwidth',[],...
            'slow',[],'slowerr',[],'decay',[],'decayerr',[],...
            'cslow',[],'cslowerr',[],'cdecay',[],'cdecayerr',[],...
            'cluster',b,'kname',[],'st',[],'ev',[],'delaz',[],...
            'corrections',[],'corrcoef',[],...
            'freq',results(a).filter.corners,'phase',results(a).phase,...
            'runname',results(a).runname,'dirname',results(a).dirname,...
            'time',datestr(now));
        
        % insert known info
        pf(total).kname=kname(idx,:);
        pf(total).st=st(idx,:);
        pf(total).ev=ev;
        pf(total).delaz=delaz(idx,:);
        
        % great circle distance and width
        pf(total).gcdist=max(delaz(idx,1))-min(delaz(idx,1));
        pf(total).azwidth=max(delaz(idx,2))-min(delaz(idx,2));
        
        % corrections
        pf(total).corrections=fixcorrstruct(results(a).corrections,idx);
        
        % correlation coefficients
        pf(total).corrcoef=...
            submat(ndsquareform(results(a).useralign.xc.cg),1:2,idx,3,1);
        
        % find slowness & decay rate
        [m,covm]=wlinem(delaz(idx,1),rtime(idx),1,diag(rtimeerr(idx).^2));
        pf(total).slow=m(2);
        pf(total).slowerr=sqrt(covm(2,2));
        [m,covm]=wlinem(delaz(idx,1),crtime(idx),1,diag(rtimeerr(idx).^2));
        pf(total).cslow=m(2);
        pf(total).cslowerr=sqrt(covm(2,2));
        [m,covm]=wlinem(delaz(idx,1),log(rampl(idx)),1,...
            diag(log(rampl(idx)+ramplerr(idx).^2)-log(rampl(idx))));
        pf(total).decay=m(2);
        pf(total).decayerr=sqrt(covm(2,2));
        [m,covm]=wlinem(delaz(idx,1),log(crampl(idx)),1,...
            diag(log(crampl(idx)+ramplerr(idx).^2)-log(crampl(idx))));
        pf(total).cdecay=m(2);
        pf(total).cdecayerr=sqrt(covm(2,2));
    end
end

end


function [s]=fixcorrstruct(s,good)
fields=fieldnames(s);
for i=1:numel(fields)
    if(isstruct(s.(fields{i})))
        s.(fields{i})=fixcorrstruct(s.(fields{i}),good);
    else
        s.(fields{i})=s.(fields{i})(good);
    end
end
end
