function [s]=tdxcmigrate(data,ll,slow,w)
%TDXCMIGRATE    Migrate time-domain cross-correlations
%
%    Usage:    s=tdxcmigrate(xcdata,latlon,slow)
%              s=tdxcmigrate(xcdata,latlon,slow,weights)
%
%    Description:
%     S=TDXCMIGRATE(XCDATA,LATLON,SLOW) performs migration of cross
%     correlograms in the time domain following Brzak et al 2009.  The
%     Hilbert transform on the correlograms is skipped as it is not a
%     straightforward operation.  This gives quite similar results in
%     comparison to GEOFSS algorithms.  The output struct S matches the
%     GEOFSS struct format and is compatible with those functions.
%
%     S=TDXCMIGRATE(XCDATA,LATLON,SLOW,WEIGHTS) specifies the relative
%     weights for each correlogram in XCDATA (must match size of XCDATA).
%
%    Notes:
%     - This impliments the migration algorithm used in:
%        Brzak et al 2009, g3, doi:10.1029/2008GC002234
%        Burtin et al 2010, GJI, doi:10.1111/j.1365-246X.2010.04701.x
%        Gu & Shen 2012, BSSA, doi:10.1785/0120100010
%       With the following modifications:
%        - no hilbert transform
%        - migration amplitudes are forced positive
%     - This gives nearly equivalent results to GEOFSSXC (mainly a static
%       shift in dB values) with some careful filters & WHITEN=FALSE.
%
%    Examples:
%     % Compare GEOFSSXC & TDXCMIGRATE:
%     xcdata=iirfilter(xcdata,'bp','b','c',[1/27 1/25],'o',4,'p',2);
%     [lat,lon]=meshgrid(-10:.25:10,-5:.25:15);
%     t=tdxcmigrate(xcdata,[lat(:) lon(:)],29);
%     f=geofssxc(xcdata,[lat(:) lon(:)],29,[1/27 1/25]);
%     plotgeofss(geofssavg(t));
%     plotgeofss(geofssavg(f));
%
%    See also: GEOFSSXC, PLOTGEOFSS, GEOFSSAVG, GEOFSSSUB

%     Version History:
%        June 12, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2012 at 14:05 GMT

% todo:
% - test the effect of hilbert transforming 2-way greens function

% check nargin
error(nargchk(3,4,nargin));

% check struct
error(seizmocheck(data,'dep'));

% number of correlograms
ncorr=numel(data);

% defaults for optionals
if(nargin<4 || isempty(w)); w=ones(numel(data),1); end

% check inputs
if(~isreal(ll) || ndims(ll)~=2 || size(ll,2)~=2)
    error('seizmo:tdxcmigrate:badInput',...
        'LATLON must be a Nx2 real matrix of [LAT LON]!');
elseif(~isreal(slow) || ~isvector(slow) || any(slow<=0))
    error('seizmo:tdxcmigrate:badInput',...
        'SLOW must be a positive real vector in sec/deg!');
elseif(numel(w)~=ncorr || any(w(:)<0) || ~isreal(w))
    error('seizmo:tdxcmigrate:badInput',...
        'WEIGHTS must be equal sized with XCDATA & be positive numbers!');
end

% column vector slownesses
slow=slow(:);
nslow=numel(slow);

% fix lat/lon
[ll(:,1),ll(:,2)]=fixlatlon(ll(:,1),ll(:,2));
nll=size(ll,1);

% convert weights to column vector and normalize
w=w(:);
w=w./sum(w);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% migrating correlograms
try
    % verbosity
    verbose=seizmoverbose;
    
    % hilbert transform
    % how to -pi/2 phase shift both forward & reverse green's functions?
    % - we could split by +/- time, hilbert & rejoin (test this!)
    %data=hilbrt(data);
    
    % get some useful info from headers
    [b,npts,delta,st,ev,a,f]=getheader(data,...
        'b','npts','delta','st','ev','a utc','f utc');
    delta=delta(1);
    fnyq=1/(2*delta);
    
    % find unique station locations
    loc=unique([st; ev],'rows');
    nsta=size(loc,1);
    
    % start/end time range (handles old correlograms)
    a=cell2mat(a); a=a(all(~isnan(a),2),:);
    f=cell2mat(f); f=f(all(~isnan(f),2),:);
    if(~isempty(a))
        [ai,ai]=min(timediff(a,a(1,:),'utc'));
        a=a(ai,:);
    else
        a=[0 0 0 0 0];
    end
    if(~isempty(f))
        [fi,fi]=max(timediff(f,f(1,:),'utc'));
        f=f(fi,:);
    else
        f=[0 0 0 0 0];
    end
    
    % array center
    [clat,clon]=arraycenter([st(:,1); ev(:,1)],[st(:,2); ev(:,2)]);
    
    % setup output
    s.nsta=nsta;
    s.stla=loc(:,1);
    s.stlo=loc(:,2);
    s.stel=loc(:,3);
    s.stdp=loc(:,4);
    s.butc=a;
    s.eutc=f;
    s.delta=delta(1);
    s.npts=npts(1);
    s.vector=[false nslow~=1];
    s.latlon=ll;
    s.slow=slow;
    s.freq=[0; fnyq];
    s.method='center';
    s.npairs=ncorr;
    s.center=[clat clon];
    s.weights=w;
    s.spectra=zeros(nll,nslow,'single');
    
    % distance difference for the phasors that "steer" the array
    % dd is NLLxNCORR
    ev=ev.'; st=st.';
    distm=sphericalinv(ll(:,ones(ncorr,1)),ll(:,2*ones(ncorr,1)),...
        ev(ones(nll,1),:),ev(2*ones(nll,1),:));
    dists=sphericalinv(ll(:,ones(ncorr,1)),ll(:,2*ones(ncorr,1)),...
        st(ones(nll,1),:),st(2*ones(nll,1),:));
    dd=dists-distm;
    
    % detail message
    if(verbose)
        fprintf('Getting Migrated Amplitudes\n');
        print_time_left(0,nslow*ncorr);
    end
    
    % loop over slownesses
    for a=1:nslow
        % convert distance difference to time
        dt=dd.*slow(a);
        
        % now convert time to sample
        x=round((dt-b(:,ones(1,nll))')/delta);
        
        % loop over correlograms
        for c=1:ncorr
            ok=x(:,c)>0 & x(:,c)<=npts(c); % handle short correlograms
            s.spectra(:,a)=s.spectra(:,a)+data(c).dep(x(ok,c))*w(c);
            if(verbose)
                print_time_left((a-1)*ncorr+c,nslow*ncorr);
            end
        end
    end
    
    % force positive
    s.spectra=abs(s.spectra);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end
