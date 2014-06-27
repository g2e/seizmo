function [data]=reverse_correlations(data)
%REVERSE_CORRELATIONS    Reverses correlations (switch master & slave)
%
%    Usage:    data=reverse_correlations(data)
%
%    Description:
%     DATA=REVERSE_CORRELATIONS(DATA) does a time reversal of correlograms
%     in SEIZMO struct DATA and switches the master and slave header field
%     info as well as the filenames.  This is only for correlations
%     generated with CORRELATE.
%
%    Notes:
%     - Requires that all records have the fields KUSER0 & KUSER1 set to
%       'MASTER' & 'SLAVE'.
%
%    Header changes:
%     Switches B & E (for timeseries records):
%      B = -E
%      E = -B
%     Master & Slave Field Switches:
%      ST   <-> EV
%      KNETWK <-> KT0
%      KSTNM  <-> KT1
%      KHOLE  <-> KT2
%      kCMPNM <-> KT3
%      USER0  <-> USER1
%      CMPINC <-> USER2
%      CMPAZ  <-> USER3
%      A,F    <-> T0,T1
%      NXSIZE <-> NYSIZE
%     Updates:
%      GCARC, AZ, BAZ, DIST, Z, T3, T4
%
%    Examples:
%     % Get all cross correlation pairs (only computing half of them):
%     corr0=correlate(data); % auto correlations
%     corr1=correlate(data,'mcxc','noauto');
%     corr2=reverse_correlations(corr1);
%
%    See also: CORRELATE, REVERSE, SPLIT_AUTO_CORRELATIONS, ISXC,
%              ROTATE_CORRELATIONS, NAME_CORRELATIONS,
%              NO_REDUNDANT_CORRELATIONS, HORZ_CORRELATIONS_SETS,
%              IS_FULL_MATRIX_OF_CORRELATIONS

%     Version History:
%        Apr. 12, 2010 - initial version
%        Apr. 15, 2010 - resets name field, includes scrollbar if verbose
%        Feb. 11, 2011 - mass nargchk fix
%        Jan. 27, 2012 - update for cmpinc/cmpaz fields, better checkheader
%                        usage
%        Oct. 21, 2012 - update for a,f,t0,t1 fields
%        Sep.  9, 2013 - stricter checkheader call
%        Sep. 20, 2013 - updated See also section
%        May  30, 2014 - use different renaming for stacks
%        June 12, 2014 - handle fd i/o
%        June 16, 2014 - bugfix: update sb for fd i/o
%        June 25, 2014 - bugfix: iftype id in getheader call
%        June 26, 2014 - now z,t3,t4 fields are updated, fixed verbosity
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 26, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
    
% hide verbosity from underlying functions
verbose=seizmoverbose(false);

% attempt reversal
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'XYZ_IFTYPE','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    
    % number of records
    nrecs=numel(data);
    
    % detail message
    if(verbose); disp('Reversing Correlogram(s)'); end
    
    % get header info
    [kuser0,kuser1,user0,user1,user2,user3,cmpinc,cmpaz,st,ev,knetwk,...
        kstnm,khole,kcmpnm,kt0,kt1,kt2,kt3,b,e,autc,futc,t0utc,t1utc,...
        a,f,t3,nx,ny,iftype,sb,sdelta,nspts]=getheader(data,'kuser0',...
        'kuser1','user0','user1','user2','user3','cmpinc','cmpaz','st',...
        'ev','knetwk','kstnm','khole','kcmpnm','kt0','kt1','kt2','kt3',...
        'b','e','a utc','f utc','t0 utc','t1 utc','a','f','t3','nxsize',...
        'nysize','iftype id','sb','sdelta','nspts');
    
    % correlogram signature
    if(~all(strcmp('MASTER',kuser0) & strcmp('SLAVE',kuser1)))
        error('seizmo:reverse_correlations:notCorrelogram',...
            'Correlograms appear to have malformed headers!');
    end
    
    % filetypes
    rlim=strcmpi(iftype,'irlim');
    amph=strcmpi(iftype,'iamph');
    time=~rlim & ~amph;
    
    % detail message
    if(verbose); print_time_left(0,nrecs); end
    
    % get slave reftime
    z=cell2mat(t0utc);
    z=mat2cell(fixtimes([z(:,1:end-1) z(:,end)-t3],'utc'),ones(nrecs,1));
    
    % fix header info
    data=changeheader(data,'st',ev,'ev',st,'user0',user1,...
        'user1',user0,'knetwk',kt0,'kstnm',kt1,'khole',kt2,'kcmpnm',kt3,...
        'kt0',knetwk,'kt1',kstnm,'kt2',khole,'kt3',kcmpnm,...
        'cmpinc',user2,'cmpaz',user3,'user2',cmpinc,'user3',cmpaz,...
        'z',z,'t3',a,'t4',f,'nxsize',ny,'nysize',nx,...
        'a utc',t0utc,'f utc',t1utc,'t0 utc',autc,'t1 utc',futc);
    
    % reverse data
    % - just complex conjugate in frequency domain
    if(any(rlim)); data(rlim)=solofun(data(rlim),@(x)[x(:,1) -x(:,2)]); end
    if(any(amph))
        data(amph)=amph2rlim(data(amph));
        data(amph)=solofun(data(amph),@(x)[x(:,1) -x(:,2)]);
        data(amph)=rlim2amph(data(amph));
    end
    if(any(rlim | amph))
        data(rlim | amph)=changeheader(data(rlim | amph),'sb',...
            -sb(rlim | amph)-sdelta(rlim | amph).*(nspts(rlim | amph)-1));
    end
    if(any(time))
        data(time)=reverse(data(time));
        data(time)=changeheader(data(time),'b',-e(time),'e',-b(time));
    end
    
    % update delaz info
    oldstate=checkheader_state(true);
    data=checkheader(data,'all','ignore','old_delaz','fix');
    checkheader_state(oldstate);
    
    % rename
    short=~strncmp('CORR_',{data.name}',5);
    digits=['%0' num2str(fix(log10(max([user0; user1])))+1) 'd'];
    name=strcat('CORR_-_MASTER_-_REC',num2str(user1,digits),'_-_',...
        knetwk,'.',kstnm,'.',khole,'.',kcmpnm,'_-_SLAVE_-_REC',...
        num2str(user0,digits),'_-_',kt0,'.',kt1,'.',kt2,'.',kt3);
    if(any(short))
        name(short)=strcat(knetwk(short),'.',kstnm(short),'.',...
            khole(short),'.',kcmpnm(short),'_-_',kt0(short),'.',...
            kt1(short),'.',kt2(short),'.',kt3(short));
    end
    [data.name]=deal(name{:});
    
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % rethrow error
    error(lasterror);
end

end
