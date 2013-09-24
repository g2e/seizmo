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
%     Switches B & E :
%      B = -E
%      E = -B
%     Master & Slave Field Switches :
%      ST   <-> EV
%      KNETWK <-> KT0
%      KSTNM  <-> KT1
%      KHOLE  <-> KT2
%      kCMPNM <-> KT3
%      USER0  <-> USER1
%      CMPINC <-> USER2
%      CMPAZ  <-> USER3
%      A,F    <-> T0,T1
%     Updates:
%      GCARC, AZ, BAZ, DIST
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 20, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt reversal
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    
    % number of records
    nrecs=numel(data);
    
    % check verbosity
    verbose=seizmoverbose;
    
    % detail message
    if(verbose); disp('Reversing Correlogram(s)'); end
    
    % get header info
    [kuser0,kuser1,user0,user1,user2,user3,cmpinc,cmpaz,st,ev,...
        knetwk,kstnm,khole,kcmpnm,kt0,kt1,kt2,kt3,...
        b,e,a,f,t0,t1]=getheader(data,...
        'kuser0','kuser1','user0','user1','user2','user3','cmpinc',...
        'cmpaz','st','ev','knetwk','kstnm','khole','kcmpnm','kt0','kt1',...
        'kt2','kt3','b','e','a','f','t0','t1');
    
    % correlogram signature
    if(~all(strcmp('MASTER',kuser0) & strcmp('SLAVE',kuser1)))
        error('seizmo:reverse_correlations:notCorrelogram',...
            'Correlograms appear to have malformed headers!');
    end
    
    % detail message
    if(verbose); print_time_left(0,nrecs); end
    
    % fix header info
    data=changeheader(data,'b',-e,'e',-b,'st',ev,'ev',st,'user0',user1,...
        'user1',user0,'knetwk',kt0,'kstnm',kt1,'khole',kt2,'kcmpnm',kt3,...
        'kt0',knetwk,'kt1',kstnm,'kt2',khole,'kt3',kcmpnm,...
        'cmpinc',user2,'cmpaz',user3,'user2',cmpinc,'user3',cmpaz,...
        'a',t0,'f',t1,'t0',a,'t1',f);
    
    % reverse data
    data=reverse(data);
    
    % update delaz info
    oldstate=checkheader_state(true);
    data=checkheader(data,'all','ignore','old_delaz','fix');
    checkheader_state(oldstate);
    
    % rename
    digits=['%0' num2str(fix(log10(max([user0; user1])))+1) 'd'];
    name=strcat('CORR_-_MASTER_-_REC',num2str(user1,digits),'_-_',...
        knetwk,'.',kstnm,'.',khole,'.',kcmpnm,'_-_SLAVE_-_REC',...
        num2str(user0,digits),'_-_',kt0,'.',kt1,'.',kt2,'.',kt3);
    [data.name]=deal(name{:});
    
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
