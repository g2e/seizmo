function [data]=reverse_correlations(data)
%REVERSE_CORRELATIONS    Reverses correlations (switch master & slave)
%
%    Usage:    data=reverse_correlations(data)
%
%    Description: DATA=REVERSE_CORRELATIONS(DATA) does a time reversal of
%     correlograms in SEIZMO struct DATA and switches the master and slave
%     header field info.  This is only for correlations generated with
%     CORRELATE.
%
%    Notes:
%     - Checks that all records have the fields KUSER0 & KUSER1 set to
%       'MASTER' & 'SLAVE'.  Warns if they are not (an indication that the
%       records are not from CORRELATE).
%
%    Header changes:
%     Switches B & E :
%      B = -E
%      E = -B
%     Master & Slave Field Switches :
%      STLA   <-> EVLA
%      STLO   <-> EVLO
%      STEL   <-> EVEL
%      STDP   <-> EVDP
%      KNETWK <-> KT0
%      KSTNM  <-> KT1
%      KHOLE  <-> KT2
%      kCMPNM <-> KT3
%      USER0  <-> USER1
%     Updates:
%      GCARC, AZ, BAZ, DIST
%
%    Examples:
%     Get all cross correlation pairs (only computing half of them):
%      corr1=correlate(data);
%      corr2=reverse_correlations(corr1);
%
%    See also: CORRELATE, REVERSE

%     Version History:
%        Apr. 12, 2010 - initial version
%        Apr. 15, 2010 - resets name field, includes scrollbar if verbose
%        Feb. 11, 2011 - mass nargchk fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data (dep)
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt reversal
try
    % number of records
    nrecs=numel(data);
    
    % check verbosity
    verbose=seizmoverbose;
    
    % detail message
    if(verbose); disp('Reversing Correlogram(s)'); end
    
    % get header info
    [kuser0,kuser1,user0,user1,st,ev,knetwk,kstnm,khole,kcmpnm,...
        kt0,kt1,kt2,kt3,b,e]=getheader(data,'kuser0','kuser1','user0',...
        'user1','st','ev','knetwk','kstnm','khole','kcmpnm','kt0','kt1',...
        'kt2','kt3','b','e');
    
    % correlogram signature
    if(~all(strcmp('MASTER',kuser0) & strcmp('SLAVE',kuser1)))
        warning('seizmo:reverse_correlations:notCorrelogram',...
            'Some records are probably not correlograms!');
    end
    
    % detail message
    if(verbose); print_time_left(0,nrecs); end
    
    % fix header info
    data=changeheader(data,'b',-e,'e',-b,'st',ev,'ev',st,'user0',user1,...
        'user1',user0,'knetwk',kt0,'kstnm',kt1,'khole',kt2,'kcmpnm',kt3,...
        'kt0',knetwk,'kt1',kstnm,'kt2',khole,'kt3',kcmpnm);
    
    % reverse data
    data=reverse(data);
    
    % update delaz info
    checkheader_state(true);
    data=checkheader(data,'all','ignore','old_delaz','fix');
    checkheader_state(false);
    
    % rename
    digits=['%0' num2str(fix(log10(max([user0; user1])))+1) 'd'];
    name=strcat('CORR_-_MASTER_-_REC',num2str(user1,digits),'_-_',knetwk,'.',kstnm,'.',khole,'.',kcmpnm,'_-_SLAVE_-_REC',num2str(user0,digits),'_-_',kt0,'.',kt1,'.',kt2,'.',kt3);
    [data.name]=deal(name{:});
    
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end
    
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
