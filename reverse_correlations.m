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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 12, 2010 at 21:45 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

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
end

% attempt reversal
try
    % detail message
    if(seizmoverbose); disp('Reversing Correlogram(s)'); end
    
    % get header info
    [kuser0,kuser1,user0,user1,st,ev,knetwk,kstnm,khole,kcmpnm,...
        kt0,kt1,kt2,kt3,b,e]=getheader(data,'kuser0','kuser1','user0',...
        'user1','st','ev','knetwk','kstnm','khole','kcmpnm','kt0','kt1',...
        'kt2','kt3','b','e');
    
    % correlogram signature
    if(~all(strcmp('MASTER',kuser0) & strcmp('SLAVE',kuser1)))
        warning('seizmo:reverse_correlations:notCorrelogram',...
            '');
    end
    
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
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
