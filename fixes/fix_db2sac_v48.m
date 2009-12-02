function [data]=fix_db2sac_v48(data,knetwk)
%FIX_DB2SAC_V48    Cleans up header info for SAC files created with DB2SAC
%
%    Usage:    data=fix_db2sac_v48(data,knetwk)
%
%    Description: DATA=FIX_DB2SAC_V48(DATA,KNETWK) fixes the headers of SAC
%     files exported using DB2SAC from ANTELOPE version 4.8.  The fixes
%     are:
%     1. adjusts DELTA field slightly (single to double precision issue)
%     2. sets IDEP field to 'IUNKN' (unknown units) - for db2sac -counts
%     3. sets IZTYPE field to 'IB'
%     4. sets KNETWK field to KNETWK (note this is the 2nd input argument)
%     5. splits KCMPNM field to set KHOLE and KCMPNM
%     6. sets empty KHOLE to '__'
%     7. sets LOVROK to TRUE (allow overwrite)
%
%    Notes:
%
%    Header changes: DELTA, IDEP, IZTYPE, KNETWK, KCMPNM, KHOLE, LOVROK
%
%    Examples:
%     Read, clean up, and overwrite some DB2SAC-made SAC files:
%      w(fix_db2sac_v48(r('*')))
%
%    See also: FIX_SOD_V222, FIX_RDSEED_V48, FIX_CAMEROON

%     Version History:
%        Dec.  1, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  1, 2009 at 22:45 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=get_checkheader_state;
    set_checkheader_state(false);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt rest
try
    % fix delta
    data=fixdelta(data);
    
    % fix kname
    kcmpnm=getheader(data,'kcmpnm');
    khole=kcmpnm;
    for i=1:numel(data)
        words=getwords(kcmpnm{i},'_');
        kcmpnm(i)=words(1);
        if(numel(words)>1)
            khole(i)=words(2);
        else
            % force '__' for empty khole
            khole(i)={'__'};
        end
    end
    
    % update header
    data=changeheader(data,'idep','iunkn','iztype','ib','lovrok',true,...
        'knetwk',knetwk,'kcmpnm',kcmpnm,'khole',khole);

    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
