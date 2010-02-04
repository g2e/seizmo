function [data]=fix_db2sac_v48(data,knetwk)
%FIX_DB2SAC_V48    Cleans up headers of SAC files created with DB2SAC
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
%     8. set USER7, USER8, NORID, NEVID to undefined
%
%    Notes:
%     - SAC records created with DB2SAC may be off by 1 millisecond
%       compared to records created with TREXCERPT.  No fix is available
%       for this discrepancy.  A check appeared to show that TREXCERPT is
%       accurate.
%
%    Header changes: DELTA, IDEP, IZTYPE, KNETWK, KCMPNM, KHOLE, LOVROK,
%                    USER7, USER8, NORID, NEVID
%
%    Examples:
%     Read, clean up, and overwrite some DB2SAC-made SAC files:
%      w(fix_db2sac_v48(r('*')))
%
%    See also: FIX_SOD_V222, FIX_RDSEED_V48, FIX_TREXCERPT_V48,
%              FIX_CAMEROON

%     Version History:
%        Dec.  1, 2009 - initial version
%        Dec.  2, 2009 - no data requirement now
%        Dec.  4, 2009 - USER7, USER8, NORID, NEVID unset
%        Jan. 30, 2010 - reduced calls to seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 30, 2010 at 20:45 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure & header
data=checkheader(data);

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% attempt fixes
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
        'user7',nan,'user8',nan,'knetwk',knetwk,'kcmpnm',kcmpnm,...
        'khole',khole,'nevid',nan,'norid',nan);

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
