function [data]=fix_db2sac_v48(data,knetwk)
%FIX_DB2SAC_V48    Cleans up headers of SAC files created with DB2SAC
%
%    Usage:    data=fix_db2sac_v48(data,knetwk)
%
%    Description:
%     DATA=FIX_DB2SAC_V48(DATA,KNETWK) fixes the headers of SAC files
%     exported using DB2SAC from ANTELOPE version 4.8.  The fixes are:
%      1. sets IDEP field to 'IUNKN' (unknown units) - for db2sac -counts
%      2. sets IZTYPE field to 'IB'
%      3. sets KNETWK field to KNETWK (note this is the 2nd input argument)
%      4. splits KCMPNM field to set KHOLE and KCMPNM
%      5. sets empty KHOLE to '__'
%      6. sets LOVROK to TRUE (allow overwrite)
%      7. set USER7, USER8, NORID, NEVID to undefined
%
%    Notes:
%     - SAC records created with DB2SAC may be off by 1 millisecond
%       compared to records created with TREXCERPT.  No fix is available
%       for this discrepancy.  A check appeared to show that TREXCERPT is
%       accurate.
%
%    Header changes: IDEP, IZTYPE, KNETWK, KCMPNM, KHOLE, LOVROK,
%                    USER7, USER8, NORID, NEVID
%
%    Examples:
%     % Read, clean up, and overwrite some DB2SAC-made SAC files:
%     w(fix_db2sac_v48(r('*')))
%
%    See also: FIX_SOD_V222, FIX_RDSEED_V48, FIX_TREXCERPT_V48,
%              FIX_CAMEROON, FIXDELTA

%     Version History:
%        Dec.  1, 2009 - initial version
%        Dec.  2, 2009 - no data requirement now
%        Dec.  4, 2009 - USER7, USER8, NORID, NEVID unset
%        Jan. 30, 2010 - reduced calls to seizmocheck
%        Mar. 24, 2010 - drop fixdelta call
%        Feb. 11, 2011 - mass nargchk fix
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure & header
data=checkheader(data);

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt fixes
try
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
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
