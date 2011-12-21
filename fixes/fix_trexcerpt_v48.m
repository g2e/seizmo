function [data]=fix_trexcerpt_v48(data,knetwk)
%FIX_TREXCERPT_V48    Cleans up headers of SAC files created with TREXCERPT
%
%    Usage:    data=fix_trexcerpt_v48(data,knetwk)
%
%    Description:
%     DATA=FIX_TREXCERPT_V48(DATA,KNETWK) fixes the headers of SAC files
%     exported using TREXCERPT from ANTELOPE version 4.8.  The fixes are:
%     1. sets IDEP field to 'IUNKN' (unknown units)
%     2. sets IZTYPE field to 'IB'
%     3. sets KNETWK field to KNETWK (note this is the 2nd input argument)
%     4. splits KCMPNM field to set KHOLE and KCMPNM
%     5. sets empty KHOLE to '__'
%     6. sets LOVROK to TRUE (allow overwrite)
%     7. set NORID, NEVID to undefined
%     8. fixes E field
%
%    Notes:
%     - DEPMIN, DEPMEN, DEPMAX are set if the data is also read in
%
%    Header changes: DEP, IZTYPE, KNETWK, KCMPNM, KHOLE, LOVROK,
%                    NORID, NEVID, E, DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     % Read, clean up, and overwrite some TREXCERPT-made SAC files:
%     w(fix_trexcerpt_v48(r('*')))
%
%    See also: FIX_SOD_V222, FIX_RDSEED_V48, FIX_CAMEROON, FIX_DB2SAC_V48,
%              FIXDELTA

%     Version History:
%        Dec.  4, 2009 - initial version
%        Jan. 30, 2010 - reduced calls to seizmocheck
%        Mar. 24, 2010 - drop fixdelta call
%        Feb. 11, 2011 - mass nargchk fix
%        Dec.  1, 2011 - make checkheader ignore undefined iztype
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  1, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure & header
data=checkheader(data,'invalid_iztype','ignore');

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
        'knetwk',knetwk,'kcmpnm',kcmpnm,'khole',khole,...
        'nevid',nan,'norid',nan);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
