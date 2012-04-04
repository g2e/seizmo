function [data]=fix_sod_v222(data)
%FIX_SOD_V222    Cleans up headers of SAC files created with SOD
%
%    Usage:    data=fix_sod_v222(data)
%
%    Description:
%     DATA=FIX_SOD_V222(DATA) fixes the headers of SAC files
%     exported using SOD v2.2.2.  The fixes are:
%     1. correctly set E field
%     2. setting IZTYPE to IO (if o field is zero)
%     3. fix KCMPNM to not include the stream code
%     4. setting EVEL to 0 (if EVLA/EVLO/EVDP defined)
%     5. setting empty KHOLE to '__'
%
%    Notes:
%
%    Header changes: E, IZTYPE, KCMPNM, KHOLE, EVEL
%
%    Examples:
%     % Read, clean up, and overwrite some SOD-made SAC files:
%     w(fix_sod_v222(r('*')))
%
%    See also: FIX_RDSEED_V48, FIX_DB2SAC_V48, FIX_TREXCERPT_V48,
%              FIX_CAMEROON, FIXDELTA

%     Version History:
%        Nov. 23, 2009 - initial version
%        Dec.  1, 2009 - minor doc update
%        Dec.  2, 2009 - no data requirement now
%        Dec.  4, 2009 - minor doc update
%        Jan. 30, 2010 - fixes for checking state functions
%        Mar. 24, 2010 - drop fixdelta call
%        Aug. 21, 2010 - nargchk fix, updated undef/nan handling
%        Jan. 31, 2011 - minor doc fixes
%        Mar. 19, 2011 - stream code fix was wrong (sigh)
%        Apr.  3, 2012 - use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 22:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % do iztype fix before checkheader to avoid warning
    zo=getheader(data,'o')==0;
    if(any(zo))
        data(zo)=changeheader(data(zo),'iztype','io');
    end
    
    % check header
    data=checkheader(data);
    
    % fix evel
    evdef=~isnan(getheader(data,'ev'));
    fixel=evdef(:,1) & evdef(:,2) & ~evdef(:,3) & evdef(:,4);
    if(any(fixel))
        data(fixel)=changeheader(data(fixel),'evel',0);
    end
    
    % drop stream code from kcmpnm
    % 'fix' empty khole to '__' (for external scripts)
    [khole,kcmpnm]=getheader(data,'khole','kcmpnm');
    khole(strcmp(khole,''))={'__'};
    for i=1:numel(kcmpnm); kcmpnm{i}=kcmpnm{i}(end-2:end); end
    data=changeheader(data,'kcmpnm',kcmpnm,'khole',khole);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
