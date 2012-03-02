function [data]=fix_rdseed_v48(data)
%FIX_RDSEED_V48    Cleans up headers of SAC files created with RDSEED
%
%    Usage:    data=fix_rdseed_v48(data)
%
%    Description:
%     DATA=FIX_RDSEED_V48(DATA) fixes the headers of SAC files
%     exported using RDSEED v4.8.  The fixes are:
%     1. synchronizing reference time to the origin time (if o defined)
%     2. setting IZTYPE to IO (if o defined)
%     3. setting LOVROK to TRUE (allow overwrite)
%     4. setting EVEL to 0 (if EVLA/EVLO/EVDP defined)
%     5. rounding MAG to the nearest hundredth (if defined)
%     6. setting empty KHOLE to '__'
%     7. changing evdp from km->m if evdp<1000
%
%    Notes:
%     - DEPMIN, DEPMEN, DEPMAX are set if the data is also read in
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX, IZTYPE, EVDP, EVEL, KHOLE,
%     LOVROK, MAG, NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC,
%     A, F, O, B, E, Tn
%
%    Examples:
%     % Read, clean up, and overwrite some RDSEED-made SAC files:
%     w(fix_rdseed_v48(r('*')))
%
%    See also: FIX_SOD_V222, FIX_DB2SAC_V48, FIX_TREXCERPT_V48,
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
%        Mar.  1, 2012 - use seizmocheck, km->m depth fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2012 at 22:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt rest
try
    % get origin info
    [o,ev,mag]=getheader(data,'o','ev','mag');
    
    % who's defined
    odef=~isnan(o);
    evdef=~isnan(ev);
    fixel=evdef(:,1) & evdef(:,2) & ~evdef(:,3) & evdef(:,4);
    fixdp=ev(:,4)<1000;
    magdef=~isnan(mag);
    
    % fix origin
    if(any(odef))
        data(odef)=timeshift(data(odef),-o,'io');
        data(odef)=changeheader(data(odef),'o',0);
    end
    if(any(fixel))
        data(fixel)=changeheader(data(fixel),'evel',0);
    end
    if(any(fixdp))
        data(fixdp)=changeheader(data(fixdp),'evdp',ev(fixdp,4)*1000);
    end
    if(any(magdef))
        data(magdef)=changeheader(data(magdef),'mag',round(mag*100)/100);
    end
    
    % lovrok to true
    % 'fix' empty khole to '__' (for external scripts)
    khole=getheader(data,'khole');
    khole(strcmp(khole,''))={'__'};
    data=changeheader(data,'lovrok',true,'khole',khole);

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
