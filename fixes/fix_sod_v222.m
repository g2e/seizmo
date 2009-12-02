function [data]=fix_sod_v222(data)
%FIX_SOD_V222    Cleans up header info for SAC files created with SOD
%
%    Usage:    data=fix_sod_v222(data)
%
%    Description: DATA=FIX_SOD_V222(DATA) fixes the headers of SAC files
%     exported using SOD v2.2.2.  The fixes are:
%     1. correctly set E field
%     2. setting IZTYPE to IO (if o field is zero)
%     3. fix KCMPNM to not include the stream code
%     4. setting EVEL to 0 (if EVLA/EVLO/EVDP defined)
%     5. adjusting DELTA slightly (single to double precision issue)
%     6. setting empty KHOLE to '__'
%
%    Notes:
%
%    Header changes: E, DELTA, IZTYPE, KCMPNM, KHOLE, EVEL
%
%    Examples:
%     Read, clean up, and overwrite some SOD-made SAC files:
%      w(fix_sod_v222(r('*)))
%
%    See also: FIX_RDSEED_V48, FIX_DB2SAC_V48, FIX_CAMEROON

%     Version History:
%        Nov. 23, 2009 - initial version
%        Dec.  1, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  1, 2009 at 22:45 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
[h,idx]=versioninfo(data,'dep');

% get undefined values
undef=getsubfield(h,'undef','ntype').';
undef=undef(idx);

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % do iztype fix before checkheader to avoid warning
    o=getheader(data,'o');
    zo=o==0;
    if(any(zo))
        data(zo)=changeheader(data(zo),'iztype','io');
    end
    
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
    
    % fix evel
    ev=getheader(data,'ev');
    evdef=ev~=undef(:,[1 1 1 1]);
    fixel=evdef(:,1) & evdef(:,2) & ~evdef(:,3) & evdef(:,4);
    if(any(fixel))
        data(fixel)=changeheader(data(fixel),'evel',0);
    end
    
    % drop stream code from kcmpnm
    % 'fix' empty khole to '__' (for external scripts)
    [khole,kcmpnm]=getheader(data,'khole','kcmpnm');
    khole(strcmp(khole,''))={'__'};
    kcmpnm=char(strnlen(kcmpnm,3));
    data=changeheader(data,'kcmpnm',kcmpnm,'khole',khole);

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
