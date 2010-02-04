function [data]=setasiday(data)
%SETASIDAY    Sets the reference times of records to the start of the day
%
%    Usage:    data=setasiday(data)
%
%    Description: DATA=SETASIDAY(DATA) sets the refence times of records in
%     SEIZMO struct DATA to the beginning of the day that the reference
%     time is currently in.  For example, if the reference time of a record
%     is 2006:145:15:33:10.885 then SETASIDAY will shift it to
%     2006:145:00:00:00.000.  All timing fields are correspondingly shifted
%     to preserve the timing of the data.  The IZTYPE field is set to IDAY.
%
%    Notes:
%     - This is mainly useful for ambient noise studies for synchronizing
%       records to the beginning of the day that the data is recorded on.
%
%    Header changes: NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC,
%                    A, B, E, F, O, Tn, IZTYPE
%
%    Examples:
%     Comparison of header info before and after SETASIDAY:
%
%      listheader(data0(1),'b','nz')
%        FILE: STEP_1-DB2SAC_LHZ_RAW/2005009135451.81.CM07.LHZ_01 - 1
%       ---------------------------
%                      B = 0
%                 NZYEAR = 2005
%                 NZJDAY = 9
%                 NZHOUR = 13
%                  NZMIN = 54
%                  NZSEC = 51
%                 NZMSEC = 814
%
%      listheader(setasiday(data0(1)),'b','nz')
%        FILE: STEP_1-DB2SAC_LHZ_RAW/2005009135451.81.CM07.LHZ_01 - 1
%       ---------------------------
%                       B = 50091.814
%                  NZYEAR = 2005
%                  NZJDAY = 9
%                  NZHOUR = 0
%                   NZMIN = 0
%                   NZSEC = 0
%                  NZMSEC = 0
%
%    See also: TIMESHIFT, SYNCHRONIZE

%     Version History:
%        Dec.  2, 2009 - initial version
%        Feb.  3, 2010 - versioninfo caching, seizmoverbose support
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  3, 2010 at 18:10 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
[h,idx]=versioninfo(data);

% get undefined values
undef=getsubfield(h,'undef','ntype').';
undef=undef(idx);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
oldversioninfocache=versioninfo_cache(true);

% attempt rest
try
    % verbosity
    verbose=seizmoverbose;
    
    % detail message
    if(verbose)
        disp('Shifting Reference Times of Record(s) to a Day Boundary');
    end
    
    % get reference time
    [nz]=getheader(data,'nz');

    % find those with reftime undefined
    bad=sum(nz==undef(:,ones(1,6)),2)>0;
    if(any(bad))
        error('seizmo:setiday:NZUndefined',...
            ['Records:\n' sprintf('%d ',find(bad)) ...
            '\nOne or more NZ fields are undefined!']);
    end

    % now shift the records
    data=timeshift(data,...
        nz(:,3)*3600+nz(:,4)*60+nz(:,5)+nz(:,6)/1000,'iday');

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
    
    % rethrow error
    error(lasterror)
end

end
