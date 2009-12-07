function [data]=fix_cameroon(data)
%FIX_CAMEROON    Makes data corrections for the Cameroon Seismic Experiment
%
%    Usage:    data=fix_cameroon(data)
%
%    Description: DATA=FIX_CAMEROON(DATA) makes simple corrections to the
%     records in DATA that correspond to problematic stations from the
%     Camerooon Seismic Experiment (deployed from 2005-2007, Network Code
%     XB, Station Names CM01-CM32).  These include amplitude and timing
%     corrections.
%
%    Notes:
%     - Look at the comments within this script for details on corrections.
%     - CM08 timing correction resets IZTYPE to 'iunkn' and all timing
%       fields except 'b' and 'e' (ie 'a' 'f' 'o' 'tn') are shifted too.
%       So you will likely need to re-synchronize the records.
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX, IZTYPE, NZYEAR, NZJDAY,
%     NZHOUR, NZMIN, NZSEC, NZMSEC, A, F, O, Tn
%
%    Examples:
%     Clean up some Cameroon records acquired via RDSEED:
%      data=fix_cameroon(fix_rdseed_v48(readseizmo('*XB*CM*')));
%
%    See also: FIX_SOD_V222, FIX_DB2SAC_V48, FIX_RDSEED_V48

%     Version History:
%        Dec.  1, 2009 - initial version
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
sundef=getsubfield(h,'undef','stype').';
undef=undef(idx);
sundef=sundef(idx);

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% get verbosity
v=seizmoverbose;

% attempt rest
try
    % get header info
    [b,kname]=getheader(data,'b utc','kname');

    % unwrap cells
    b=cell2mat(b);
    
    % check that b, kname are defined
    bad=sum(b==undef(:,ones(1,5)),2)>0 ...
        | sum(strcmpi(kname,sundef(:,ones(1,4))),2)>0;
    if(any(bad))
        v=seizmoverbose;
        if(v)
            warning('seizmo:fix_cameroon:badHeader',...
                ['Records:\n' sprintf('%d ',find(bad)) ...
                '\nOne or more of the following fields are undefined:\n'...
                'B NZ KNAME\nThe records will likely not be corrected!']);
        end
        b(bad,:)=nan;
    end

    % fix CM02 - AMPLITUDE OF ALL CHANNELS ARE OFF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % divide by 32
    %
    % applies to all data recorded after May 15, 2006
    %
    % Caused by a corruption of the DAS memory which caused a multiplier of
    % 32 to be added to the DAS parameters (can be seen in the logfiles).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fixme=strcmpi(kname(:,1),'XB') & strcmpi(kname(:,2),'CM02') ...
        & ((b(:,1)==2006 & b(:,2)>135) | b(:,1)==2007);
    if(any(fixme))
        if(v)
            warning('seizmo:fix_cameroon:CM02',...
                ['Records:\n' sprintf('%d ',find(fixme)) ...
                '\nXB.CM02 amplitudes divided by 32!']);
        end
        data(fixme)=divide(data(fixme),32);
    end

    % fix CM08 - ARRIVALS ARE EARLY BY 10+ MINUTES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TIMING_GOOD = TIMING_BAD + CORRECTION
    %
    % where correction is:
    %
    % from 2006 Jan 1st 00:00 GMT to Jun 19 16:00 GMT
    % 689.00-seconds_since_start_of_2006*0.0000012594763
    %
    % from 2006 Jun 19 16:01 GMT to July 4 00:00 GMT
    % 799.55-seconds_since_start_of_2006*0.0000012594763
    %
    % 0.0000012594763 ~ drift of roughly 1sec/10days
    %
    % Empirical correction was found from logfile entries of clock status
    % where the clock error in integer seconds was given.  A best fit line
    % through all the errors gives the correction slope.  Intercept was
    % chosen to best match how clock error integer was reported (truncation
    % was assumed based on how RT-130 systems handle leap second data) so
    % the best fit line "rests" on top of the step-like integer dataset.
    %
    % Honestly this should be redone/checked more thoroughly.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fixme=strcmpi(kname(:,1),'XB') & strcmpi(kname(:,2),'CM08') ...
        & b(:,1)==2006 & b(:,2)<185;
    if(any(fixme))
        if(v)
            warning('seizmo:fix_cameroon:CM08',...
                ['Records:\n' sprintf('%d ',find(fixme)) ...
                '\nXB.CM08 timing adjusted significantly (>600s)!']);
        end
        % get time since 2006.001 00:00:00.000
        % multiply by 0.0000012594763 to get drift
        % add dc-shift of 689 or 799.55
        z=[2006 1 0 0 0];
        zoffset=timediff(z,b(fixme,:));
        zoffset=zoffset*0.0000012594763;
        s689=b(fixme,2)<170 | (b(fixme,2)==170 & b(fixme,3)<16);
        s689=s689*689; s689(s689==0)=799.55;
        cm08corr=s689-zoffset;
        % now shift reference time and relative fields
        data(fixme)=timeshift(data(fixme),-cm08corr);
        % unshift b & e fields (but not the rest)
        data(fixme)=timeshift(data(fixme),cm08corr,[],[],'user','b','e');
    end

    % fix CM14 - AMPLITUDE OF THE EW CHANNEL IS LOW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multiply by 78.3
    %
    % applies for entire time station is running
    %
    % Empirical correction was derived from comparing the calibration pulse
    % amplitudes of the east and north channels.  High frequencies are
    % corrupted by electronic noise - filtering out high frequencies is
    % necessary.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fixme=strcmpi(kname(:,1),'XB') & strcmpi(kname(:,2),'CM14') ...
        & (strcmpi(kname(:,4),'BHE') | strcmpi(kname(:,4),'LHE')) ...
        & (b(:,1)==2006 | b(:,1)==2007);
    if(any(fixme))
        if(v)
            warning('seizmo:fix_cameroon:CM14',...
                ['Records:\n' sprintf('%d ',find(fixme)) ...
                '\nXB.CM14.0?.?HE amplitudes multiplied by 78.3!']);
        end
        data(fixme)=multiply(data(fixme),78.3);
    end

    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
