function [data]=fix_cameroon(data)
%FIX_CAMEROON    Makes data corrections for the Cameroon Seismic Experiment
%
%    Usage:    data=fix_cameroon(data)
%
%    Description:
%     DATA=FIX_CAMEROON(DATA) makes simple corrections to the records in
%     DATA that correspond to problematic stations from the Camerooon
%     Seismic Experiment (deployed from 2005-2007, Network Code XB, Station
%     Names CM01-CM32).  These include amplitude and timing corrections.
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
%     % Clean up some Cameroon records acquired via RDSEED:
%     data=fix_cameroon(fix_rdseed_v48(readseizmo('*XB*CM*')));
%
%    See also: FIX_SOD_V222, FIX_DB2SAC_V48, FIX_RDSEED_V48

%     Version History:
%        Dec.  1, 2009 - initial version
%        Jan. 30, 2010 - reduce CHECKHEADER calls
%        Aug. 21, 2010 - nargchk fix, updated undef/nan handling, fixed
%                        warnings
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 21:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

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
    % get verbosity
    verbose=seizmoverbose;
    
    % get header info
    [b,kname]=getheader(data,'b utc','kname');
    b=cell2mat(b);

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
        if(verbose)
            warning('seizmo:fix_cameroon:CM02',...
                ['XB.CM02 amplitudes divided by 32!' ...
                '\nRecord(s):\n' sprintf('%d ',find(fixme))]);
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
        if(verbose)
            warning('seizmo:fix_cameroon:CM08',...
                ['XB.CM08 timing adjusted significantly (>600s)!' ...
                '\nRecord(s):\n' sprintf('%d ',find(fixme))]);
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
        if(verbose)
            warning('seizmo:fix_cameroon:CM14',...
                ['XB.CM14.0?.?HE amplitudes multiplied by 78.3!' ...
                '\nRecord(s):\n' sprintf('%d ',find(fixme))]);
        end
        data(fixme)=multiply(data(fixme),78.3);
    end

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
