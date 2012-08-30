function [data,fo,nyq]=iirfilter(data,varargin)
%IIRFILTER    Apply an IIR filter to SEIZMO records
%
%    Usage:    data=iirfilter(data,type,style,'option1',value1,...)
%              [data,fo,nyq]=iirfilter(...)
%              data=iirfilter(data,fo,nyq,'option1',value1,...)
%
%    Description:
%     DATA=IIRFILTER(DATA,TYPE,STYLE,'OPTION1',VALUE1,...) utilizes the set
%     of supplied parameters to build IIR filters to the desired specs and
%     implements them on the records in SEIZMO struct DATA.  The design
%     process includes the ability to automatically determine the best
%     fitting order (aka number of poles) when the filter is over designed.
%     Usually this is when the specs constain parameters like transition
%     bandwidth, stopband corners, attenuation, etc.  See the options below
%     for details/defaults.
%
%      Valid TYPE strings (NO DEFAULT):
%       'LO' - Low-pass filter
%       'HI' - High-pass filter
%       'BP' - Band-pass filter
%       'BS' - Band-stop filter
%
%      Valid STYLE strings (NO DEFAULT):
%       'BUTTER':  Filter is monotonic in passband and stopband.
%             Any specs beyond passband corner(s) and the filter order
%             triggers automatic filter order calculation to attempt to
%             implement them.  Order calculation will fit passband corners
%             exactly and will always have a transition bandwidth less than
%             or equal to that designed.  This may lead to inaccurate
%             stopband corner positioning.
%
%       'CHEBY1':  Filter is equiripple in the passband and monotonic in
%             the stopband.  Any specs beyond passband corners, maximum
%             passband attenuation, and the filter order triggers automatic
%             filter order calculation.  Order calculation will fit
%             passband corners exactly and will always have a transition
%             bandwidth less than or equal to that designed.  This may lead
%             to inaccurate stopband corner positioning.
%
%       'CHEBY2':  Filter is monotonic in the passband and equiripple in
%             the stopband.  Any specs beyond stopband corners, minimum
%             stopband attenuation, and the filter order triggers automatic
%             filter order calculation.  Order calculation will fit 
%             stopband corners exactly and will always have a transition 
%             bandwidth less than or equal to that designed.  This will 
%             give a wider passband than desired.
%
%       'ELLIP':  Filter is equiripple in passband and stopband.
%             Automatic filter order determination is only required for
%             specs defining the stopband corner positioning.  Order 
%             calculation will fit passband corners exactly and will always
%             have a transition bandwidth less than or equal to that 
%             designed.  This may lead to inaccurate stopband corner 
%             positioning.
%
%      AT LEAST ONE OF THE FOLLOWING OPTIONS IS REQUIRED:
%       IIRFILTER(...,'PCORNER',CORNERS,...)
%       IIRFILTER(...,'SCORNER',CORNERS,...)
%
%      THE REST ARE OPTIONAL:
%       IIRFILTER(...,'NPOLES',NUMBER_OF_POLES,...)
%       IIRFILTER(...,'PASSES',PASSES_OPTION,...)
%       IIRFILTER(...,'PATTEN',PASSBAND_ATTENUATION,...)
%       IIRFILTER(...,'SATTEN',STOPBAND_ATTENUATION,...)
%       IIRFILTER(...,'MIRROR',TRUE_OR_FALSE,...)
%       IIRFILTER(...,'TRANBW',FRACTION_OF_PASSBAND,...)
%
%      PCORNER: Array of passband corner(s). Must be 1 element for low/high
%       pass filters and 2 elements for bandpass/stop filters.  There is no
%       default value for PCORNER, but if SCORNER is defined then PCORNER
%       may be automatically defined using SCORNER and TRANBW.  PCORNER may
%       be set to an empty array if such behavior is desired.  Note that
%       PCORNER will be honored in auto-order determination for 'BUTTER',
%       'CHEBY1' & 'ELLIP' style filters.  'CHEBY2'-style filters only
%       honors the stopband corners (the passband will thus be wider than
%       requested).
%
%      SCORNER: Array of stopband corner(s). Must be 1 element for low/high
%       pass filters and 2 elements for bandpass/stop filters.  There is no
%       default value for SCORNER, but if PCORNER is defined then SCORNER
%       may be automatically defined using PCORNER and TRANBW.  SCORNER may
%       be set to an empty array if such behavior is desired.  Note that
%       SCORNER will be honored only for 'CHEBY2'-style filters (the
%       attenuation will be stronger than requested for the other styles).
%
%      NPOLES: Optional but SUGGESTED!  If undefined or 0, automatic filter
%       order calculations are made (utilizes the supplied/default corners
%       and attenuation parameters to define a suitable order).  Filter
%       order is analogous to the sharpness/steepness/slope of the
%       transition from the passband to the stopband.  Higher orders
%       provide better frequency resolution at the cost time resolution
%       (aka ringing).  Note that automatic order determination MAY ALTER
%       CORNER POSITIONING SIGNIFICANTLY (a warning is issued if so).
%
%      PASSES: Default is 1.  Accepts 1, 2, -1 or -2.  Option 1 is a simple
%       forward filter (single pass).  Option 2 filters forwards then
%       backwards to perform zero-phase filtering.  Option -1 is for
%       backwards filtering only and option -2 forward filters after
%       backwards filtering - aka a 2-pass in reverse order).  Two-pass
%       filtering is convenient for zero-phase filtering as a single-pass
%       IIR filter causes phase/group dispersion.  However, backwards
%       filtering makes the resulting signal acausal and should not be done
%       in studies where the onset of a phase needs to be preserved.
%
%      PATTEN: Default is 3dB.  Indicates the maximum permissible passband
%       amplitude loss/ripple.
%
%      SATTEN: Default is 30dB.  Indicates the minimum permissible stopband
%       loss or the maximum ripple (ie. nothing in the stopband is allowed
%       above this value).
%
%      MIRROR: Default is TRUE.  Indicates whether the filter is
%       implemented on record(s) with a mirror-flip leading the record
%       (TRUE) or if filtering is done in the normal fashion (FALSE).
%       Attaching a mirror-flip to the beginning of a record prior to
%       filtering is a method of reducing edge-effects by essentially
%       removing the jump from zero at the start of the record.  The cost
%       of this operation is that of filtering a record approximately twice
%       as long.
%
%      TRANBW: Default is 1/3.  Defines the width of the transition from
%       the stopband to the passband in terms of fraction of the passband
%       width.  Defining TRANBW will indicated the location of the
%       complimentary corner(s).  If both the passband and stopband
%       corner(s) are defined, then TRANBW is ignored.  Setting TRANBW
%       (even to the default value) will trigger automatic filter order
%       calculation.  To avoid this, set TRANBW to an empty array or do not
%       call the option.
%
%
%     [DATA,FO,NYQ]=IIRFILTER(...) outputs the filter object(s) in FO.
%     The filter objects can be visualized and checked using FVTOOL.
%     Remember that the filter is tied to the sampling frequency of the
%     record(s), so a dataset with records of varying samplerates will
%     produce a cell array of filter objects for FO.  The associated
%     Nyquist frequencies can be found in NYQ.  For example, NYQ(1) gives
%     the Nyquist frequency of filter object FO{1} and so on.
%
%     DATA=IIRFILTER(DATA,FO,NYQ,'OPTION1',VALUE1,...) allows reuse of
%     filters (see the Usage form directly above) or external design of
%     filters (requires a dfilt.df2tsos filter object!).  Only options
%     'MIRROR' and 'PASSES' are implemented (ie. all design options are
%     ignored in this case).  FO and NYQ must be the same size!
%
%    Notes: For narrow band filters, PARTICULARLY AT LOW Hz, consider 
%           decimating the signal (using SQUISH or SYNCRATES) prior to
%           filtering to aid in the accuracy of the filter design.  This is
%           because the analog to digital transformation of the filter
%           design will lose accuracy if the features of the filter are
%           too narrow.  Basically try to keep the filter width above a few
%           percent of the total bandwidth of the signal.  
%
%           Remember, length of the trace does not affect filter design,
%           but the sampling/Nyquist frequency does.
%
%    Header changes: DEPMIN, DEPMEN, DEPMAX
%
%    Examples:
%     % 4th order lowpass butter filter with a passband corner of 10s:
%     data=iirfilter(data,'low','butter','c',1/10,'np',4)
%
%     % Bandpass elliptic filter with passband corners of 15s and 12s using
%     % default transition bandwidth & attenuation to get the filter order:
%     data=iirfilter(data,'bp','ellip','c',[1/15 1/12])
%
%     % 2-pass 4th-order bandpass elliptic filter w/ 0.1dB passband ripple:
%     data=iirfilter(data,'bp','e','c',[1/20 1/15],'o',4,'p',2,'pr',0.1)
%
%     % Verify filter(s) with the FVTOOL (Signal Processing Toolbox):
%     [data,fo,nyq]=iirfilter(data,'bp','e','c',[1/20 1/15],'o',4);
%     fvtool(fo{:}) % note that this will plot all the filters together
%
%    See also: INTERPOLATE, SYNCRATES, SQUISH, STRETCH, DFT, IDFT,
%              IIRDESIGN, BUTTER, CHEBY1, CHEBY2, ELLIP, OMEGAGAUSSIAN

%     Version History:
%        Feb. 16, 2008 - initial version
%        Feb. 23, 2008 - minor improvements to corner adjustment
%        Feb. 29, 2008 - improved messages, comments, documentation
%        Mar.  6, 2008 - fixed bug to allow for all filter types, fixed bug
%                        in corner adjustment, auto-order docs/warning
%        May   5, 2008 - minor doc update
%        Oct.  7, 2008 - .x to .dep update
%        Nov. 16, 2008 - update warnings/errors for name change
%        Nov. 25, 2008 - improved use of checking, update for name changes,
%                        brought records-matrix-records out of subfunction
%        Apr. 23, 2009 - fixes for octave compatibility
%        Jun. 25, 2009 - update for name changes
%        Oct. 12, 2009 - minor doc update
%        Dec.  7, 2009 - minor doc update
%        Feb.  2, 2010 - complete code rehash, filter design now in
%                        IIRDESIGN (highly updated), filters record by
%                        record, seizmoverbose support, proper SEIZMO
%                        handling, versioninfo caching, better checks,
%                        added version history, complete doc reformat
%        Sep. 20, 2010 - passes arg now uses -1,-2 rather than 3,4
%        Feb. 11, 2011 - mass nargchk fix, dropped versioninfo caching
%        Jan. 28, 2012 - doc update, better checkheader usage
%        Feb.  5, 2012 - added omegagaussian to see also
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  5, 2012 at 15:05 GMT

% todo:

% check number of inputs
error(nargchk(2,inf,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt filtering
try
    % check headers
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);
    
    % get header info
    delta=getheader(data,'delta');
    
    % get nyquist frequencies
    [dt,idx,idx]=unique(delta);
    nyq=1./(2*dt);
    if(numel(nyq)>1)
        warning('seizmo:iirfilter:multipleSampleRates',...
            ['Dataset has multiple sample rates!\n' ...
            'Filters may vary significantly with sample rate!']);
    end
    
    % get filter design
    [fo,passes,mirror]=iirdesign(nyq,varargin{:});
    
    % detail message
    if(verbose)
        disp('Filtering Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end
        
        % get data class, convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % filter
        if(passes==1)
            data(i).dep=impfilt(data(i).dep,fo{idx(i)},mirror,true);
        elseif(passes==2)
            data(i).dep=impfilt(data(i).dep,fo{idx(i)},mirror,true);
            data(i).dep=impfilt(data(i).dep,fo{idx(i)},mirror,false);
        elseif(passes==-1)
            data(i).dep=impfilt(data(i).dep,fo{idx(i)},mirror,false);
        elseif(passes==-2)
            data(i).dep=impfilt(data(i).dep,fo{idx(i)},mirror,false);
            data(i).dep=impfilt(data(i).dep,fo{idx(i)},mirror,true);
        end
        
        % change class back
        data(i).dep=oclass(data(i).dep);
        
        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % update header
    data=changeheader(data,...
        'depmax',depmax,'depmin',depmin,'depmen',depmen);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end

function [x]=impfilt(x,fo,mirror,forward)
%IMPFILT   Implements filter
%   Takes a mirror option which does pseudo-IC to limit edge effects.
%   Forward logical allows for forward or reverse filtering.

% reverse
if(~forward); x=x(end:-1:1,:); end

% mirror-flip
if(mirror)
    % prepend a mirror-flip of the series to limit edge effects
    x=filter(fo,[2*x(ones(end-1,1),:)-x(end:-1:2,:); x],1);
    x=x(ceil(end/2):end,:);
else
    % normal filtering
    x=filter(fo,x,1);
end

% reverse
if(~forward); x=x(end:-1:1,:); end

end
