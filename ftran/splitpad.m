function [data]=splitpad(data,pow2pad)
%SPLITPAD    Splits and zero pads SEIZMO data records
%
%    Usage:    data=splitpad(data)
%              data=splitpad(data,pow2pad)
%
%    Description:
%     DATA=SPLITPAD(DATA) splits records so the starting point is at 0 and
%     all points prior to that come at the end of the record(s).  The
%     record(s) are zeropadded to a power of 2.  Basically the layout is:
%     positive time data, zeros, negative time data.  This is to preserve
%     the phase information when transforming correlations of seismic noise
%     to the frequency domain.
%
%     DATA=SPLITPAD(DATA,POW2PAD) lets the power of 2 zero-padding be
%     adjusted using an integer POW2PAD according to the formula:
%                fftlength=2^(nextpow2(NPTS)+POW2PAD)
%     The default value is 0 and is a good choice for most problems.
%     Setting POW2PAD to <0 is not recommended and will truncate the time
%     series.  Setting POW2PAD >0 will increase the frequency resolution of
%     the spectrogram at the cost of increased computation time and array
%     size.
%
%    Notes:
%     - The record's must have a point at zero time.  An error is issued if
%       this is not the case.
%
%    Header changes: B, E, NPTS, DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Do NOT adjust the POW2PAD option to DFT when using SPLITPAD:
%     fdata=dft(splitpad(data));
%
%    See also: CUT, DFT

%     Version History:
%        Mar. 19, 2010 - initial version
%        Feb. 11, 2011 - mass nargchk fix
%        Feb. 15, 2011 - minor doc update
%        Aug. 26, 2011 - use checkheader more effectively, added some
%                        tolerance to the point at 0 requirement
%        Dec. 21, 2011 - doc update
%        May  29, 2012 - pow2pad=0 by default
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 10:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt to split & pad
try
    % check headers
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);
    
    % check pow2pad
    if(nargin<2 || isempty(pow2pad)); pow2pad=0; end
    if(~isnumeric(pow2pad) || ~isreal(pow2pad) ...
            || any(pow2pad~=fix(pow2pad)) ...
            || ~any(numel(pow2pad)==[1 nrecs]))
        error('seizmo:splitpad:badInput',...
            ['POW2PAD must be an integer or an array\n' ...
            'of integers (one per record)!']);
    end
    if(isscalar(pow2pad)); pow2pad=pow2pad(ones(nrecs,1),1); end
    pow2pad=pow2pad(:);
    
    % pull relevant timing info
    [b,delta,npts,ncmp]=getheader(data,'b','delta','npts','ncmp');
    nnpts=2.^(nextpow2n(npts)+pow2pad);
    e=b+(nnpts-1).*delta;
    
    % detail message
    if(verbose)
        disp('Splitting & Zero-padding Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % check has point at 0 (within some tolerance)
        idx=find(...
            abs(b(i)+(0:delta(i):delta(i)*(npts(i)-1)))<delta(i)/10,1);
        if(isempty(idx))
            error('seizmo:splitpad:noZeroPoint',...
                'Record %d does not have a point at zero!',i);
        end
        
        % get data class (for zeros)
        oclass=class(data(i).dep);
        
        % split & pad
        data(i).dep=[data(i).dep(idx:end,:);
            zeros(nnpts(i)-npts(i),ncmp(i),oclass); 
            data(i).dep(1:idx-1,:)];
        
        % dep*
        if(nnpts(i))
            depmen(i)=nanmean(data(i).dep(:));
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % update headers
    data=changeheader(data,'b',0,'e',e,'npts',nnpts,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);

    % rethrow error
    error(lasterror);
end

end
