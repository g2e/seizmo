function [data]=powerspectraldensity(data,method,pow2pad)
%POWERSPECTRALDENSITY    Returns power spectral density of SEIZMO records
%
%    Usage:    data=powerspectraldensity(data)
%              data=powerspectraldensity(data,method)
%              data=powerspectraldensity(data,method,pow2pad)
%
%    Description:
%     DATA=POWERSPECTRALDENSITY(DATA) computes the power spectral density
%     for each record in SEIZMO struct DATA using the fft approach.  Output
%     records are in decibels (dBs).
%
%     DATA=POWERSPECTRALDENSITY(DATA,METHOD) changes the method used to
%     compute the power spectra.  Currently there is only 'fft'.  Other
%     methods such as multitaper, music, and welch may be added later.
%     Giving [] will use the default of 'fft'.
%
%     DATA=POWERSPECTRALDENSITY(DATA,METHOD,POW2PAD) lets the power of 2
%     zero-padding be adjusted by an integer POW2PAD in the formula:
%                fftlength=2^(nextpow2(NPTS)+POW2PAD)
%     The default value is 0 and provides good frequency resolution.
%     Using a value of >0 will improve the frequency resolution.
%
%    Notes:
%     - Setting POW2PAD<0 will wrap the data (summing the segments
%       together) to preserve the power spectra.  This differs from DFT
%       which truncates the data.
%
%    Header Changes: B, SB, E, DELTA, SDELTA, NPTS, NSPTS, IFTYPE
%                    DEPMEN, DEPMIN, DEPMAX
%     B, E, DELTA, and NPTS are changed to the beginning frequency, nyquist
%     frequency, frequency interval, and number of unique frequency points
%     (we drop the negative frequencies here).  The original values 
%     of B, DELTA, and NPTS are saved in the header as SB, SDELTA, and 
%     NSPTS and are restored when the IDFT command is performed.  Using
%     POW2PAD less than 0 will change the value stored in NSPTS.  IFTYPE is
%     altered to indicate an XY filetype.
%
%    Examples:
%     % These should be the same:
%     plot2(powerspectradensity(data),keeppw(dft(data)),'xscale','log');
%
%     % See the affect of the POW2PAD option:
%     plot2(powerspectraldensity(data([1 1 1]),[],-1:1),'xscale','log');
%
%    See also: DFT, IDFT, KEEPPW

%     Version History:
%        Dec.  1, 2011 - initial version
%        Dec.  5, 2011 - more comments
%        Dec. 14, 2011 - non-object version (yay!)
%        Dec. 21, 2011 - wrapping for pow2pad<0, apply delta so spectra is
%                        independent of sample rate, pow2pad<0 wrapping fix
%                        (required division by original npts), examples, no
%                        need to turn off checkheader
%        Feb.  6, 2012 - output is in dBs now
%        May  29, 2012 - pow2pad=0 by default
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  29, 2012 at 13:30 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% get psd
try
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % get number of records
    nrecs=numel(data);
    
    % method
    valid={'fft'};
    if(nargin<2 || isempty(method)); method='fft'; end
    if(ischar(method)); method=cellstr(method); end
    if(~iscellstr(method) || any(~ismember(lower(method),valid)) ...
            || ~any(numel(method)==[1 nrecs]))
        error('seizmo:powerspectraldensity:badInput',...
            ['METHOD must be one of the following:\n' ...
            sprintf('%s ',valid{:})]);
    end
    if(isscalar(method)); method(1:nrecs,1)=method; end
    method=method(:);
    
    % pow2pad
    if(nargin<3 || isempty(pow2pad)); pow2pad=0; end
    if(~isnumeric(pow2pad) || ~isreal(pow2pad) ...
            || any(pow2pad~=fix(pow2pad)) ...
            || ~any(numel(pow2pad)==[1 nrecs]))
        error('seizmo:powerspectraldensity:badInput',...
            ['POW2PAD must be an integer or an array\n' ...
            'of integers (one per record)!']);
    end
    if(isscalar(pow2pad)); pow2pad(1:nrecs,1)=pow2pad; end
    pow2pad=pow2pad(:);
    
    % get sample interval & number of points
    [delta,npts,b]=getheader(data,'delta','npts','b');
    
    % get frequency info
    nspts=2.^(nextpow2n(npts)+pow2pad);
    sb=0; se=1./(2*delta); sdelta=1./(delta.*nspts);
    
    % detail message
    if(verbose)
        disp('Getting Power Spectral Density of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over each record
    depmen=nan(nrecs,1); depmin=depmen; depmax=depmen;
    for i=1:nrecs
        % skip dataless
        if(isempty(data(i).dep))
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end
        
        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);
        
        % wrap to preserve power if pow2pad<0
        if(pow2pad(i)<0)
            for j=1:nspts(i)
                data(i).dep(j,:)=sum(data(i).dep(j:nspts(i):end,:));
            end
            data(i).dep=data(i).dep(1:nspts(i),:);
        end
        
        % which method?
        switch method{i}
            case 'fft' % periodogram
                data(i).dep=fft(data(i).dep,nspts(i),1);
                data(i).dep=2*data(i).dep.*conj(data(i).dep)/npts(i);
                data(i).dep=delta(i)*data(i).dep(1:(nspts(i)/2+1),:);
        end
        
        % convert to dB & change class back
        data(i).dep=oclass(10*log10(data(i).dep));
        
        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % truncate npts if pow2pad<0
    if(any(pow2pad<0))
        npts(pow2pad<0)=nspts(pow2pad<0);
    end
    
    % update header
    data=changeheader(data,'b',sb,'e',se,'delta',sdelta,'sb',b,...
        'sdelta',delta,'nspts',npts,'npts',nspts/2+1,'iftype','ixy',...
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
