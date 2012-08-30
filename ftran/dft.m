function [data]=dft(data,format,pow2pad)
%DFT    Performs a discrete fourier transform on SEIZMO data records
%
%    Usage:    data=dft(data)
%              data=dft(data,format)
%              data=dft(data,format,pow2pad)
%
%    Description:
%     DFT(DATA,FORMAT) converts SEIZMO records from the time domain to the
%     frequency domain using a discrete fourier transform (Matlab's fft).
%     Following SAC formatting, an option FORMAT can be given to select
%     between real-imaginary and amplitude-phase formats (FORMAT can be
%     either 'amph' or 'rlim' - the default is 'amph').
%
%     DFT(DATA,FORMAT,POW2PAD) lets the power of 2 zero-padding be adjusted
%     using an integer POW2PAD according to the formula:
%                fftlength=2^(nextpow2(NPTS)+POW2PAD)
%     The default value is 0 and is a good choice for most problems.
%     Setting POW2PAD to <0 is not recommended and will truncate the time
%     series.  Setting POW2PAD >0 will increase the frequency resolution of
%     the spectrogram at the cost of increased computation time and array
%     size.
%
%    Notes:
%     - SAC (and thus SEIZMO for sanity) calculates spectral data according
%       to Parseval's theorem.  Dividing records (except for phase!) by 
%       npts*delta/2 gives results that may better match the amplitudes
%       of sinusoid functions.  
%
%    Header Changes: B, SB, E, DELTA, SDELTA, NPTS, NSPTS, IFTYPE
%                    DEPMEN, DEPMIN, DEPMAX
%     In the frequency domain B, DELTA, and NPTS are changed to the 
%     beginning frequency, sampling frequency, and number of data points in
%     the transform (includes negative frequencies).  The original values 
%     of B, DELTA, and NPTS are saved in the header as SB, SDELTA, and 
%     NSPTS and are restored when the IDFT command is performed.
%
%    Examples:
%     % To take the derivative of a time-series in the frequency domain:
%     data=idft(omegamultiply(dft(data)))
%
%    See also: IDFT, AMPH2RLIM, RLIM2AMPH, OMEGADIVIDE, OMEGAMULTIPLY,
%              OMEGASHIFT, OMEGAHILBERT, OMEGAANALYTIC, OMEGAGAUSSIAN

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 28, 2008 - seischk support and code cleanup
%        Mar.  2, 2008 - readibility change
%        Mar.  3, 2008 - now compatible with SAC
%        Mar.  4, 2008 - cleaned up errors and warnings
%        May  12, 2008 - name changed to dft
%        June 11, 2008 - added example
%        June 29, 2008 - doc update, .dep rather than .x
%        July 18, 2008 - dataless support, one ch call, updates DEP* fields
%        Oct.  6, 2008 - minor code cleaning
%        Nov. 22, 2008 - update for new name schema
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June  3, 2009 - allow individual record options, warn on
%                        pow2pad<0, better checking, global option access
%        June 11, 2009 - fix bug that assigned delta for all file the same
%        Oct. 15, 2009 - force fft down columns
%        Jan. 29, 2010 - seizmoverbose support, proper SEIZMO handling,
%                        improved messaging, drop GETNCMP usage
%        Feb. 11, 2011 - mass nargchk fix, mass seizmocheck fix
%        Dec. 21, 2011 - doc update, better checkheader usage, drop globals
%        Feb.  4, 2012 - minor doc update
%        May  29, 2012 - pow2pad=0 by default
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  29, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt fast fourier transform
try
    % check headers
    data=checkheader(data,...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR');
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);

    % valid format values
    valid={'amph' 'rlim'};

    % defaults
    if(nargin<2 || isempty(format)); format='amph'; end
    if(nargin<3 || isempty(pow2pad)); pow2pad=0; end

    % check options
    if(ischar(format)); format=cellstr(format); end
    if(~iscellstr(format) || ~any(numel(format)==[1 nrecs]) ...
            || any(~ismember(lower(format),valid)))
        error('seizmo:dft:badInput',...
            ['FORMAT option must be one of the following:\n'...
            sprintf('%s ',valid{:})]);
    end
    if(isscalar(format)); format(1:nrecs,1)=format; end
    format=format(:);
    if(~isnumeric(pow2pad) || ~isreal(pow2pad) ...
            || any(pow2pad~=fix(pow2pad)) ...
            || ~any(numel(pow2pad)==[1 nrecs]))
        error('seizmo:dft:badInput',...
            ['POW2PAD must be an integer or an array\n' ...
            'of integers (one per record)!']);
    end
    if(isscalar(pow2pad)); pow2pad(1:nrecs,1)=pow2pad; end
    pow2pad=pow2pad(:);

    % special warning for POW2PAD
    if(any(pow2pad<0))
        warning('seizmo:dft:dataTruncation',...
            ['Setting option POW2PAD < 0 will\n' ...
            'truncate data from the records!\n' ...
            'Record(s):\n' sprintf('%d ',find(pow2pad<0))])
    end

    % retreive header info
    [b,delta,npts,ncmp]=getheader(data,'b','delta','npts','ncmp');

    % output type
    iftype(1:nrecs,1)={'iamph'};
    rlim=strcmpi(format,'rlim');
    if(any(rlim)); iftype(rlim)={'irlim'}; end
    
    % detail message
    if(verbose)
        disp('Transforming Record(s) to the Frequency Domain');
        print_time_left(0,nrecs);
    end

    % loop through records
    sb=nan(nrecs,1); se=sb; sdelta=sb; nspts=sb;
    depmen=sb; depmin=sb; depmax=sb;
    for i=1:nrecs
        % skip dataless (but increase ncmp)
        if(isempty(data(i).dep))
            data(i).dep=zeros(0,2*ncmp);
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % save class and convert to double precision
        oclass=str2func(class(data(i).dep));
        data(i).dep=double(data(i).dep);

        % get frequency info
        nspts(i)=2^(nextpow2(npts(i))+pow2pad(i));
        sb(i)=0; se(i)=1/(delta(i)*2); sdelta(i)=2*se(i)/nspts(i);

        % truncate npts if POW2PAD<0
        if(pow2pad(i)<0)
            npts(i)=2^(nextpow2(npts(i))+pow2pad(i));
        end

        % fft
        data(i).dep=delta(i)*fft(data(i).dep,nspts(i),1);    % SAC compatible
        %data(i).dep=2*fft(data(i).dep,nspts(i),1)/npts(i);  % for sinusoid amplitudes

        % expand data to make room for split
        data(i).dep(:,(1:ncmp)*2)=data(i).dep;

        % split complex by desired filetype
        switch lower(format{i})
            case 'rlim'
                data(i).dep(:,1:2:end)=real(data(i).dep(:,2:2:end));
                data(i).dep(:,2:2:end)=imag(data(i).dep(:,2:2:end));
            case 'amph'
                data(i).dep(:,1:2:end)=abs(data(i).dep(:,2:2:end));
                data(i).dep(:,2:2:end)=angle(data(i).dep(:,2:2:end));
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

    % update header (note there is no field 'se')
    data=changeheader(data,'b',sb,'e',se,'delta',sdelta,'sb',b,...
        'sdelta',delta,'nspts',npts,'npts',nspts,'iftype',iftype,...
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
