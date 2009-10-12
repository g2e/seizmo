function [data]=dft(data,format,pow2pad)
%DFT    Performs a discrete fourier transform on SEIZMO data records
%
%    Usage:    data=dft(data)
%              data=dft(data,format)
%              data=dft(data,format,pow2pad)
%
%    Description: DFT(DATA,FORMAT) converts SEIZMO records from the time 
%     domain to the frequency domain using a discrete fourier transform
%     (Matlab's fft).  Following SAC formatting, an option FORMAT can be
%     given to select between real-imaginary and amplitude-phase formats
%     (FORMAT can be either 'amph' or 'rlim' - the default is 'amph').  
%
%     DFT(DATA,FORMAT,POW2PAD) lets the power of 2 zero-padding be adjusted
%     using an integer POW2PAD according to the formula:
%                fftlength=2^(nextpow2(NPTS)+POW2PAD)
%     The default value is 1 and is a good choice for most problems.
%     Setting POW2PAD to <1 is not recommended and setting <0 will truncate
%     the time series.  Setting POW2PAD >1 will increase the frequency
%     resolution of the spectrogram at the cost of increased computation
%     time and array size.
%
%    Notes:
%     - SAC (and thus SEIZMO for sanity) calculates spectral data according
%       to Parseval's theorem.  Dividing records (except for phase!) by 
%       npts*delta/2 gives results that may better match the amplitudes
%       of sinusoid functions.  
%
%    Header Changes: B, SB, E, DELTA, SDELTA, NPTS, NSPTS
%                    DEPMEN, DEPMIN, DEPMAX
%     In the frequency domain B, DELTA, and NPTS are changed to the 
%     beginning frequency, sampling frequency, and number of data points in
%     the transform (includes negative frequencies).  The original values 
%     of B, DELTA, and NPTS are saved in the header as SB, SDELTA, and 
%     NSNPTS and are restored when the IDFT command is performed.
%
%    Examples:
%     To take the derivative of a time-series in the frequency domain:
%      data=idft(multiplyomega(dft(data)))
%
%    See also: IDFT, AMPH2RLIM, RLIM2AMPH, DIVIDEOMEGA, MULTIPLYOMEGA

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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:05 GMT

% todo:

% check nargin
msg=nargchk(1,3,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% valid option values
valid.FORMAT={'amph' 'rlim'};

% defaults
option.FORMAT='amph';
option.POW2PAD=1;

% get options from SEIZMO global
global SEIZMO
try
    fields=fieldnames(SEIZMO.DFT);
    for i=1:numel(fields)
        if(~isempty(SEIZMO.DFT.(fields{i})))
            option.(fields{i})=SEIZMO.DFT.(fields{i});
        end
    end
catch
end

% get option from command line
if(nargin>=2 && ~isempty(format)); option.FORMAT=format; end
if(nargin>=3 && ~isempty(pow2pad)); option.POW2PAD=pow2pad; end

% check options
nrecs=numel(data);
fields=fieldnames(option);
for i=1:numel(fields)
    % specific checks
    switch lower(fields{i})
        case 'format'
            if(iscellstr(option.(fields{i})))
                option.(fields{i})=char(option.(fields{i}));
            end
            if(isempty(option.(fields{i})) ...
                    || ~ischar(option.(fields{i})) ...
                    || ~any(size(option.(fields{i}),1)==[1 nrecs]) ...
                    || ~isempty(setdiff(lower(option.(fields{i})),...
                    valid.(fields{i}))))
                error('seizmo:dft:badInput',...
                    ['%s option must be one of the following:\n'...
                    sprintf('%s ',valid.(fields{i}){:})],fields{i});
            end
            if(size(option.(fields{i}),1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
            option.(fields{i})=cellstr(option.(fields{i}));
        case 'pow2pad'
            if(isempty(option.(fields{i})) || ...
                    any(fix(option.(fields{i}))~=option.(fields{i})) ...
                    || ~any(numel(option.(fields{i}))==[1 nrecs]))
                error('seizmo:dft:badInput',...
                    ['%s option must be an integer or an array\n'...
                    'of integers with one option per record.'],fields{i});
            end
            if(isscalar(option.(fields{i})))
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),1);
            end
    end
end

% special warning for POW2PAD
if(any(option.POW2PAD<0))
    warning('seizmo:dft:dataTruncation',...
        ['Records: ' sprintf('%d ',find(option.POW2PAD<0))...
        '\nSetting option POW2PAD < 0 will\n'...
        'truncate data from the records!'])
end
    
% retreive header info
leven=getlgc(data,'leven');
iftype=getenumdesc(data,'iftype');
[b,delta,npts]=getheader(data,'b','delta','npts');
ncmp=getncmp(data);

% check leven,iftype
if(any(~strcmpi(leven,'true')))
    error('seizmo:dft:illegalOperation',...
        'Illegal operation on unevenly spaced record!')
elseif(any(~strcmpi(iftype,'Time Series File')...
        & ~strcmpi(iftype,'General X vs Y file')))
    error('seizmo:dft:illegalOperation',...
        'Illegal operation on spectral/xyz record!')
end

% output type
iftype(1:nrecs,1)={'Spectral File-Ampl/Phase'};
rlim=strcmpi(option.FORMAT,'rlim');
if(any(rlim)); iftype(rlim)={'Spectral File-Real/Imag'}; end

% loop through records
sb=nan(nrecs,1); se=sb; sdelta=sb; nspts=sb;
depmen=sb; depmin=sb; depmax=sb;
for i=1:nrecs
    % skip dataless (but increase ncmp)
    if(isempty(data(i).dep)); data(i).dep=zeros(0,2*ncmp); continue; end
    
    % save class and convert to double precision
    oclass=str2func(class(data(i).dep));
    data(i).dep=double(data(i).dep);
    
    % get frequency info
    nspts(i)=2^(nextpow2(npts(i))+option.POW2PAD(i));
    sb(i)=0; se(i)=1/(delta(i)*2); sdelta(i)=2*se(i)/nspts(i);
    
    % truncate npts if POW2PAD<0
    if(option.POW2PAD(i)<0)
        npts(i)=2^(nextpow2(npts(i))+option.POW2PAD(i));
    end
    
    % fft
    data(i).dep=delta(i)*fft(data(i).dep,nspts(i));    % SAC compatible
    %data(i).dep=2*fft(data(i).dep,nspts(i))/npts(i);  % for sinusoid amplitudes
    
    % expand data to make room for split
    data(i).dep(:,(1:ncmp)*2)=data(i).dep;
    
    % split complex by desired filetype
    switch lower(option.FORMAT{i})
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
    depmen(i)=mean(data(i).dep(:)); 
    depmin(i)=min(data(i).dep(:)); 
    depmax(i)=max(data(i).dep(:));
end

% update header (note there is no field 'se')
data=changeheader(data,'b',sb,'e',se,'delta',sdelta,'sb',b,...
    'sdelta',delta,'nspts',npts,'npts',nspts,'iftype',iftype,...
    'depmen',depmen,'depmin',depmin,'depmax',depmax);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
