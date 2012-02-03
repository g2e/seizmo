function []=seizmo2wav(data,factor,varargin)
%SEIZMO2WAV    Writes SEIZMO records as WAV files
%
%    Usage:    seizmo2wav(data)
%              seizmo2wav(data,factor)
%              seizmo2wav(data,factor,'option',value,...)
%
%    Description:
%     SEIZMO2WAV(DATA) writes the records in SEIZMO struct DATA as .wav
%     files.  Note that each record is written as a separate wave file (a
%     .wav extension is added to the record names).  The records are
%     compressed time-wise by a factor of 1000 so that 1000s in the record
%     corresponds to 1s of audio.  This is done so that seismic frequencies
%     (20mHz-20Hz) are shifted to the audible frequency range (20Hz-20KHz).
%     To target other seismic frequencies see the next usage form below.
%
%     SEIZMO2WAV(DATA,FACTOR) time-compresses records by FACTOR.  See the
%     Notes section below for help when using FACTOR to target a specific
%     frequency range.  After scaling by FACTOR, all records are resampled
%     to 44.1KHz to avoid playback issues.
%
%     SEIZMO2WAV(DATA,FACTOR,'OPTION',VALUE,...) provides access to pass
%     options to CHANGENAME & CHANGEPATH.  See WRITESEIZMO for a list of
%     the options.
%
%    Notes:
%     - FACTOR vs Audible Frequency Range (Periods)
%       10000      2mHz-2Hz      (500s-0.5s) -- good for global eq
%       1000       20mHz-20Hz    (50s-0.05s) -- good for regional eq
%       100        200mHz-200Hz  (5s-0.005s) -- good for local eq
%
%    Examples:
%     % Lowpass filter at 20s (0.05Hz) and play (5000-5s audible):
%     data=iirfilter(data,'lp','b','c',1/20,'o',2);
%     seizmo2wav(data,1e5);
%
%     % The same but de-emphasize strong signals:
%     data=iirfilter(data,'lp','b','c',1/20,'o',2);
%     data=raise(data,1/4);
%     seizmo2wav(data,1e5);
%
%    See also: PLAYSEIZMO, WAVWRITE, WAVREAD

%     Version History:
%        Apr. 25, 2010 - initial version
%        July 15, 2010 - nargchk fixes
%        Jan. 28, 2012 - doc update, seizmocheck_state bugfix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2012 at 12:25 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>2 && mod(nargin,2))
    error('seizmo:seizmo2wav:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% try playback
try
    % allow filename and filepath changes
    data=writeparameters(data,varargin{:},'nameappend','.wav');
    
    % number of records
    nrecs=numel(data);
    
    % verbosity
    verbose=seizmoverbose;
    
    % check factor
    if(nargin==1 || isempty(nargin)); factor=1000; end
    if(~isreal(factor) || any(factor<=0) || ~any(numel(factor)==[1 nrecs]))
        error('seizmo:seizmo2wav:badFactor',...
            ['FACTOR must be a positve real scalar or a vector of\n' ...
            'positive reals with one value per record in DATA!']);
    end
    if(isscalar(factor)); factor(1:nrecs,1)=factor; end
    
    % get up/down resampling factors to achieve 44.1KHz
    delta=getheader(data,'delta');
    [n,d]=rrat(delta./factor*44100,1e-6);
    
    % detail message
    if(verbose)
        disp('Writing Record(s) as WAVE Audio File(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    for i=1:nrecs
        % extract, resample, play
        try
            % scaling between +/-1.1 b/c it seems to avoid popping 4 me
            wavwrite(...
                rescale(resample(double(data(i).dep),n(i),d(i)),1.1),...
                44100,fullfile(data(i).path,data(i).name));
        catch
            error(lasterror);
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end

function [record]=rescale(record,scale)
% record scaled between +/-scale (to avoid popping)
record=record/(scale*max(abs(record)));
end
