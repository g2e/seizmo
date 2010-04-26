function []=playseizmo(data,factor)
%PLAYSEIZMO    Plays SEIZMO records as sound
%
%    Usage:    playseizmo(data)
%              playseizmo(data,factor)
%
%    Description: PLAYSEIZMO(DATA) plays the records in SEIZMO struct DATA
%     as sound files.  Note that each record is played separately, so this
%     may take a while.  The records are compressed time-wise by a factor
%     of 1000 so that 20Hz corresponds to 20KHz and 20mHz corresponds to
%     20Hz.  We do this because the range 20Hz-20KHz is roughly the range
%     of human hearing and will make much of the seismic spectrum audible.
%
%     PLAYSEIZMO(DATA,FACTOR) time-compresses records by FACTOR.  See the
%     Notes section below for help when using FACTOR to target a specific
%     frequency range.  After scaling by FACTOR, all records are resampled
%     to 44.1KHz to avoid playback issues.  This step unfortunately takes
%     some time.
%
%    Notes:
%     - FACTOR vs Audible Frequency Range (Periods)
%       10000      2mHz-2Hz      (500s-0.5s) -- good for global eq
%       1000       20mHz-20Hz    (50s-0.05s) -- good for regional eq
%       100        200mHz-200Hz  (5s-0.005s) -- good for local eq
%
%    Examples:
%     Lowpass filter at 20s (0.05Hz) and play (5000-5s audible):
%       data=iirfilter(data,'lp','b','c',1/20,'o',2);
%       playseizmo(data,1e5);
%
%     The same but de-emphasize strong signals:
%      data=iirfilter(data,'lp','b','c',1/20,'o',2);
%      data=raise(data,1/4);
%      playseizmo(data,1e5);
%
%    See also: SOUND, SOUNDSC, SEIZMO2WAV

%     Version History:
%        Apr. 25, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 25, 2010 at 12:25 GMT

% todo:

% check nargin
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% check data (dep)
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% try playback
try
    % number of records
    nrecs=numel(data);
    
    % verbosity
    verbose=seizmoverbose;
    
    % check factor
    if(nargin==1 || isempty(nargin)); factor=1000; end
    if(~isreal(factor) || any(factor<=0) || ~any(numel(factor)==[1 nrecs]))
        error('seizmo:playseizmo:badFactor',...
            ['FACTOR must be a positve real scalar or a vector of\n' ...
            'positive reals with one value per record in DATA!']);
    end
    if(isscalar(factor)); factor(1:nrecs,1)=factor; end
    
    % get up/down resampling factors to achieve 44.1KHz
    delta=getheader(data,'delta');
    [n,d]=rrat(delta./factor*44100,1e-6);
    
    % detail message
    if(verbose)
        disp('Playing Record(s) as Audio');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    for i=1:nrecs
        % extract, resample, play
        try
            % scaling between +/-1.1 b/c it seems to avoid popping 4 me
            playsound(resample(double(data(i).dep),n(i),d(i)),44100,1.1);
        catch
            error(lasterror)
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end

function []=playsound(record,rate,scale)
% play record scaled between +/-scale (to avoid popping)
sound(record/(scale*max(abs(record(:)))),rate);
end
