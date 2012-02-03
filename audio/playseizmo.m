function []=playseizmo(data,speedup,playrate)
%PLAYSEIZMO    Plays SEIZMO records as sound
%
%    Usage:    playseizmo(data)
%              playseizmo(data,speedup)
%              playseizmo(data,speedup,playrate)
%
%    Description:
%     PLAYSEIZMO(DATA) plays the records in SEIZMO struct DATA as sound
%     files.  Note that each record is played separately, so this may take
%     a while.  The records are sped up by a factor of 1000 so that 20Hz
%     corresponds to 20KHz and 20mHz corresponds to 20Hz.  We do this
%     because the range 20Hz-20KHz is roughly the range of human hearing
%     and will make much of the seismic spectrum audible.
%
%     PLAYSEIZMO(DATA,SPEEDUP) speeds up records by SPEEDUP.  See the
%     Notes section below for help when using SPEEDUP to target a specific
%     frequency range.  After scaling by SPEEDUP, all records are resampled
%     to 44.1KHz to avoid playback issues.  This step unfortunately takes
%     some time.
%
%     PLAYSEIZMO(DATA,SPEEDUP,PLAYRATE) adjusts the sample rate of the
%     sound output.  The default is 44.1KHz.  Other rates that might be
%     more compatible: 8192Hz, 48KHz.  PLAYRATE is in Hz.
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
%     playseizmo(data,1e5);
%
%     % The same but de-emphasize strong signals:
%     data=iirfilter(data,'lp','b','c',1/20,'o',2);
%     data=raise(data,1/4);
%     playseizmo(data,1e5);
%
%    See also: SOUND, SOUNDSC, SEIZMO2WAV

%     Version History:
%        Apr. 25, 2010 - initial version
%        July 15, 2010 - add playrate arg, doc update
%        Jan. 28, 2012 - doc update, seizmocheck_state bugfix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2012 at 09:25 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% try playback
try
    % number of records
    nrecs=numel(data);
    
    % verbosity
    verbose=seizmoverbose;
    
    % check speedup
    if(nargin<2 || isempty(speedup)); speedup=1000; end
    if(nargin<3 || isempty(playrate)); playrate=44100; end
    if(~isreal(speedup) || any(speedup<=0) ...
            || ~any(numel(speedup)==[1 nrecs]))
        error('seizmo:playseizmo:badSpeedUp',...
            ['SPEEDUP must be a positve real scalar or a vector of\n' ...
            'positive reals with one value per record in DATA!']);
    end
    if(~isreal(playrate) || any(playrate<=0) ...
            || ~any(numel(playrate)==[1 nrecs]))
        error('seizmo:playseizmo:badPlayRate',...
            ['PLAYRATE must be a positve real scalar or a vector of\n' ...
            'positive reals with one value per record in DATA!']);
    end
    if(isscalar(speedup)); speedup(1:nrecs,1)=speedup; end
    if(isscalar(playrate)); playrate(1:nrecs,1)=playrate; end
    
    % get up/down resampling factors to achieve 44.1KHz
    delta=getheader(data,'delta');
    [n,d]=rrat(delta./speedup*playrate,1e-6);
    
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
            playsound(resample(double(data(i).dep),...
                n(i),d(i)),playrate(i),1.1);
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

function []=playsound(record,rate,scale)
% play record scaled between +/-scale (to avoid popping)
playblocking(audioplayer(record/(scale*max(abs(record(:)))),rate,16));
%sound(record/(scale*max(abs(record(:)))),rate);
end
