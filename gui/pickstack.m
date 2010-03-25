function [onset,fh]=pickstack(data,varargin)
%PICKSTACK    Interactive picking of signal onset in SEIZMO record stack
%
%    Usage:    onset=pickstack(data)
%              onset=pickstack(data,...)
%              [onset,fh]=pickstack(...)
%
%    Description: ONSET=PICKSTACK(DATA) presents an interactive menu and
%     plot to facilitate picking the onset of a signal in the stack of
%     records contained in SEIZMO struct DATA.  The returned value ONSET is
%     a relative time value.  Note that the records are essentially passed
%     to ADDRECORDS which performs a sample-by-sample addition without
%     regard to timing.  MAKE SURE ALL SAMPLES ARE TIME ALIGNED IN THIS
%     USAGE CASE or your stack likely will not be useful.  You will likely
%     need to call CUT and/or INTERPOLATE first to make this happen.
%
%     ONSET=PICKSTACK(DATA,...) passes extra arguments on to STACK (and
%     thus on to INTERPOLATE).  This is the typical usage form.  See
%     INTERPOLATE for argument details.
%
%     [ONSET,FH]=PICKSTACK(...) also returns the plot's figure handle as
%     the second output argument FH.
%
%    Notes:
%
%    Examples:
%     Pick a normalized stack sampled at 5sps from -200 to 200 seconds:
%      onset=pickstack(normalize(data),5,[],-200,200);
%
%    See also: STACK, USERALIGN

%     Version History:
%        Mar. 24, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2010 at 23:55 GMT

% todo:

% check nargin
msg=nargchk(1,6,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt interpolation
try
    % stack dataset
    data=stack(data,varargin{:});
    
    % default offset
    onset=0;
    
    % outer loop - only breaks free on user command
    happy_user=false; fh=-1;
    while(~happy_user)
        % plot stack if necessary
        if(~ishandle(fh))
            fh=plot1(data,'title','STACKED RECORD');
            span=ylim;
            hold on
            goh=plot([0 0],span,'y--','linewidth',2);
            hold off
        end
        
        % get user choice
        choice=menu(...
            ['ADJUST STACKED SIGNAL ONSET FROM ' num2str(onset) ...
            ' SECONDS?'],'YES','NO');
        
        % act on user choice
        switch choice
            case 1 % adjust
                % plot stack if necessary
                if(~ishandle(fh))
                    fh=plot1(data,'title','STACKED RECORD');
                    span=ylim;
                    hold on
                    goh=plot([0 0],span,'y--','linewidth',2);
                    hold off
                end
                
                % menu telling user how to interactively adjust onset
                prompt={'+-------------------------------------------------------+'
                    '|                 Welcome to SEIZMO''s interactive picking function                |'
                    '+-------------------------------------------------------+'
                    '|                                                                                                               |'
                    '|                                            MOUSE USAGE                                             |'
                    '|                                                                                                               |'
                    '|    LEFT CLICK                      MIDDLE CLICK                      RIGHT CLICK   |'
                    '+-------------+          +-------------+          +--------------+'
                    '|      Set Signal                      Finalize Onset                                              |'
                    '|        Onset                                                                                              |'
                    '+-------------------------------------------------------+'};
                % way cooler menu -- if only matlab gui's used fixed width
                if(strcmpi(getapplication,'OCTAVE'))
                    prompt={'+-------------------------------------------------------+'
                        '|   Welcome to SEIZMO''s interactive picking function   |'
                        '+-------------------------------------------------------+'
                        '|                                                       |'
                        '|                     MOUSE USAGE:                      |'
                        '|                                                       |'
                        '|   LEFT CLICK        MIDDLE CLICK         RIGHT CLICK  |'
                        '+---------------+   +--------------+    +---------------+'
                        '|   Set Signal       Finalize Onset                     |'
                        '|     Onset                                             |'
                        '+-------------------------------------------------------+'};
                end
                menu(prompt,'I''M READY!');
                
                % loop until right click
                button=1;
                while (button~=2)
                    % bring plot to focus (redraw if closed)
                    if(~ishandle(fh))
                        % redraw figure
                        fh=plot1(data,'title','STACKED RECORD');
                        span=ylim;
                        hold on
                        goh=plot([0 0],span,'y--','linewidth',2);
                        hold off
                    else
                        % bring figure to focus
                        figure(fh);
                    end

                    % get user click/key
                    try
                        [x,y,button]=ginput(1);
                    catch
                        % user closed window - break from loop
                        button=2;
                    end

                    % action on left click
                    if (button==1)
                        onset=x;
                        set(goh,'xdata',[x x]);
                    end
                end
            case 2 % exit
                happy_user=true;
        end
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
