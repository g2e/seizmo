function [onset,ax]=pickstack(data,varargin)
%PICKSTACK    Interactive picking of signal onset in SEIZMO record stack
%
%    Usage:    onset=pickstack(data)
%              onset=pickstack(data,...)
%              [onset,ax]=pickstack(...)
%              [...]=pickstack(ax,...)
%
%    Description:
%     ONSET=PICKSTACK(DATA) presents an interactive menu and plot to
%     facilitate picking the onset of a signal in the stack of records
%     contained in SEIZMO struct DATA.  The returned value ONSET is a
%     relative time value.  Note that the records are essentially passed to
%     ADDRECORDS which performs a sample-by-sample addition without regard
%     to timing.  MAKE SURE ALL SAMPLES ARE TIME ALIGNED IN THIS USAGE CASE
%     or your stack likely will not be useful.  You will likely need to
%     call CUT and/or INTERPOLATE first to make this happen.
%
%     ONSET=PICKSTACK(DATA,...) passes extra arguments on to STACK (and
%     thus on to INTERPOLATE).  This is the typical usage form.  See
%     INTERPOLATE for argument details.
%
%     [ONSET,FH]=PICKSTACK(...) also returns the plot's axis handle as
%     the second output argument AX.
%
%     [...]=PICKSTACK(AX,...) specifies the axes to plot in using the
%     handle AX.
%
%    Notes:
%
%    Examples:
%     % Pick a normalized stack sampled at 5sps from -200 to 200 seconds:
%     onset=pickstack(normalize(data),5,[],-200,200);
%
%    See also: STACK, USERALIGN

%     Version History:
%        Mar. 24, 2010 - initial version
%        Aug. 26, 2010 - nargchk fix, update for new plotting functions,
%                        onset line plots in correct place on redraw
%        Jan.  6, 2011 - use key2zoompan
%        Apr.  3, 2012 - minor doc update, use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 23:55 GMT

% todo:

% check nargin
error(nargchk(1,7,nargin));

% check data structure
if(isseizmo(data,'dep'))
    ax=-1;
else
    ax=data;
    data=varargin{1};
    varargin(1)=[];
    error(seizmocheck(data,'dep'));
end

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt interpolation
try
    % stack dataset
    data=stack(data,varargin{:});
    
    % default offset
    onset=0;
    
    % outer loop - only breaks free on user command
    happy_user=false; firsttime=true;
    while(~happy_user)
        % plot stack if necessary
        if(~ishandle(ax) || firsttime)
            ax=plot1(data,'title','STACKED RECORD','ax',ax);
            span=ylim(ax);
            hold(ax,'on');
            goh=plot(ax,[onset onset],span,'y--','linewidth',2);
            hold(ax,'off');
            firsttime=false;
        end
        
        % get user choice
        choice=menu(...
            ['ADJUST STACKED SIGNAL ONSET FROM ' num2str(onset) ...
            ' SECONDS?'],'YES','NO');
        
        % act on user choice
        switch choice
            case 1 % adjust
                % plot stack if necessary
                if(~ishandle(ax))
                    ax=plot1(data,'title','STACKED RECORD');
                    span=ylim(ax);
                    hold(ax,'on');
                    goh=plot(ax,[onset onset],span,'y--','linewidth',2);
                    hold(ax,'off');
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
                    if(~ishandle(ax))
                        % redraw figure
                        ax=plot1(data,'title','STACKED RECORD');
                        span=ylim(ax);
                        hold(ax,'on');
                        goh=plot(ax,[onset onset],span,'y--',...
                            'linewidth',2);
                        hold(ax,'off');
                    else
                        % bring plot to focus
                        axes(ax);
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
                    else
                        key2zoompan(button,ax);
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
    error(lasterror);
end

end
