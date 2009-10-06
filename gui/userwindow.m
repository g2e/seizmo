function [data,win,fh]=userwindow(data,fill,func,varargin)
%USERWINDOW    Interactively window SEIZMO records
%
%    Usage:    data=userwindow(data)
%              data=userwindow(data,fill)
%              data=userwindow(data,fill,func)
%              data=userwindow(data,fill,func,'field',value,...)
%              [data,win]=userwindow(...)
%              [data,win,fh]=userwindow(...)
%
%    Description: DATA=USERWINDOW(DATA) presents an interactive menu and
%     plot to facilitate windowing a dataset with a few mouse clicks.  By
%     default, the windowed data is not padded with zeros nor is the mean
%     or trend removed.
%
%     DATA=USERWINDOW(DATA,FILL) toggles zero-padding of records in DATA.
%     If FILL is FALSE or empty (the default), no zero-padding is done.
%     Set FILL to TRUE to pad incomplete records with zeros so that they
%     extend across the window.
%
%     DATA=USERWINDOW(DATA,FILL,FUNC) applies function FUNC to records in
%     DATA after windowing and before the confirmation window.  FUNC must
%     be a function handle.  Some common function handles for this are
%     @removemean and @removetrend.
%
%     DATA=USERWINDOW(DATA,FILL,FUNC,'FIELD',VALUE,...) passes field/value
%     pairs to the plotting function, to allow further customization.
%
%     [DATA,WIN]=USERWINDOW(...) also returns the window limits in WIN as
%     [start end].
%
%     [DATA,WIN,FH]=USERWINDOW(...) returns the figure handles in FH.
%
%    Notes:
%     - automatically switches start and end window limits if they are not
%       proper (ie if the end comes before the start)
%
%    Header changes: B, E, NPTS, DELTA, DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Let the user window the data, adding zeros to fill any incomplete
%     records in the window and remove the trend after windowing:
%      data=userwindow(data,true,@removetrend);
%
%    See also: usertaper, usercluster, selectrecords

%     Version History:
%        Sep.  5, 2009 - rewrite
%        Sep.  9, 2009 - added documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  9, 2009 at 03:45 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% check fill
if(nargin<2 || isempty(fill)); fill=false; end

% check function handle
if(nargin<3 || isempty(func))
    func=@deal;
elseif(~isa(func,'function_handle'))
    error('seizmo:userwindow:badInput','FUNC must be a function handle!');
end

% outer loop - only breaks free by user command
happy_user=false; fh=[-1 -1];
while(~happy_user)
    % explain to the user how this works with a little prompt
    % and make them decide what kind of plot to use for the
    % windowing, offering them back out options.  This prompt
    % looks like trash because of default menu fonts.
    prompt={'+-------------------------------------------------------+'
            '|                Welcome to SEIZMO''s interactive windowing function             |'
            '+-------------------------------------------------------+'
            '|                                                                                                                |'
            '|                                             MOUSE USAGE                                             |'
            '|                                                                                                                |'
            '|    LEFT CLICK                       MIDDLE CLICK                      RIGHT CLICK   |'
            '+-----------+                +------------+               +------------+'
            '|   Mark Window                     Finalize Marks                     Mark Window   |'
            '|          Start                                                                              End           |'
            '+-------------------------------------------------------+'
            '|                                                                                                                |'
            '|                                                   NOTES                                                   |'
            '|                                                                                                                |'
            '|          + You may refine window marks until you finalize                        |'
            '|          + When finalized, a new plot with the windowed                           |'
            '|              waveforms will appear, as well as a confirmation                      |'
            '|              prompt.  You will have the option to re-window.                      |'
            '|                                                                                                                |'
            '+-------------------------------------------------------+'
            '|                                                                                                                |'
            '|                 PLEASE CHOOSE AN OPTION BELOW TO PROCEED!                 |'
            '|                                                                                                                |'
            '+-------------------------------------------------------+'};
    
    
    % way cooler menu -- if only matlab gui's used fixed width
    %{
    prompt={'+-------------------------------------------------------+'
            '|  Welcome to SEIZMO''s interactive windowing function   |'
            '+-------------------------------------------------------+'
            '|                                                       |'
            '|                     MOUSE USAGE:                      |'
            '|                                                       |'
            '| LEFT CLICK          MIDDLE CLICK          RIGHT CLICK |'
            '+------------+      +--------------+      +-------------+'
            '| Mark Window        Finalize Marks         Mark Window |'
            '|   Start                                      End      |'
            '+-------------------------------------------------------+'
            '|                                                       |'
            '|                        NOTES:                         |'
            '|                                                       |'
            '|  + You may refine window marks until you finalize     |'
            '|  + When finalized, a new plot with the windowed       |'
            '|    waveforms will appear, as well as a confirmation   |'
            '|    prompt.  You will have the option to re-window.    |'
            '|                                                       |'
            '+-------------------------------------------------------+'
            '|                                                       |'
            '|        PLEASE CHOOSE AN OPTION BELOW TO PROCEED!      |'
            '|                                                       |'
            '+-------------------------------------------------------+'};
    %}
    
    % display prompt and get user choice
    choice=menu(prompt,'OVERLAY PLOT','EVENLY SPACED PLOT',...
        'DISTANCE SPACED PLOT','DO NOT WINDOW','DIE!');
    
    % proceed by user choice
    switch choice
        case 1 % overlay
            fh(1)=plot2(data,varargin{:});
        case 2 % evenly spaced
            fh(1)=plot0(data,varargin{:});
        case 3 % distance spaced
            fh(1)=recordsection(data,varargin{:});
        case 4 % no window
            win=[];
            return
        case 5 % immediate death
            error('seizmo:userwindow:killYourSelf',...
                'User demanded Seppuku!')
    end
    
    % add window limit markers
    figure(fh(1));
    span=ylim;
    win=xlim;
    hold on
    goh(1)=plot([win(1) win(1)],span,'g','linewidth',4);
    goh(2)=plot([win(2) win(2)],span,'r','linewidth',4);
    hold off
    
    % loop until user finalizes markers
    final=false;
    while(~final)
        % focus plot and let user pick
        figure(fh(1));
        [x,y,button]=ginput(1);
        
        % which mouse button?
        switch button
            case 1
                % left click - update window start
                win(1)=x;
                set(goh(1),'xdata',[x x])
            case 3
                % right click - update window end
                win(2)=x;
                set(goh(2),'xdata',[x x])
            case 2
            % middle click - finalize markers
            if (win(1)>win(2))
                % start and end reversed - fix
                win=win([2 1]);
                set(goh(1),'xdata',[win(1) win(1)])
                set(goh(2),'xdata',[win(2) win(2)])
            end
            final=true;
        end
    end
    
    % get windowed data
    data2=cut(data,'z',win(1),win(2),'fill',fill);
    
    % apply function post cut
    data2=func(data2);

    % proceed by user choice
    switch choice
        case 1 % overlay
            fh(2)=plot2(data2,varargin{:});
        case 2 % evenly spaced
            fh(2)=plot0(data2,varargin{:});
        case 3 % distance spaced
            fh(2)=recordsection(data2,varargin{:});
    end
    
    % confirm results
    choice=menu('KEEP WINDOW?','YES','NO - TRY AGAIN','NO - DIE!');
    switch choice
        case 1 % all done!
            data=data2;
            happy_user=true;
        case 2 % please try again
            close(fh);
            fh=[-1 -1];
        case 3 % i bear too great a shame to go on
            error('seizmo:userwindow:killYourSelf',...
                'User demanded Seppuku!')
    end
end

end
