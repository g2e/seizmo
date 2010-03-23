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
%     [DATA,WIN]=USERWINDOW(...) returns a struct WIN with the following
%     fields:
%       WIN.limits  --  limits of the window applied as [START END]
%       WIN.fill    --  fill gaps in data window with zeros (TRUE/FALSE)
%       WIN.func    --  post-windowing function ran on the data
%     Note that the .limits field will be an empty array if no windowing is
%     performed.
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
%    See also: CUT, USERTAPER, USERCLUSTER, SELECTRECORDS

%     Version History:
%        Sep.  5, 2009 - rewrite
%        Sep.  9, 2009 - added documentation
%        Mar.  1, 2010 - updated for newer checking methods
%        Mar. 12, 2010 - pretty text menu for Octave
%        Mar. 15, 2010 - added graphical selection/entry of fill & func,
%                        win is now a struct
%        Mar. 18, 2010 - robust to menu/figure closing
%        Mar. 23, 2010 - preserve last window limits
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 23, 2010 at 01:25 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
end

% attempt windowing
try
    % check fill
    if(nargin<2 || isempty(fill)); fill=false; end

    % check function handle
    if(nargin<3 || isempty(func))
        func=@deal;
    elseif(~isa(func,'function_handle'))
        error('seizmo:userwindow:badInput','FUNC must be a function handle!');
    end
    
    % window parameters
    win.limits=[];
    win.func=func;
    win.fill=fill;

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
        if(strcmpi(getapplication,'OCTAVE'))
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
        end

        % display prompt and get user choice
        choice=menu(prompt,'RESELECT FILL OPTION',...
            'RESELECT POST-WINDOW FUNCTION',...
            'OVERLAY PLOT','EVENLY SPACED PLOT',...
            'DISTANCE SPACED PLOT','DO NOT WINDOW','CRASH!');

        % proceed by user choice
        switch choice
            case 0 % closed menu
                continue;
            case 1 % reselect fill option
                fillstr='NO';
                if(win.fill); fillstr='YES'; end
                j=menu('FILL MISSING DATA WITH ZEROS?',...
                    ['CURRENT (' fillstr ')'],'YES','NO');
                
                % set function
                switch j
                    case 1 % current
                        % leave fill alone
                    case 2 % fill
                        win.fill=true;
                    case 3 % no fill
                        win.fill=false;
                end
                
                % go back to main menu
                continue;
            case 2 % reselect post-window function
                j=menu('SELECT A FUNCTION TO APPLY POST-WINDOW',...
                    ['CURRENT (' upper(func2str(win.func)) ')'],...
                    'NONE (DEAL)','REMOVEMEAN','REMOVETREND','CUSTOM');
                
                % set function
                switch j
                    case 1 % current
                        % leave function alone
                    case 2 % none (use deal)
                        win.func=@deal;
                    case 3 % rmean
                        win.func=@removemean;
                    case 4 % rtrend
                        win.func=@removetrend;
                    case 5 % custom/cmdline
                        tmp=inputdlg(...
                            ['Custom Post-Window Function? [' ...
                            func2str(win.func) ']:'],...
                            'Custom Post-Window Function',1,...
                            {func2str(win.func)});
                        if(~isempty(tmp))
                            try
                                win.func=str2func(tmp{:});
                            catch
                                % do not change win.func
                            end
                        end
                end
                
                % go back to main menu
                continue;
            case 3 % overlay
                fh(1)=plot2(data,varargin{:});
            case 4 % evenly spaced
                fh(1)=plot0(data,varargin{:});
            case 5 % distance spaced
                fh(1)=recordsection(data,varargin{:});
            case 6 % no window
                win.limits=[];
                return;
            case 7 % immediate death
                error('seizmo:userwindow:killYourSelf',...
                    'User demanded Seppuku!');
        end

        % add window limit markers
        figure(fh(1));
        span=ylim;
        if(isempty(win.limits)); win.limits=xlim; end
        hold on
        goh(1)=plot([win.limits(1) win.limits(1)],span,'g','linewidth',4);
        goh(2)=plot([win.limits(2) win.limits(2)],span,'r','linewidth',4);
        hold off

        % loop until user finalizes markers
        final=false;
        while(~final)
            % bring plot to focus (redraw if closed)
            if(~ishandle(fh(1)))
                % redraw (pretty rare to get here)
                switch choice
                    case 3 % overlay
                        fh(1)=plot2(data,varargin{:});
                    case 4 % evenly spaced
                        fh(1)=plot0(data,varargin{:});
                    case 5 % distance spaced
                        fh(1)=recordsection(data,varargin{:});
                end
                span=ylim;
                win.width=xlim;
                hold on
                goh(1)=plot([win.width(1) win.width(1)],span,'g',...
                    'linewidth',4);
                goh(2)=plot([win.width(2) win.width(2)],span,'r',...
                    'linewidth',4);
                hold off
            else
                % bring figure to focus
                figure(fh(1));
            end
            
            % get user click/key
            try
                [x,y,button]=ginput(1);
            catch
                % user closed window - break from loop
                button=2;
            end

            % which mouse button?
            switch button
                case 1
                    % left click - update window start
                    win.limits(1)=x;
                    set(goh(1),'xdata',[x x])
                case 3
                    % right click - update window end
                    win.limits(2)=x;
                    set(goh(2),'xdata',[x x])
                case 2
                    % middle click - finalize markers
                    if (win.limits(1)>win.limits(2))
                        % start and end reversed - fix
                        win.limits=win.limits([2 1]);
                        set(goh(1),'xdata',[win.limits(1) win.limits(1)])
                        set(goh(2),'xdata',[win.limits(2) win.limits(2)])
                    end
                    final=true;
            end
        end

        % get windowed data
        data2=cut(data,'z',win.limits(1),win.limits(2),'fill',win.fill);

        % apply function post cut
        data2=win.func(data2);

        % proceed by user choice
        switch choice
            case 3 % overlay
                fh(2)=plot2(data2,varargin{:});
            case 4 % evenly spaced
                fh(2)=plot0(data2,varargin{:});
            case 5 % distance spaced
                fh(2)=recordsection(data2,varargin{:});
        end

        % confirm results
        choice=0;
        while(~choice)
            choice=menu('KEEP WINDOW?',...
                'YES','NO - TRY AGAIN','NO - CRASH!');
            switch choice
                case 1 % rainbow's end
                    data=data2;
                    happy_user=true;
                case 2 % never never quit!
                    close(fh(ishandle(fh)));
                    fh=[-1 -1];
                case 3 % i bear too great a shame to go on
                    error('seizmo:userwindow:killYourSelf',...
                        'User demanded Seppuku!');
            end
        end
    end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
