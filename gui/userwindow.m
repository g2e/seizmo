function [data,win,ax]=userwindow(data,limits,fill,func,varargin)
%USERWINDOW    Interactively window SEIZMO records
%
%    Usage:    data=userwindow(data)
%              data=userwindow(data,initwin)
%              data=userwindow(data,initwin,fill)
%              data=userwindow(data,initwin,fill,func)
%              data=userwindow(data,initwin,fill,func,'field',value,...)
%              [data,win]=userwindow(...)
%              [data,win,ax]=userwindow(...)
%
%    Description:
%     DATA=USERWINDOW(DATA) presents an interactive menu and plot to
%     facilitate windowing a dataset with a few mouse clicks.  By default,
%     the windowed data is not padded with zeros nor is the mean or trend
%     removed.
%
%     DATA=USERWINDOW(DATA,INITWIN) sets the initial window limits to
%     INITWIN.  INITWIN should be a 1x2 vector of real numbers ordered as
%     [START END].  The values are in relative time.
%
%     DATA=USERWINDOW(DATA,INITWIN,FILL) toggles zero-padding of records in
%     DATA.  If FILL is FALSE or empty (the default), no zero-padding is
%     done.  Set FILL to TRUE to pad incomplete records with zeros so that
%     they extend across the window.
%
%     DATA=USERWINDOW(DATA,INITWIN,FILL,FUNC) applies function FUNC to
%     records in DATA after windowing and before the confirmation menu.
%     FUNC must be a function handle.  Some common function handles for
%     this are @removemean and @removetrend.
%
%     DATA=USERWINDOW(DATA,INITWIN,FILL,FUNC,'FIELD',VALUE,...) passes
%     field/value pairs to the plotting function, to allow further
%     customization.
%
%     [DATA,WIN]=USERWINDOW(...) returns a struct WIN with the following
%     fields:
%       WIN.limits  --  limits of the window applied as [START END]
%       WIN.fill    --  fill gaps in data window with zeros (TRUE/FALSE)
%       WIN.func    --  post-windowing function ran on the data
%     Note that the .limits field will be an empty array if no windowing is
%     performed.
%
%     [DATA,WIN,AX]=USERWINDOW(...) returns the axes handle in AX.
%
%    Notes:
%     - automatically switches start and end window limits if they are not
%       proper (ie if the end comes before the start)
%
%    Header changes: B, E, NPTS, DELTA, DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Let the user window the data, adding zeros to fill any incomplete
%     % records in the window and remove the trend after windowing:
%     data=userwindow(data,true,@removetrend);
%
%    See also: CUT, USERTAPER, USERCLUSTER, SELECTRECORDS, USERWINNOW,
%              USERMOVEOUT

%     Version History:
%        Sep.  5, 2009 - rewrite
%        Sep.  9, 2009 - added documentation
%        Mar.  1, 2010 - updated for newer checking methods
%        Mar. 12, 2010 - pretty text menu for Octave
%        Mar. 15, 2010 - added graphical selection/entry of fill & func,
%                        win is now a struct
%        Mar. 18, 2010 - robust to menu/figure closing
%        Mar. 23, 2010 - preserve last window limits
%        Apr. 22, 2010 - replace crash with exit (but still crash)
%        Aug. 26, 2010 - update for axes plotting output, checkheader fix
%        Nov.  4, 2010 - couple minor bugfixes
%        Nov. 12, 2010 - added initial window argument
%        Jan.  6, 2011 - use key2zoompan
%        Jan. 17, 2011 - altered the menus (no more crashing exit)
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 11:00 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>4 && mod(nargin,2))
    error('seizmo:userwindow:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data,'dep'));

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
    
    % rethrow error
    error(lasterror);
end

% attempt windowing
try
    % check limits
    if(nargin<2)
        limits=[];
    elseif(~isempty(limits) && (numel(limits)~=2 || ~isreal(limits)))
        error('seizmo:userwindow:badInput',...
            'INITWIN must be 1x2 vector as [START END]!');
    end
    
    % check fill
    if(nargin<3 || isempty(fill))
        fill=false;
    elseif(~isscalar(fill) || (~islogical(fill) && ~isnumeric(fill)))
        error('seizmo:userwindow:badInput',...
            'FILL must be TRUE or FALSE!');
    end

    % check function handle
    if(nargin<4 || isempty(func))
        func=@deal;
    elseif(~isa(func,'function_handle'))
        error('seizmo:userwindow:badInput','FUNC must be a function handle!');
    end
    
    % window parameters
    win.limits=limits;
    win.func=func;
    win.fill=fill;

    % outer loop - only breaks free by user command
    happy_user=false; ax=-1; reax={};
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
        choice=menu(prompt,...
            'RESELECT FILL OPTION',...
            'RESELECT POST-WINDOW FUNCTION',...
            'SELECTION IN AN OVERLAY PLOT',...
            'SELECTION IN AN EVENLY SPACED PLOT',...
            'SELECTION IN AN DISTANCE SPACED PLOT',...
            'DO NOT WINDOW');

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
                ax=plot2(data,varargin{:},reax{:});
            case 4 % evenly spaced
                ax=plot0(data,varargin{:},reax{:});
            case 5 % distance spaced
                ax=recordsection(data,varargin{:},reax{:});
            case 6 % no window
                win.limits=[];
                return;
        end
        
        % use this axis
        reax={'ax' ax};

        % add window limit markers
        span=ylim(ax);
        if(isempty(win.limits)); win.limits=xlim(ax); end
        hold(ax,'on');
        goh(1)=plot(ax,[win.limits(1) win.limits(1)],span,'g',...
            'linewidth',4);
        goh(2)=plot(ax,[win.limits(2) win.limits(2)],span,'r',...
            'linewidth',4);
        hold(ax,'off');

        % loop until user finalizes markers
        final=false;
        while(~final)
            % bring plot to focus (redraw if closed)
            if(~ishandle(ax))
                % redraw (pretty rare to get here)
                switch choice
                    case 3 % overlay
                        ax=plot2(data,varargin{:});
                    case 4 % evenly spaced
                        ax=plot0(data,varargin{:});
                    case 5 % distance spaced
                        ax=recordsection(data,varargin{:});
                end
                reax={'ax' ax};
                span=ylim(ax);
                win.limits=xlim(ax);
                hold(ax,'on');
                goh(1)=plot(ax,[win.limits(1) win.limits(1)],span,'g',...
                    'linewidth',4);
                goh(2)=plot(ax,[win.limits(2) win.limits(2)],span,'r',...
                    'linewidth',4);
                hold(ax,'off');
            else
                % bring axes to focus
                axes(ax);
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
                otherwise
                    key2zoompan(button,ax);
            end
        end

        % get windowed data
        data2=cut(data,'z',win.limits(1),win.limits(2),'fill',win.fill);

        % apply function post cut
        data2=win.func(data2);

        % proceed by user choice
        switch choice
            case 3 % overlay
                ax=plot2(data2,varargin{:},reax{:});
            case 4 % evenly spaced
                ax=plot0(data2,varargin{:},reax{:});
            case 5 % distance spaced
                ax=recordsection(data2,varargin{:},reax{:});
        end

        % confirm results
        choice=0;
        while(~choice)
            choice=menu('KEEP WINDOW?',...
                'YES','NO - TRY AGAIN');
            switch choice
                case 1 % rainbow's end
                    data=data2;
                    happy_user=true;
                case 2 % never never quit!
                    if(ishandle(ax))
                        reax={'ax' ax};
                    else
                        reax={};
                        ax=-1;
                    end
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
    error(lasterror);
end

end
