function [data,tpr,ax]=usertaper(data,width,func,varargin)
%USERTAPER    Interactively taper SEIZMO records
%
%    Usage:    data=usertaper(data)
%              data=usertaper(data,width)
%              data=usertaper(data,width,func)
%              data=usertaper(data,width,func,'field',value,...)
%              [data,tpr]=usertaper(...)
%              [data,tpr,ax]=usertaper(...)
%
%    Description:
%     DATA=USERTAPER(DATA) presents an interactive menu and plot interface
%     to taper records in a dataset with a few mouse clicks.  The default
%     taper type is that set by function TAPER.  This may be modified using
%     the menu presented.  By default no mean or trend removal is done
%     after tapering.
%
%     DATA=USERTAPER(DATA,WIDTH) alters the default taper width.  Note that
%     the default here is [0 0] which means no leading or trailing tapers. 
%     This differs from the default WIDTH in the function TAPER.  For more
%     details on how WIDTH is specified see there.
%
%     DATA=USERTAPER(DATA,WIDTH,FUNC) applies function FUNC to records in
%     DATA after tapering and before the confirmation window.  FUNC must be
%     a function handle.  Some common function handles for this are
%     @removemean and @removetrend.
%
%     DATA=USERTAPER(DATA,WIDTH,FUNC,'FIELD',VALUE,...) passes field/value
%     pairs to the plotting function, to allow further customization.
%
%     [DATA,TPR]=USERTAPER(...) returns a struct TPR with the following
%     fields:
%      TPR.type    --  type of taper utilized (string)
%      TPR.width   --  width of tapers applied as [LEADING TRAILING]
%      TPR.option  --  taper option if applicible
%      TPR.func    --  post-taper function ran on the data
%     Note that the .type & .option will be empty unless set interactively
%     (see TAPER for defaults).  The .width field will be an empty array if
%     no tapering is performed.
%
%     [DATA,TPR,AX]=USERTAPER(...) returns the axes handle in AX.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     % Taper and remove the trend afterwards (before confirmation window):
%     data=usertaper(data,[],@removetrend);
%
%    See also: TAPER, USERWINDOW, USERCLUSTER, SELECTRECORDS, USERWINNOW,
%              USERMOVEOUT

%     Version History:
%        Sep.  9, 2009 - rewrite and added documentation
%        Sep. 23, 2009 - updated for taper changes
%        Mar.  1, 2010 - updated for newer checking methods
%        Mar. 12, 2010 - pretty text menu for Octave
%        Mar. 15, 2010 - added graphical selection/entry of func
%        Mar. 18, 2010 - robust to menu/figure closing
%        Mar. 23, 2010 - preserve last taper widths
%        Apr. 22, 2010 - replace crash with exit (but still crash)
%        Aug. 26, 2010 - update for axes plotting output, checkheader fix
%        Jan.  6, 2011 - use key2zoompan
%        Jan. 17, 2011 - allow specifying the default taper width, altered
%                        the menus (no more crashing exit)
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 11:00 GMT

% todo:
% - subplot showing taper
% - option for wvtool

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>3 && ~mod(nargin,2))
    error('seizmo:usertaper:badNumInputs',...
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

% attempt tapering
try
    % check taper width
    if(nargin<2 || isempty(width))
        width=[0 0];
    elseif(numel(width)>2 || ~isreal(width))
        error('seizmo:usertaper:badInput',...
            'WIDTH must be a scalar or 1x2 vector of taper widths!');
    end
    
    % force into untapered limits
    % - this is how I work with it internally
    if(isscalar(width))
        width=[width 1-width];
    else
        width(2)=1-width(2);
    end
    
    % check function handle
    if(nargin<3 || isempty(func))
        func=@deal;
    elseif(~isa(func,'function_handle'))
        error('seizmo:usertaper:badInput',...
            'FUNC must be a function handle!');
    end

    % taper types
    types={[],'barthannwin','bartlett','blackman','blackmanharris',...
        'bohmanwin','chebwin','flattopwin','gausswin','hamming','hann',...
        'kaiser','nuttallwin','parzenwin','rectwin','triang','tukeywin'};

    % taper parameters
    tpr.type=[];
    tpr.width=width;
    tpr.option=[];
    tpr.func=func;
    
    % fix x axis label
    varargin=[{'xlabel' 'Normalized Length'} varargin];

    % length normalization
    [b,e,npts,delta]=getheader(data,'b','e','npts','delta');
    data=changeheader(data,'b',0,'e',1,'delta',1./(npts-1));

    % outer loop - only breaks free by user command
    happy_user=false; ax=-1; reax={};
    while(~happy_user)
        % explain to the user how this works with a little prompt
        % and make them decide what kind of plot to use for the
        % windowing, offering them back out options.  This prompt
        % looks like trash because of default menu fonts.
        prompt={'+-------------------------------------------------------+'
                '|                Welcome to SEIZMO''s interactive tapering function                |'
                '+-------------------------------------------------------+'
                '|                                                                                                               |'
                '|                                            MOUSE USAGE                                             |'
                '|                                                                                                               |'
                '|    LEFT CLICK                      MIDDLE CLICK                      RIGHT CLICK   |'
                '+-------------+          +-------------+          +--------------+'
                '|  Mark Untapered                 Finalize Marks                  Mark Untapered  |'
                '|          Start                                                                             End           |'
                '+-------------------------------------------------------+'
                '|                                                                                                               |'
                '|                                                  NOTES                                                   |'
                '|                                                                                                               |'
                '|          + You may refine tapering limits until you finalize                       |'
                '|          + When finalized, a new plot with the tapered                              |'
                '|              waveforms will appear, as well as a confirmation                     |'
                '|              prompt.  You will have the option to re-taper.                         |'
                '|                                                                                                               |'
                '+-------------------------------------------------------+'
                '|                                                                                                               |'
                '|                PLEASE CHOOSE AN OPTION BELOW TO PROCEED!                 |'
                '|                                                                                                               |'
                '+-------------------------------------------------------+'};


        % way cooler menu -- if only matlab gui's used fixed width
        if(strcmpi(getapplication,'OCTAVE'))
            prompt={'+-------------------------------------------------------+'
                    '|   Welcome to SEIZMO''s interactive tapering function   |'
                    '+-------------------------------------------------------+'
                    '|                                                       |'
                    '|                     MOUSE USAGE:                      |'
                    '|                                                       |'
                    '|   LEFT CLICK        MIDDLE CLICK         RIGHT CLICK  |'
                    '+---------------+   +--------------+    +---------------+'
                    '| Mark Untapered     Finalize Marks      Mark Untapered |'
                    '|     Start                                    End      |'
                    '+-------------------------------------------------------+'
                    '|                                                       |'
                    '|                        NOTES:                         |'
                    '|                                                       |'
                    '|  + You may refine tapering limits until you finalize  |'
                    '|  + When finalized, a new plot with the tapered        |'
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
            'RESELECT TAPER TYPE',...
            'RESELECT POST-TAPER FUNCTION',...
            'SELECTION IN AN OVERLAY PLOT',...
            'SELECTION IN AN EVENLY SPACED PLOT',...
            'SELECTION IN AN DISTANCE SPACED PLOT',...
            'DO NOT TAPER');

        % proceed by user choice
        switch choice
            case 0 % closed menu
                continue;
            case 1 % reselect taper type
                j=menu('SELECT A TAPER TYPE','DEFAULT','BARTHANN',...
                    'BARTLETT','BLACKMAN','BLACKMAN-HARRIS','BOHMAN',...
                    'CHEBYCHEV','FLAT TOP','GAUSSIAN','HAMMING','HANN',...
                    'KAISER','NUTTALL','PARZEN','RECTANGULAR',...
                    'TRIANGULAR','TUKEY');
                if(j); tpr.type=types{j}; end

                % taper options
                switch j
                    case 7 % cheb
                        tpr.option=inputdlg(['Chebychev - Stopband '...
                            'Attenuation (in dB)? [100]:'],...
                            'Chebychev Taper Option',1,{'100'});
                        if(isempty(tpr.option))
                            tpr.option=100;
                        else
                            tpr.option=str2double(tpr.option{:});
                        end
                        if(isempty(tpr.option) || isnan(tpr.option))
                            tpr.option=100;
                        end
                    case 9 % gauss
                        tpr.option=inputdlg(...
                            'Gaussian - Number of std dev? [3.5]:',...
                            'Gaussian Taper Option',1,{'3.5'});
                        if(isempty(tpr.option))
                            tpr.option=3.5;
                        else
                            tpr.option=str2double(tpr.option{:});
                        end
                        if(isempty(tpr.option) || isnan(tpr.option))
                            tpr.option=3.5;
                        end
                    case 12 % kaiser
                        tpr.option=inputdlg(['Kaiser - Stopband '...
                            'Attenuation Factor Beta? [7.5]:'],...
                            'Kaiser Taper Option',1,{'7.5'});
                        if(isempty(tpr.option))
                            tpr.option=7.5;
                        else
                            tpr.option=str2double(tpr.option{:});
                        end
                        if(isempty(tpr.option) || isnan(tpr.option))
                            tpr.option=7.5;
                        end
                    case 17 % tukey
                        tpr.option=inputdlg(...
                            'Tukey - Taper to Constant Ratio? [1]:',...
                            'Tukey Taper Option',1,{'1'});
                        if(isempty(tpr.option))
                            tpr.option=1;
                        else
                            tpr.option=str2double(tpr.option{:});
                        end
                        if(isempty(tpr.option) || isnan(tpr.option))
                            tpr.option=1;
                        end
                end
                % go back to main menu
                continue;
            case 2 % reselect post-taper function
                j=menu('SELECT A FUNCTION TO APPLY POST-TAPER',...
                    ['CURRENT (' upper(func2str(tpr.func)) ')'],...
                    'NONE (DEAL)','REMOVEMEAN','REMOVETREND','CUSTOM');
                
                % set function
                switch j
                    case 1 % current
                        % leave function alone
                    case 2 % none (use deal)
                        tpr.func=@deal;
                    case 3 % rmean
                        tpr.func=@removemean;
                    case 4 % rtrend
                        tpr.func=@removetrend;
                    case 5 % custom/cmdline
                        tmp=inputdlg(...
                            ['Custom Post-Taper Function? [' ...
                            func2str(tpr.func) ']:'],...
                            'Custom Post-Taper Function',1,...
                            {func2str(tpr.func)});
                        if(~isempty(tmp))
                            try
                                tpr.func=str2func(tmp{:});
                            catch
                                % do not change tpr.func
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
            case 6 % no taper
                tpr.type=[];
                tpr.width=[];
                tpr.option=[];
                data=changeheader(data,'b',b,'e',e,'delta',delta);
                return;
        end
        
        % use this axis
        reax={'ax' ax};

        % add taper limit markers
        span=ylim(ax);
        if(isempty(tpr.width)); tpr.width=xlim; end
        hold(ax,'on');
        goh(1)=plot(ax,[tpr.width(1) tpr.width(1)],span,'g',...
            'linewidth',4);
        goh(2)=plot(ax,[tpr.width(2) tpr.width(2)],span,'r',...
            'linewidth',4);
        hold(ax,'off')

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
                tpr.width=xlim(ax);
                hold(ax,'on');
                goh(1)=plot(ax,[tpr.width(1) tpr.width(1)],span,'g',...
                    'linewidth',4);
                goh(2)=plot(ax,[tpr.width(2) tpr.width(2)],span,'r',...
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
                    % left click - update untapered start
                    tpr.width(1)=x;
                    set(goh(1),'xdata',[x x])
                case 3
                    % right click - update untapered end
                    tpr.width(2)=x;
                    set(goh(2),'xdata',[x x])
                case 2
                    % middle click - finalize markers
                    if (tpr.width(1)>tpr.width(2))
                        % start and end reversed - fix
                        tpr.width=tpr.width([2 1]);
                        set(goh(1),'xdata',[tpr.width(1) tpr.width(1)])
                        set(goh(2),'xdata',[tpr.width(2) tpr.width(2)])
                    end
                    final=true;
                otherwise
                    key2zoompan(button,ax);
            end
        end

        % get tapered data
        tpr.width(2)=1-tpr.width(2);
        data2=taper(data,tpr.width,0,tpr.type,tpr.option);
        tpr.width(2)=1-tpr.width(2);

        % apply function post cut
        data2=tpr.func(data2);

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
            choice=menu('KEEP TAPER?',...
                'YES','NO - TRY AGAIN');
            switch choice
                case 1 % rainbow's end
                    data=changeheader(data2,'b',b,'e',e,'delta',delta);
                    happy_user=true;
                case 2 % never, never quit!
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
