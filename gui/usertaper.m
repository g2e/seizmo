function [data,tpr,fh]=usertaper(data,func,varargin)
%USERTAPER    Interactively taper SEIZMO records
%
%    Usage:    data=usertaper(data)
%              data=usertaper(data,func)
%              data=usertaper(data,func,'field',value,...)
%              [data,tpr]=usertaper(...)
%              [data,tpr,fh]=usertaper(...)
%
%    Description: DATA=USERTAPER(DATA) presents an interactive menu and
%     plot interface to taper records in a dataset with a few mouse clicks.
%     The default taper type is that set by function TAPER.  This may be
%     modified using the menu presented.  By default no mean or trend
%     removal is done after tapering.
%
%     DATA=USERTAPER(DATA,FUNC) applies function FUNC to records in DATA
%     after tapering and before the confirmation window.  FUNC must be a
%     function handle.  Some common function handles for this are
%     @removemean and @removetrend.
%
%     DATA=USERTAPER(DATA,FUNC,'FIELD',VALUE,...) passes field/value pairs
%     to the plotting function, to allow further customization.
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
%     [DATA,TPR,FH]=USERTAPER(...) returns the figure handles in FH.
%
%    Notes:
%
%    Header changes: DEPMEN, DEPMIN, DEPMAX
%
%    Examples:
%     Taper and remove the trend afterwards (before confirmation window):
%      data=usertaper(data,@removetrend);
%
%    See also: TAPER, USERWINDOW, USERCLUSTER, SELECTRECORDS

%     Version History:
%        Sep.  9, 2009 - rewrite and added documentation
%        Sep. 23, 2009 - updated for taper changes
%        Mar.  1, 2010 - updated for newer checking methods
%        Mar. 12, 2010 - pretty text menu for Octave
%        Mar. 15, 2010 - added graphical selection/entry of func
%        Mar. 18, 2010 - robust to menu/figure closing
%        Mar. 23, 2010 - preserve last taper widths
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 23, 2010 at 01:25 GMT

% todo:
% - subplot showing taper
% - option for wvtool

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

% attempt tapering
try
    % check function handle
    if(nargin<2 || isempty(func))
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
    tpr.width=[0 0];
    tpr.option=[];
    tpr.func=func;

    % length normalization
    [b,e,npts,delta]=getheader(data,'b','e','npts','delta');
    data=changeheader(data,'b',0,'e',1,'delta',1./(npts-1));

    % outer loop - only breaks free by user command
    happy_user=false; fh=[-1 -1];
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
        choice=menu(prompt,'RESELECT TAPER TYPE',...
            'RESELECT POST-TAPER FUNCTION','OVERLAY PLOT',...
            'EVENLY SPACED PLOT','DISTANCE SPACED PLOT',...
            'DO NOT TAPER','CRASH!');

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
                fh(1)=plot2(data,varargin{:});
            case 4 % evenly spaced
                fh(1)=plot0(data,varargin{:});
            case 5 % distance spaced
                fh(1)=recordsection(data,varargin{:});
            case 6 % no taper
                tpr.type=[];
                tpr.width=[];
                tpr.option=[];
                data=changeheader(data,'b',b,'e',e,'delta',delta);
                return;
            case 7 % immediate death
                error('seizmo:usertaper:killYourSelf',...
                    'User demanded Seppuku!');
        end

        % add taper limit markers
        figure(fh(1));
        span=ylim;
        if(isempty(tpr.width)); tpr.width=xlim; end
        hold on
        goh(1)=plot([tpr.width(1) tpr.width(1)],span,'g','linewidth',4);
        goh(2)=plot([tpr.width(2) tpr.width(2)],span,'r','linewidth',4);
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
                tpr.width=xlim;
                hold on
                goh(1)=plot([tpr.width(1) tpr.width(1)],span,'g',...
                    'linewidth',4);
                goh(2)=plot([tpr.width(2) tpr.width(2)],span,'r',...
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
            end
        end

        % fix end taper
        tpr.width(2)=1-tpr.width(2);

        % get windowed data
        data2=taper(data,tpr.width,0,tpr.type,tpr.option);

        % apply function post cut
        data2=tpr.func(data2);

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
            choice=menu('KEEP TAPER?',...
                'YES','NO - TRY AGAIN','NO - CRASH!');
            switch choice
                case 1 % rainbow's end
                    data=changeheader(data2,'b',b,'e',e,'delta',delta);
                    happy_user=true;
                case 2 % never, never quit!
                    close(fh(ishandle(fh)));
                    fh=[-1 -1];
                case 3 % i bear too great a shame to go on
                    error('seizmo:usertaper:killYourSelf',...
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
