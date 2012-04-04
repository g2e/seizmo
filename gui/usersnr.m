function [snr,s,ax]=usersnr(data,nwin,swin,method,varargin)
%USERSNR    Interactively select windows for SNR estimation
%
%    Usage:    snr=usersnr(data)
%              snr=usersnr(data,noisewin,signalwin)
%              snr=usersnr(data,noisewin,signalwin,method)
%              snr=usersnr(data,noisewin,signalwin,method,'param',val,...)
%              [snr,s]=usersnr(...)
%              [snr,s,ax]=usersnr(...)
%
%    Description:
%     SNR=USERSNR(DATA) presents an interactive menu and plot to facilitate
%     signal-to-noise ratio estimation of records in SEIZMO struct DATA.
%     The estimates are output in an Nx1 array SNR where N is the number of
%     records in DATA.  The default window limits extend across the timing
%     of the entire dataset, so it is necessary to adjust these limits
%     (otherwise all SNR estimates will be 1).  The default method is
%     'PEAK2PEAK'.  The defaults may be adjusted using the alternative
%     usage forms below.
%
%     SNR=USERSNR(DATA,NOISEWIN,SIGNALWIN) specifies the noise and signal
%     window limits.  Both NOISEWIN & SIGNALWIN must be 1x2 real-valued
%     arrays indicating [start_time end_time] relative to the reference
%     time(s) of the record(s) in DATA.  The defaults extend across the
%     timing of the entire dataset.  These may be changed by the user (they
%     are just the initial values).
%
%     SNR=USERSNR(DATA,NOISEWIN,SIGNALWIN,METHOD) specifies the SNR
%     estimation method to use.  METHOD must be 'PEAK2PEAK', 'RMS',
%     'ROBUSTRMS', 'PEAK2RMS', or 'PEAK2ROBUSTRMS'.  The default is
%     'PEAK2PEAK'.  See QUICKSNR for more details.  The method may be
%     changed by the user (this is the initial value).
%
%     SNR=USERSNR(DATA,NOISEWIN,SIGNALWIN,METHOD,'PARAM',VAL,...) passes
%     parameter/value pairs to the underlying plotting function to allow
%     customization.
%
%     [SNR,S]=USERSNR(...) returns a struct S with the following fields:
%      S.noisewin   --  time limits of the noise window (relative times) 
%      S.signalwin  --  time limits of the signal window (relative times)
%      S.method     --  SNR estimation method
%      S.plottype   --  function handle of plotting function
%
%     [SNR,S,AX]=USERSNR(...) returns the axes handle, AX, of the plot.
%
%    Notes:
%
%    Examples:
%     % Specify the default window limits and method
%     % and let the user modify them if they desire:
%     [snr,s,ax]=usersnr(data,[-90 -15],[-15 60],'robustrms');
%
%    See also: QUICKSNR, USERWINDOW, CUT

%     Version History:
%        Mar. 18, 2010 - initial version
%        Apr. 21, 2010 - set plot name & add additional methods
%        Aug. 26, 2010 - update for axes plotting output, checkheader fix
%        Jan.  6, 2011 - use key2zoompan
%        Jan. 18, 2011 - remove exit button
%        Jan. 23, 2011 - redraw plot if plot type selected
%        Apr.  3, 2012 - minor doc update, use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 11:00 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>4 && mod(nargin,2))
    error('seizmo:usersnr:badNumInputs',...
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

% attempt snr
try
    % get header values
    [b,e]=getheader(data,'b','e');
    
    % default window limits
    b=min(b);
    e=max(e);
    
    % default options
    if(nargin<2 || isempty(nwin)); nwin=[b e]; end
    if(nargin<3 || isempty(swin)); swin=[b e]; end
    if(nargin<4 || isempty(method)); method='peak2peak'; end
    
    % check options
    snrmethods={'PEAK2PEAK' 'RMS' 'ROBUSTRMS' 'PEAK2RMS' 'PEAK2ROBUSTRMS'};
    if(~isequal(size(nwin),[1 2]) || ~isreal(nwin))
        error('seizmo:usersnr:badInput',...
            'NOISEWIN must be a 1x2 real-valued array!');
    elseif(~isequal(size(swin),[1 2]) || ~isreal(swin))
        error('seizmo:usersnr:badInput',...
            'SIGNALWIN must be a 1x2 real-valued array!');
    elseif(~ischar(method) || size(method,1)~=1 || ...
            ~any(strcmpi(method,snrmethods)))
        error('seizmo:usersnr:badInput',...
            ['METHOD must be one of the following:\n' ...
            sprintf('%s ',snrmethods{:})]);
    end
    
    % assign to output struct
    s.noisewin=nwin;
    s.signalwin=swin;
    s.method=method;
    s.plottype=@plot0;
    
    % loop until user is satisfied
    happy_user=false; ax=-1; reax={}; newplot=false;
    while(~happy_user)
        % create prompt to explain to the user how window selection works
        % This prompt looks like trash because of default menu fonts.
        prompt={'+-------------------------------------------------------+'
                '|            Welcome to SEIZMO''s interactive SNR estimation function         |'
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
                '|                                                                                                                |'
                '+-------------------------------------------------------+'
                '|                                                                                                                |'
                '|                 PLEASE CHOOSE AN OPTION BELOW TO PROCEED!                 |'
                '|                                                                                                                |'
                '+-------------------------------------------------------+'};

        % way cooler menu -- if only matlab gui's used fixed width
        if(strcmpi(getapplication,'OCTAVE'))
            prompt={'+----------------------------------------------------------+'
                    '| Welcome to SEIZMO''s interactive SNR estimation function |'
                    '+----------------------------------------------------------+'
                    '|                                                          |'
                    '|                       MOUSE USAGE:                       |'
                    '|                                                          |'
                    '| LEFT CLICK            MIDDLE CLICK           RIGHT CLICK |'
                    '+------------+        +--------------+       +-------------+'
                    '| Mark Window          Finalize Marks          Mark Window |'
                    '|    Start                                        End      |'
                    '+----------------------------------------------------------+'
                    '|                                                          |'
                    '|                          NOTES:                          |'
                    '|                                                          |'
                    '|   + You may refine window marks until you finalize       |'
                    '|                                                          |'
                    '+----------------------------------------------------------+'
                    '|                                                          |'
                    '|          PLEASE CHOOSE AN OPTION BELOW TO PROCEED!       |'
                    '|                                                          |'
                    '+----------------------------------------------------------+'};
        end
        
        % plot only if new plot type or no plot
        if(~ishandle(ax) || newplot)
            % make title
            ptitle={['SNR Estimation Method: ' upper(s.method)]
                'Yellow Dashed Line --  Noise Window Start'
                '  Blue Dashed Line --  Noise Window End  '
                '  Green Solid Line -- Signal Window Start'
                '    Red Solid Line -- Signal Window End  '};
            
            % plot records
            ax=s.plottype(data,'title',ptitle,...
                varargin{:},reax{:});
            
            % use this axis
            reax={'ax' ax};
            
            % add window limit markers
            % - yellow/blue dashed == noise window
            % - red/green == signal window
            span=ylim(ax);
            hold(ax,'on');
            gh(1)=plot(ax,[s.noisewin(1) s.noisewin(1)],span,'--y',...
                'linewidth',4);
            gh(2)=plot(ax,[s.noisewin(2) s.noisewin(2)],span,'--b',...
                'linewidth',4);
            gh(3)=plot(ax,[s.signalwin(1) s.signalwin(1)],span,'g',...
                'linewidth',4);
            gh(4)=plot(ax,[s.signalwin(2) s.signalwin(2)],span,'r',...
                'linewidth',4);
            hold(ax,'off');
            newplot=false;
        end
        
        % get user choice
        choice=menu(prompt,...
            ['SELECT SNR METHOD (' upper(s.method) ')'],...
            ['ADJUST NOISE WINDOW ([' num2str(s.noisewin) '])'],...
            ['ADJUST SIGNAL WINDOW ([' num2str(s.signalwin) '])'],...
            ['SELECT PLOT TYPE (' upper(func2str(s.plottype)) ')'],...
            'RETRIEVE SNR WITH THESE SETTINGS');
        
        % proceed by user choice
        switch choice
            case 1 % change snr method
                choice=menu('SELECT SNR ESTIMATION METHOD',...
                    ['CURRENT (' upper(s.method) ')'],...
                    snrmethods{:});
                
                % set method
                switch choice
                    case {2 3 4 5 6}
                        s.method=snrmethods{choice-1};
                end
                
                % update title
                if(ishandle(ax))
                    % make title
                    ptitle={['SNR Estimation Method: ' upper(s.method)]
                        'Yellow Dashed Line --  Noise Window Start'
                        '  Blue Dashed Line --  Noise Window End  '
                        '  Green Solid Line -- Signal Window Start'
                        '    Red Solid Line -- Signal Window End  '};
                    
                    % update title
                    set(get(ax,'Title'),'string',ptitle);
                end
            case 2 % change noise window
                % loop until user finalizes markers
                final=false;
                while(~final)
                    % bring plot to focus (redraw if closed)
                    if(~ishandle(ax))
                        % redraw (pretty rare to get here)
                        ax=s.plottype(data,'title',ptitle,varargin{:});
                        reax={'ax' ax};

                        % add window limit markers
                        % - yellow/blue dashed == noise window
                        % - red/green == signal window
                        span=ylim(ax);
                        hold(ax,'on');
                        gh(1)=plot(ax,[s.noisewin(1) s.noisewin(1)],span,'--y','linewidth',4);
                        gh(2)=plot(ax,[s.noisewin(2) s.noisewin(2)],span,'--b','linewidth',4);
                        gh(3)=plot(ax,[s.signalwin(1) s.signalwin(1)],span,'g','linewidth',4);
                        gh(4)=plot(ax,[s.signalwin(2) s.signalwin(2)],span,'r','linewidth',4);
                        hold(ax,'off');
                    else
                        % bring figure to focus
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
                            % left click - update noise window start
                            s.noisewin(1)=x;
                            set(gh(1),'xdata',[x x])
                        case 3
                            % right click - update noise window end
                            s.noisewin(2)=x;
                            set(gh(2),'xdata',[x x])
                        case 2
                            % middle click - finalize markers
                            if (s.noisewin(1)>s.noisewin(2))
                                % start and end reversed - fix
                                s.noisewin=s.noisewin([2 1]);
                                set(gh(1),'xdata',[s.noisewin(1) s.noisewin(1)])
                                set(gh(2),'xdata',[s.noisewin(2) s.noisewin(2)])
                            end
                            final=true;
                        otherwise
                            key2zoompan(button,ax);
                    end
                end
            case 3 % change signal window
                % loop until user finalizes markers
                final=false;
                while(~final)
                    % bring plot to focus (redraw if closed)
                    if(~ishandle(ax))
                        % redraw (pretty rare to get here)
                        ax=s.plottype(data,'title',ptitle,varargin{:});
                        reax={'ax' ax};

                        % add window limit markers
                        % - yellow/blue dashed == noise window
                        % - red/green == signal window
                        span=ylim(ax);
                        hold(ax,'on');
                        gh(1)=plot(ax,[s.noisewin(1) s.noisewin(1)],span,'--y','linewidth',4);
                        gh(2)=plot(ax,[s.noisewin(2) s.noisewin(2)],span,'--b','linewidth',4);
                        gh(3)=plot(ax,[s.signalwin(1) s.signalwin(1)],span,'g','linewidth',4);
                        gh(4)=plot(ax,[s.signalwin(2) s.signalwin(2)],span,'r','linewidth',4);
                        hold(ax,'off');
                    else
                        % bring figure to focus
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
                            % left click - update noise window start
                            s.signalwin(1)=x;
                            set(gh(3),'xdata',[x x])
                        case 3
                            % right click - update noise window end
                            s.signalwin(2)=x;
                            set(gh(4),'xdata',[x x])
                        case 2
                            % middle click - finalize markers
                            if (s.signalwin(1)>s.signalwin(2))
                                % start and end reversed - fix
                                s.signalwin=s.signalwin([2 1]);
                                set(gh(3),'xdata',[s.signalwin(1) s.signalwin(1)])
                                set(gh(4),'xdata',[s.signalwin(2) s.signalwin(2)])
                            end
                            final=true;
                        otherwise
                            key2zoompan(button,ax);
                    end
                end
            case 4 % change plot type
                % at some point we should allow custom
                choice=menu('SELECT PLOT TYPE',...
                    ['CURRENT (' upper(func2str(s.plottype)) ')'],...
                    'OVERLAY PLOT (PLOT2)',...
                    'EVENLY SPACED PLOT (PLOT0)',...
                    'DISTANCE SPACED PLOT (RECORDSECTION)');
                
                % set plot type
                switch choice
                    case 2
                        s.plottype=@plot2;
                        newplot=true;
                    case 3
                        s.plottype=@plot0;
                        newplot=true;
                    case 4
                        s.plottype=@recordsection;
                        newplot=true;
                end
                
                % handle disappearing axes
                if(ishandle(ax))
                    reax={'ax' ax};
                else
                    reax={};
                    ax=-1;
                end
            case 5 % get snr
                happy_user=true;
        end
    end
    
    % get snr
    snr=quicksnr(data,s.noisewin,s.signalwin,s.method);
    
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
