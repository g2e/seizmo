function [data,win,ax]=userwinnow(data,limits,varargin)
%USERWINNOW    Interactively winnow SEIZMO records by distance
%
%    Usage:    data=userwinnow(data)
%              data=userwinnow(data,limits)
%              data=userwinnow(data,limits,'field',value,...)
%              [data,win,ax]=userwinnow(...)
%
%    Description:
%     DATA=USERWINNOW(DATA) presents an interactive menu and plot to
%     facilitate data winnowing by degree distance with a few mouse clicks.
%     The green bar indicates the lower limit while the red bar indicates
%     the upper limit.  SWITCHING THE ORDER (RED LOWER THAN GREEN) WILL
%     EXCLUDE THE RANGE RATHER THAN INCLUDE IT, SO BE CAREFUL.
%
%     DATA=USERWINNOW(DATA,LIMITS) specifies the initial limits of the
%     winnow.  Note that LIMITS should be a 1x2 vector of [START END].  If
%     LIMITS is specified as [END START], then data within the range is
%     winnowed out by default rather than in (this is useful when winnowing
%     by azimuth and the range desired wraps the axes limits).
%
%     DATA=USERWINNOW(DATA,LIMITS,'FIELD',VALUE,...) passes field/value
%     pairs to the plotting function, to allow further customization.  In
%     particular, using the 'yfield' option allows winnowing by other
%     header fields than GCARC (such as AZ or BAZ).
%
%     [DATA,DIST,AX]=USERWINNOW(...) returns a struct WIN with the
%     following fields:
%       WIN.yfield  --  header field used in winnowing
%       WIN.limits  --  limits of the window applied
%       WIN.cut     --  records winnowed out of the dataset
%     Note that the .limits field will be an empty array if no winnowing is
%     performed.
%
%    Notes:
%     - RED above GREEN == winnow in
%       GREEN above RED == winnow out
%
%    Examples:
%     % Winnowing by azimuth:
%     data=userwinnow(data,[0 360],'yfield','az');
%
%    See also: USERWINDOW, SELECTRECORDS, GETHEADER, RECORDSECTION,
%              USERTAPER, USERMOVEOUT

%     Version History:
%        Nov.  4, 2010 - initial version
%        Nov.  5, 2010 - cut field added to DIST struct
%        Nov.  9, 2010 - minor bug fixes
%        Jan.  6, 2011 - use key2zoompan
%        Jan. 17, 2011 - initial limits arg, allow exclusion range, yfield
%                        fully encouraged now
%        Jan. 29, 2011 - fix empty yfield bug (default to gcarc)
%        Apr.  3, 2012 - use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 11:00 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>2 && mod(nargin,2))
    error('seizmo:userwinnow:badNumInputs',...
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

% attempt distance winnowing
try
    % check limits
    if(nargin<2)
        limits=[];
    elseif(~isempty(limits) && (numel(limits)~=2 || ~isreal(limits)))
        error('seizmo:userwinnow:badInput',...
            'INITWIN must be 1x2 vector as [START END]!');
    end
    
    % get yfield
    yfield='gcarc';
    if(nargin>2 && ~iscellstr(varargin(1:2:end)))
        error('seizmo:userwinnow:badInput',...
            'Plot options must be specified as ''field''/value pairs!');
    end
    for i=1:2:nargin-2
        switch lower(varargin{i})
            case 'yfield'
                if(~isempty(varargin{i+1}))
                    yfield=varargin{i+1};
                end
        end
    end
    
    % winnow parameters
    win.yfield=yfield;
    win.limits=limits;
    win.cut=[];
    
    % get header field (gcarc by default)
    yvalue=getheader(data,yfield);

    % outer loop - only breaks free by user command
    happy_user=false; ax=-1; reax={};
    while(~happy_user)
        % explain to the user how this works with a little prompt
        % and make them decide what kind of plot to use for the
        % winowing, offering them back out options.  This prompt
        % looks like trash because of default menu fonts.
        prompt={'+-------------------------------------------------------+'
                '|                Welcome to SEIZMO''s interactive winnowing function             |'
                '+-------------------------------------------------------+'
                '|                                                                                                                |'
                '|                                             MOUSE USAGE                                             |'
                '|                                                                                                                |'
                '|    LEFT CLICK                       MIDDLE CLICK                      RIGHT CLICK   |'
                '+-----------+                +------------+               +------------+'
                '|   Mark Winnow                     Finalize Marks                     Mark Winnow   |'
                '|          Start                                                                              End           |'
                '+-------------------------------------------------------+'
                '|                                                                                                                |'
                '|                                                   NOTES                                                   |'
                '|                                                                                                                |'
                '|          + You may refine winnow marks until you finalize                        |'
                '|          + When finalized, a new plot with the winnowed                           |'
                '|              waveforms will appear, as well as a confirmation                      |'
                '|              prompt.  You will have the option to re-winnow.                      |'
                '|                                                                                                                |'
                '+-------------------------------------------------------+'
                '|                                                                                                                |'
                '|                 PLEASE CHOOSE AN OPTION BELOW TO PROCEED!                 |'
                '|                                                                                                                |'
                '+-------------------------------------------------------+'};


        % way cooler menu -- if only matlab gui's used fixed width
        if(strcmpi(getapplication,'OCTAVE'))
            prompt={'+-------------------------------------------------------+'
                    '|  Welcome to SEIZMO''s interactive winnowing function   |'
                    '+-------------------------------------------------------+'
                    '|                                                       |'
                    '|                     MOUSE USAGE:                      |'
                    '|                                                       |'
                    '| LEFT CLICK          MIDDLE CLICK          RIGHT CLICK |'
                    '+------------+      +--------------+      +-------------+'
                    '| Mark Winnow        Finalize Marks         Mark Winnow |'
                    '|   Start                                      End      |'
                    '+-------------------------------------------------------+'
                    '|                                                       |'
                    '|                        NOTES:                         |'
                    '|                                                       |'
                    '|  + You may refine winnow marks until you finalize     |'
                    '|  + When finalized, a new plot with the winnowed       |'
                    '|    waveforms will appear, as well as a confirmation   |'
                    '|    prompt.  You will have the option to re-winnow.    |'
                    '|                                                       |'
                    '+-------------------------------------------------------+'
                    '|                                                       |'
                    '|        PLEASE CHOOSE AN OPTION BELOW TO PROCEED!      |'
                    '|                                                       |'
                    '+-------------------------------------------------------+'};
        end

        % display prompt and get user choice
        choice=menu(prompt,...
            ['APPLY ' upper(yfield) ' WINNOW'],...
            'DO NOT WINNOW');

        % proceed by user choice
        switch choice
            case 0 % closed menu
                continue;
            case 1 % distance spaced
                ax=recordsection(data,varargin{:},reax{:});
            case 2 % no winnow
                win.limits=[];
                return;
        end
        
        % use this axis
        reax={'ax' ax};

        % add winnow limit markers
        span=xlim(ax);
        if(isempty(win.limits)); win.limits=ylim(ax); end
        hold(ax,'on');
        goh(1)=plot(ax,span,[win.limits(1) win.limits(1)],'g',...
            'linewidth',4);
        goh(2)=plot(ax,span,[win.limits(2) win.limits(2)],'r',...
            'linewidth',4);
        hold(ax,'off');

        % loop until user finalizes markers
        final=false;
        while(~final)
            % bring plot to focus (redraw if closed)
            if(~ishandle(ax))
                % redraw (pretty rare to get here)
                ax=recordsection(data,varargin{:});
                reax={'ax' ax};
                span=xlim(ax);
                win.limits=ylim(ax);
                hold(ax,'on');
                goh(1)=plot(ax,span,[win.limits(1) win.limits(1)],'g',...
                    'linewidth',4);
                goh(2)=plot(ax,span,[win.limits(2) win.limits(2)],'r',...
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
                % user closed winnow - break from loop
                button=2;
            end

            % which mouse button?
            switch button
                case 1
                    % left click - update winnow start
                    win.limits(1)=y;
                    set(goh(1),'ydata',[y y])
                case 3
                    % right click - update winnow end
                    win.limits(2)=y;
                    set(goh(2),'ydata',[y y])
                case 2
                    % middle click - finalize markers
                    final=true;
                otherwise
                    key2zoompan(button,ax);
            end
        end

        % get winnowed data
        % - note that the win.cut starts out as a logical array of the
        %   kept records -- we change it a few lines down to those cut
        if(win.limits(2)>win.limits(1))
            % keep those in range
            win.cut=yvalue>=win.limits(1) & yvalue<=win.limits(2);
        else
            % delete those in range
            win.cut=yvalue<win.limits(1) & yvalue>win.limits(2);
        end
        data2=data(win.cut);    % data2 are the kept records
        win.cut=find(~win.cut); % switch from kept to cut now

        % proceed by user choice
        ax=recordsection(data2,varargin{:},reax{:});

        % confirm results
        choice=0;
        while(~choice)
            choice=menu('KEEP WINNOW?',...
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
