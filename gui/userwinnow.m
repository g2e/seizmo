function [data,dist,ax]=userwinnow(data,varargin)
%USERWINNOW    Interactively winnow SEIZMO records by distance
%
%    Usage:    data=userwinnow(data)
%              data=userwinnow(data,'field',value,...)
%              [data,dist,ax]=userwinnow(...)
%
%    Description:
%     DATA=USERWINNOW(DATA) presents an interactive menu and plot to
%     facilitate data winnowing by degree distance with a few mouse clicks.
%
%     DATA=USERWINNOW(DATA,'FIELD',VALUE,...) passes field/value pairs to
%     the plotting function, to allow further customization.
%
%     [DATA,DIST,AX]=USERWINNOW(...) returns a struct DIST with the
%     following fields:
%       DIST.limits  --  limits of the window applied as [START END]
%       DIST.cut     --  records winnowed out of the dataset
%     Note that the .limits field will be an empty array if no winnowing is
%     performed.
%
%    Notes:
%
%    Examples:
%     % Winnow the data as a way of pre-selecting the data before
%     % performing subsequent operations:
%     data=userwinnow(data);
%     ... some other operations ...
%
%    See also: USERWINDOW, SELECTRECORDS, GETHEADER

%     Version History:
%        Nov.  4, 2010 - initial version
%        Nov.  5, 2010 - cut field added to DIST struct
%        Nov.  9, 2010 - minor bug fixes
%        Jan.  6, 2011 - use key2zoompan
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  6, 2011 at 11:00 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>1 && ~mod(nargin,2))
    error('seizmo:userwinnow:badNumInputs',...
        'Bad number of arguments!');
end

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
    
    % rethrow error
    error(lasterror);
end

% attempt distance winnowing
try
    % winnow parameters
    dist.limits=[];
    dist.cut=[];
    
    % get distance
    gcarc=getheader(data,'gcarc');

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
        choice=menu(prompt,'DISTANCE WINNOW','DO NOT DISTANCE WINNOW');

        % proceed by user choice
        switch choice
            case 0 % closed menu
                continue;
            case 1 % distance spaced
                ax=recordsection(data,varargin{:},reax{:});
            case 2 % no winnow
                dist.limits=[];
                return;
        end
        
        % use this axis
        reax={'ax' ax};

        % add winnow limit markers
        span=xlim(ax);
        if(isempty(dist.limits)); dist.limits=ylim(ax); end
        hold(ax,'on');
        goh(1)=plot(ax,span,[dist.limits(1) dist.limits(1)],'g',...
            'linewidth',4);
        goh(2)=plot(ax,span,[dist.limits(2) dist.limits(2)],'r',...
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
                dist.limits=ylim(ax);
                hold(ax,'on');
                goh(1)=plot(ax,span,[dist.limits(1) dist.limits(1)],'g',...
                    'linewidth',4);
                goh(2)=plot(ax,span,[dist.limits(2) dist.limits(2)],'r',...
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
                    dist.limits(1)=y;
                    set(goh(1),'ydata',[y y])
                case 3
                    % right click - update winnow end
                    dist.limits(2)=y;
                    set(goh(2),'ydata',[y y])
                case 2
                    % middle click - finalize markers
                    if (dist.limits(1)>dist.limits(2))
                        % start and end reversed - fix
                        dist.limits=dist.limits([2 1]);
                        set(goh(1),'ydata',[dist.limits(1) dist.limits(1)])
                        set(goh(2),'ydata',[dist.limits(2) dist.limits(2)])
                    end
                    final=true;
                otherwise
                    key2zoompan(button,ax);
            end
        end

        % get winnowed data
        % - note that the dist.cut starts out as a logical array of the
        %   kept records -- we change it 2 lines down
        dist.cut=gcarc>=dist.limits(1) & gcarc<=dist.limits(2);
        data2=data(dist.cut);
        dist.cut=find(~dist.cut);

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
