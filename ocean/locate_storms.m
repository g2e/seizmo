function [d,t0,bt]=locate_storms(t,T,varargin)
%LOCATE_STORMS    Interactively determine swell source (storms) info
%
%    Usage:    [x,t0]=locate_storms(t,T)
%              [x,t0]=locate_storms(t,T,...)
%
%    Description:
%     [x,t0]=LOCATE_STORMS(t,T) plots the input period T vs time t and lets
%     the user go through the data and pick limits for the data to be used
%     in estimating a particular storm's distance and origin time.  The
%     user left clicks to select the boundaries, right clicks to undo, and
%     middle clicks to exit.  Once the boundaries are selected a prediction
%     is plotted for that storm and the process may be repeated.
%
%     [x,t0]=LOCATE_STORMS(t,T,...) passes on any additional inputs to PLOT
%     when the input data is being drawn.
%
%    Notes:
%
%    Examples:
%     % Read in some WW3 primary wave period hindcast data at a specific
%     % location and then plot it up and guide the user to find storms:
%     s=ww3struct('nww3.tp.200608.grb',[],[],[],[0 0],[0 1]);
%     [x,t0]=locate_storms(s.time,s.data{1});
%
%    See also: SWELL_BACKPROJ, SWELL_FORWARD, PLOTWW3TS

%     Version History:
%        May  11, 2012 - initial version
%        May  31, 2012 - better checking
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  31, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% require equal number of elements
if(~isnumeric(t) || ~isnumeric(T))
    error('seizmo:locate_storms:badInput',...
        't & T must be numeric!');
elseif(numel(t)~=numel(T) || numel(t)<2 || ~isvector(t))
    error('seizmo:locate_storms:badInput',...
        't & T must be equal-sized vectors and have at least 2 elements!');
elseif(any(diff(t)<0))
    error('seizmo:locate_storms:badInput',...
        't must be monotonically increasing!');
end

% create plot
fh=figure('color','w');
ax=axes('parent',fh);
h=plot(ax,t(:),T(:),varargin{:});

% line color of boundaries
bcolor=[1 0 0];
blw=2;
bplace='back';

% line color of prediction & where it is drawn
pcolor=[.5 .5 .5];
plw=2;
pplace='back';

% output
d=[]; t0=[];

% period limits for boundary extent
Tlim=[min(T(:)) max(T(:))];

% make selections
none=true; cnt=0; blh=[];
while(true)
    % get user input
    try
        % force axes to front first
        axes(ax);
        [x,y,button]=ginput(1);
    catch
        % plot closed so finish
        return;
    end
    
    % require ax
    if(gca~=ax); continue; end
    
    % act by button
    switch button
        case 1 % left or ctrl+shift
            % 1st or 2nd boundary?
            if(none) % 1st = draw line
                % increment
                cnt=cnt+1;
                none=false;
                
                % save t
                bt(cnt,1)=x;
                
                % draw line
                hold(ax,'on');
                blh(cnt,1)=plot(ax,[x x],Tlim,'-',...
                    'color',bcolor,'linewidth',blw);
                movekids(blh(cnt,1),bplace);
                hold(ax,'off');
            else % 2nd = draw line, get storm info, draw prediction
                % save t
                bt(cnt,2)=x;
                
                % draw line
                hold(ax,'on');
                blh(cnt,2)=plot(ax,[x x],Tlim,'-',...
                    'color',bcolor,'linewidth',blw);
                movekids(blh(cnt,2),bplace);
                hold(ax,'off');
                
                % get storm info
                bts=sort(bt(cnt,:));
                [d(cnt,1),t0(cnt,1)]=swell_backproj(...
                    T(t>=bts(1) & t<=bts(2)),t(t>=bts(1) & t<=bts(2)));
                
                % print summary
                fprintf(['\nSTORM #' num2str(cnt) ' VALUES:\n']);
                fprintf('Distance:  %gm\n',d(cnt,1));
                fprintf('Time: %gs\n',t0(cnt,1));
                
                % make prediction and draw it
                pT=5:55;
                pt=swell_forward(pT,d(cnt,1))+t0(cnt,1);
                hold(ax,'on');
                plh(cnt)=plot(ax,pt,pT,'-','color',pcolor,'linewidth',plw);
                movekids(plh(cnt),pplace);
                hold(ax,'off');
                
                % no boundaries
                none=true;
            end
        case 2 % middle or shift
            return;
        case 3 % right or ctrl
            % undo last action
            if(~cnt) % skip if nothing done
                continue;
            elseif(none) % delete last storm & boundary
                % delete prediction & last boundary
                delete([plh(cnt) blh(cnt,2)]);
                
                % remove storm info
                d(cnt)=[];
                t0(cnt)=[];
                
                % there is 1 pnt
                none=false;
            else % delete line & decrement
                delete(blh(cnt,1));
                blh(cnt,:)=[];
                bt(cnt,:)=[];
                cnt=cnt-1;
                none=true;
            end
        otherwise
            key2zoompan(button,ax);
    end
end

end
