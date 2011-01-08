function [bad,varargout]=snrcut(dd,snr,cutoff,ax)
%SNRCUT    Cut low SNR values interactively
%
%    Usage:    bad=snrcut(dd,err)
%              bad=snrcut(dd,err,cutoff)
%              bad=snrcut(dd,err,cutoff,ax)
%              [bad,cutoff]=snrcut(...)
%              [bad,cutoff,ax]=snrcut(...)
%
%    Description:
%     BAD=SNRCUT(DD,SNR) creates an interactive plot of signal to noise
%     ratio SNR versus degree distance DD where the user is expected to
%     indicate the cutoff for removing low values.  The default cutoff is
%     2.5 and is altered when the user left clicks the plot.  To complete
%     the cut and exit, the user must middle-click the plot.  The indices
%     of the cut are returned in BAD.
%
%     BAD=SNRCUT(DD,SNR,CUTOFF) adjusts the default cutoff to CUTOFF.  The
%     default SNR cutoff is 2.5.
%
%     BAD=SNRCUT(DD,SNR,CUTOFF,AX) draws the plot in the axes given by AX.
%
%     [BAD,CUTOFF]=SNRCUT(...) returns the final cutoff used for removing
%     the outliers.
%
%     [BAD,CUTOFF,AX]=SNRCUT(...) returns the axes of the plot
%
%    Notes:
%
%    Examples:
%     % Trim some random values:
%     dd=5*rand(1,100);
%     snrcut(dd,abs(randn(1,100)));
%
%    See also: ARRCUT, AMPCUT, ERRCUT, POPCUT

%     Version History:
%        Sep. 17, 2010 - initial version
%        Dec. 12, 2010 - fixed several plotting bugs
%        Jan.  6, 2011 - proper ginput handling, use key2zoompan
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  6, 2011 at 23:55 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% defaults
if(nargin<3 || isempty(cutoff)); cutoff=2.5; end
if(nargin<4); ax=[]; end

% check inputs
if(~isreal(dd) || ~isreal(snr) || isempty(dd) ...
        || ~isequal(size(dd),size(snr)))
    error('seizmo:snrcut:badInput',...
        'DD & SNR must be equal sized nonempty real valued arrays!');
end
if(~isreal(cutoff) || ~isscalar(cutoff) || cutoff<0)
    error('seizmo:snrcut:badInput',...
        'CUTOFF must be a positive value!');
end
if(~isempty(ax) && (~isscalar(ax) || ~isreal(ax)))
    error('seizmo:snrcut:badInput',...
        'AX must be a handle to a single axis!');
end

% check the axes
if(isempty(ax) || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    % new figure
    fh=figure('color','w');
    ax=axes('parent',fh);
else
    cla(ax,'reset');
    axes(ax);
end

% make stem plot
hpts=stem(ax,dd,snr,'fill');
set(hpts,'markerfacecolor','g','linewidth',1,'color','k');

% add cutoff bar
hold(ax,'on');
xlimits=get(ax,'xlim');
hcut=plot(ax,xlimits,[cutoff cutoff],'r--','linewidth',2);
hold(ax,'off');

% label plot
set(ax,'fontsize',10,'fontweight','bold');
xlabel(ax,'Distance (^o)','fontsize',10,'fontweight','bold');
ylabel(ax,'SNR','fontsize',10,'fontweight','bold');
title(ax,{'Left Click = Change Cutoff';
    'Middle Click = Implement Cut';
    ['Cutoff = ' num2str(cutoff)]},...
    'fontsize',10,'fontweight','bold');

% let user adjust the limits
unhappy=true;
while(unhappy)
    % get snr cutoff
    axis(ax);
    try
        [x,y,button]=ginput(1);
    catch
        % plot closed - break out
        break;
    end
    
    % skip if not correct axis
    if(ax~=gca); continue; end
    
    % act based on button
    switch button
        case 1
            % get cutoff
            cutoff=abs(y);
            
            % adjust limits
            set(hcut,'ydata',[cutoff cutoff]);
            
            % reset title
            set(get(ax,'title'),'string',...
                {'Left Click = Change Cutoff';
                'Middle Click = Implement Cut';
                ['Cutoff = ' num2str(cutoff)]});
        case 2
            unhappy=false;
        otherwise
            key2zoompan(button,ax);
    end
end

% find bad
bad=find(snr>cutoff);

% output if desired
if(nargout>1); varargout={cutoff ax}; end

end
