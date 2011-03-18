function [bad,varargout]=arrcut(dd,arr,cutoff,pow,err,w,color,ax)
%ARRCUT    Cut arrival time outliers interactively
%
%    Usage:    bad=arrcut(dd,arr)
%              bad=arrcut(dd,arr,cutoff)
%              bad=arrcut(dd,arr,cutoff,pow)
%              bad=arrcut(dd,arr,cutoff,pow,err)
%              bad=arrcut(dd,arr,cutoff,pow,err,w)
%              bad=arrcut(dd,arr,cutoff,pow,err,w,color)
%              bad=arrcut(dd,arr,cutoff,pow,err,w,color,ax)
%              [bad,cutoff]=arrcut(...)
%              [bad,cutoff,ax]=arrcut(...)
%
%    Description:
%     BAD=ARRCUT(DD,ARR) creates an interactive plot of arrival times ARR
%     versus degree distance DD where the user is expected to indicate the
%     limits from the best fitting trend for removing outliers.  The
%     default cutoff is 3 standard deviations from the best fitting line
%     and is altered when the user left clicks on the plot.  The user can
%     perform the cut and exit by middle clicking the plot.  The trend by
%     default is linear (see the following usage form to adjust this).  The
%     indices of the cut arrivals are returned in BAD.
%
%     BAD=ARRCUT(DD,ARR,CUTOFF) adjusts the default cutoff to CUTOFF in
%     seconds.  The default is 3 standard deviations from the best fit line
%     (note this is NOT in seconds).
%
%     BAD=ARRCUT(DD,ARR,CUTOFF,POW) alters the power of the polynomial fit.
%     The default value is 1.
%
%     BAD=ARRCUT(DD,ARR,CUTOFF,POW,ERR) includes the arrival time errors in
%     ERR by drawing whiskers around the points.  This does not effect the
%     fit (use weighting to do that).
%
%     BAD=ARRCUT(DD,ARR,CUTOFF,POW,ERR,W) uses the weights given in W for
%     determining the fit.  W should be the same size as ARR.
%
%     BAD=ARRCUT(DD,ARR,CUTOFF,POW,ERR,W,COLOR) sets the facecolor of the
%     points in the plots.  COLOR may be a color name, a rgb triplet, or an
%     Nx3 array of triplets where N is the number of points.  The default
%     is 'w'.
%
%     BAD=ARRCUT(DD,ARR,CUTOFF,POW,ERR,W,COLOR,AX) draws the plot in the
%     axes given by AX.  AX should be 2 axes handles (the 1st is dist vs
%     time, the 2nd is dist vs residual).
%
%     [BAD,CUTOFF]=ARRCUT(...) also returns the final outlier cutoff in
%     seconds from the best fit.
%
%     [BAD,CUTOFF,AX]=ARRCUT(...) returns the axes of the 2 plots.
%
%    Notes:
%
%    Examples:
%     % Fit a sin curve with a 4th power polynomial and remove outliers:
%     dd=5*rand(1,100);
%     arr=sin(dd)+randn(1,100);
%     arrcut(dd,arr,2,4,rand(1,100));
%
%    See also: AMPCUT, ERRCUT, SNRCUT, POPCUT

%     Version History:
%        Sep. 17, 2010 - initial version
%        Dec. 12, 2010 - fixed several plotting bugs, no error input to
%                        WLINEM (improperly used anyway)
%        Jan.  6, 2011 - proper ginput handling, use key2zoompan
%        Jan.  7, 2011 - using errorbar now
%        Jan. 26, 2011 - draw small residual plot below dist vs time plot
%        Mar.  6, 2011 - coloring of marker faces
%        Mar. 10, 2011 - fix color issues
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 10, 2011 at 23:55 GMT

% todo:

% check nargin
error(nargchk(2,8,nargin));

% defaults
if(nargin<3 || isempty(cutoff)); cutoff=nan; end
if(nargin<4 || isempty(pow)); pow=1; end
if(nargin<5); err=[]; end
if(nargin<6); w=[]; end
if(nargin<7 || isempty(color)); color='w'; end
if(nargin<8); ax=[]; end

% check inputs
if(~isreal(dd) || ~isreal(arr) || isempty(dd) ...
        || ~isequal(size(dd),size(arr)))
    error('seizmo:arrcut:badInput',...
        'DD & ARR must be equal sized nonempty real valued arrays!');
end
if(~isreal(cutoff) || ~isscalar(cutoff) || cutoff<0)
    error('seizmo:arrcut:badInput',...
        'CUTOFF must be a positive value in seconds!');
end
if(~isreal(pow) || ~isscalar(pow) || pow<0 || pow~=fix(pow))
    error('seizmo:arrcut:badInput',...
        'POW must be a positive integer >=0');
end
if(~isempty(err) && (~isreal(err) || ~isequalsizeorscalar(arr,err)))
    error('seizmo:arrcut:badInput',...
        'ERR & ARR must be equal sized real valued arrays!');
end
if(~isempty(w) && (~isreal(w) || ~isequalsizeorscalar(arr,w)))
    error('seizmo:arrcut:badInput',...
        'W & ARR must be equal sized real valued arrays!');
end
if(ischar(color))
    % keep 'none' or try name2rgb (it errors if not valid)
    if(strcmpi(color,'none'))
        color='w';
    else
        color=name2rgb(color);
    end
elseif(~isreal(color) || ndims(color)~=2 ...
        || any(color(:)<0 | color(:)>1) || size(color,2)~=3 ...
        || ~any(size(color,1)~=[1 numel(dd)]))
    error('seizmo:arrcut:badInput',...
        'Numeric COLOR must be a valid rgb triplet!');
end
if(~isempty(ax) && (numel(ax)~=2 || ~isreal(ax)))
    error('seizmo:arrcut:badInput',...
        'AX must be a handle to 2 axes!');
end

% expand scalars
if(isscalar(err)); err=expandscalars(err,arr); end
if(isscalar(w)); w=expandscalars(w,arr); end

% get the fit
m=wlinem(dd,arr,pow,[],diag(w))';

% residuals
resid=arr-polyval(fliplr(m),dd);
std=sqrt(var(resid));
aresid=abs(resid);

% default cutoff
if(isnan(cutoff)); cutoff=3*std; end

% check the axes
if(isempty(ax) || any(~ishandle(ax)) ...
        || any(~strcmp('axes',get(ax,'type'))))
    % new figure
    fh=figure('color','w');
    ax(1)=subplot(5,1,2:3,'parent',fh);
    ax(2)=subplot(5,1,4:5,'parent',fh);
else
    cla(ax(1),'reset');
    cla(ax(2),'reset');
    axes(ax(1));
end

% draw the fit and cutoff limits
maxdd=max(dd(:));
mindd=min(dd(:));
pdd=[mindd-0.1*(maxdd-mindd); maxdd+0.1*(maxdd-mindd)];
parr=polyval(fliplr(m),pdd);
hfit=plot(ax(1),pdd,parr,'b','linewidth',2);
set(hfit,'tag','fit');
hold(ax(1),'on');
hcut1=plot(ax(1),pdd,[parr+cutoff parr-cutoff],'r--','linewidth',2);
hcut2=plot(ax(2),pdd,[0*parr+cutoff 0*parr-cutoff],'r--','linewidth',2);
set(hcut1,'tag','cut');
set(hcut2,'tag','cut');
linkaxes(ax,'x');
hold(ax(2),'on');

% draw the points (w/ or w/o errorbars)
if(isempty(err))
    hpnts1=scatter(ax(1),dd,arr,[],color,'filled','markeredgecolor','k');
    drawnow;
    set(hpnts1,'tag','points');
    hpnts2=scatter(ax(2),dd,resid,[],color,'filled','markeredgecolor','k');
    drawnow;
    set(hpnts2,'tag','points');
else
    h1=errorbar(ax(1),dd,arr,err,'k.');
    set(h1,'tag','errorbars');
    hpnts1=scatter(ax(1),dd,arr,[],color,'filled','markeredgecolor','k');
    drawnow;
    set(hpnts1,'tag','points');
    h2=errorbar(ax(2),dd,resid,err,'k.');
    set(h2,'tag','errorbars');
    hpnts2=scatter(ax(2),dd,resid,[],color,'filled','markeredgecolor','k');
    drawnow;
    set(hpnts2,'tag','points');
end
hold(ax(1),'off');
hold(ax(2),'off');

% label plot
set(ax(1),'fontsize',10,'fontweight','bold','xticklabel',[]);
set(ax(2),'fontsize',10,'fontweight','bold');
xlabel(ax(2),'Distance (^o)','fontsize',10,'fontweight','bold');
ylabel(ax(1),'Arrival Time (s)','fontsize',10,'fontweight','bold');
ylabel(ax(2),'Residual (s)','fontsize',10,'fontweight','bold');
title(ax(1),{'Left Click = Change Cutoff';
    'Middle Click = Implement Cut';
    ['Cutoff = ' num2str(cutoff) 's (' num2str(cutoff/std) ' stddev)'];
    '';
    ['t_{arr} = ' polystr(fliplr(m),'\Delta')]},...
    'fontsize',10,'fontweight','bold');

% let user adjust the limits
unhappy=true;
while(unhappy)
    % get arrival time outlier limits
    axis(ax(1));
    try
        [x,y,button]=ginput(1);
    catch
        % plot closed - break out
        break;
    end
    
    % skip if not correct axis
    if(~any(gca==ax)); continue; end
    idx=find(gca==ax);
    
    % act based on button
    switch button
        case 1
            % get cutoff
            switch idx
                case 1
                    cutoff=abs(y-polyval(fliplr(m),x));
                case 2
                    cutoff=abs(y);
            end
            
            % adjust limits
            set(hcut1(1),'ydata',parr+cutoff);
            set(hcut1(2),'ydata',parr-cutoff);
            set(hcut2(1),'ydata',0*parr+cutoff);
            set(hcut2(2),'ydata',0*parr-cutoff);
            
            % reset title
            set(get(ax(1),'title'),'string',...
                {'Left Click = Change Cutoff';
                'Middle Click = Implement Cut';
                ['Cutoff = ' num2str(cutoff) 's (' ...
                num2str(cutoff/std) ' stddev)'];
                '';
                ['t_{arr} = ' polystr(fliplr(m),'\Delta')]});
        case 2
            unhappy=false;
        otherwise
            key2zoompan(button,ax(idx));
    end
end

% find bad
bad=find(aresid>cutoff);

% output if desired
if(nargout>1); varargout={cutoff ax}; end

end
