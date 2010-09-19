function [bad,varargout]=arrcut(dd,arr,cutoff,pow,err,w,ax)
%ARRCUT    Cut arrival time outliers interactively
%
%    Usage:    bad=arrcut(dd,arr)
%              bad=arrcut(dd,arr,cutoff)
%              bad=arrcut(dd,arr,cutoff,pow)
%              bad=arrcut(dd,arr,cutoff,pow,err)
%              bad=arrcut(dd,arr,cutoff,pow,err,w)
%              bad=arrcut(dd,arr,cutoff,pow,err,w,ax)
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
%     BAD=ARRCUT(DD,ARR,CUTOFF,POW,ERR,W,AX) draws the plot in the axes
%     given by AX.
%
%     [BAD,CUTOFF]=ARRCUT(...) also returns the final outlier cutoff in
%     seconds from the best fit.
%
%     [BAD,CUTOFF,AX]=ARRCUT(...) returns the axes of the plot.
%
%    Notes:
%
%    Examples:
%     % Fit a sin curve with a 4th power polynomial and remove outliers:
%     dd=5*rand(1,100);
%     arr=sin(x)+randn(1,100);
%     arrcut(dd,arr,2,4,rand(1,100));
%
%    See also: AMPCUT, ERRCUT, SNRCUT, POPCUT

%     Version History:
%        Sep. 17, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 17, 2010 at 20:00 GMT

% todo:

% check nargin
error(nargchk(2,7,nargin));

% defaults
if(nargin<3 || isempty(cutoff)); cutoff=nan; end
if(nargin<4 || isempty(pow)); pow=1; end
if(nargin<5); err=[]; end
if(nargin<6); w=[]; end
if(nargin<7); ax=[]; end

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
if(~isempty(ax) && (~isscalar(ax) || ~isreal(ax)))
    error('seizmo:arrcut:badInput',...
        'AX must be a handle to a single axis!');
end

% expand scalars
if(isscalar(err)); err=expandscalars(err,arr); end
if(isscalar(w)); w=expandscalars(w,arr); end

% get the fit
m=wlinem(dd,arr,pow,diag(err),diag(w))';

% residuals
resid=arr-polyval(fliplr(m),dd);
std=sqrt(var(resid));
resid=abs(resid);

% default cutoff
if(isnan(cutoff)); cutoff=3*std; end

% check the axes
if(isempty(ax) || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    % new figure
    fh=figure('color','w');
    ax=axes('parent',fh);
else
    cla(ax,'reset');
    axes(ax);
end

% draw the fit and cutoff limits
maxdd=max(dd(:));
mindd=min(dd(:));
pdd=linspace(mindd-0.1*(maxdd-mindd),maxdd+0.1*(maxdd-mindd),100)';
parr=polyval(fliplr(m),pdd);
hfit=plot(pdd,parr,'b','linewidth',2);
hold(ax,'on');
hcut=plot(pdd,[parr+cutoff parr-cutoff],'r--','linewidth',2);

% draw the points (w/ or w/o errorbars)
if(isempty(err))
    hpnts=plot(dd,arr,'ko');
else
    h=ploterr(dd,arr,[],err,'ko');
    hpnts=h(1);
    hy=h(2);
end
hold(ax,'off');

% label plot
set(ax,'fontsize',10,'fontweight','bold');
xlabel('Distance (^o)','fontsize',10,'fontweight','bold');
ylabel('Arrival Time (s)','fontsize',10,'fontweight','bold');
title({'Left Click = Change Cutoff';
    'Middle Click = Implement Cut';
    ['Cutoff = ' num2str(cutoff) 's (' num2str(cutoff/std) ' stddev)'];
    '';
    ['t_{arr} = ' polystr(fliplr(m),'\Delta')]},...
    'fontsize',10,'fontweight','bold');

% let user adjust the limits
unhappy=true;
while(unhappy)
    axis(ax);
    [x,y,button]=ginput(1);
    switch button
        case 1
            % get cutoff
            cutoff=abs(y-polyval(fliplr(m),x));
            
            % adjust limits
            set(hcut(1),'ydata',parr+cutoff);
            set(hcut(2),'ydata',parr-cutoff);
            
            % reset title
            set(get(ax,'title'),'string',...
                {'Left Click = Change Cutoff';
                'Middle Click = Implement Cut';
                ['Cutoff = ' num2str(cutoff) 's (' ...
                num2str(cutoff/std) ' stddev)'];
                '';
                ['t_{arr} = ' polystr(fliplr(m),'\Delta')]});
        case 2
            bad=find(resid>cutoff);
            unhappy=false;
    end
end

% output if desired
if(nargout>1); varargout={cutoff ax}; end

end
