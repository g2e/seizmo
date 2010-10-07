function [bad,varargout]=ampcut(dd,amp,cutoff,pow,err,w,ax)
%AMPCUT    Cut amplitude outliers interactively
%
%    Usage:    bad=ampcut(dd,amp)
%              bad=ampcut(dd,amp,cutoff)
%              bad=ampcut(dd,amp,cutoff,pow)
%              bad=ampcut(dd,amp,cutoff,pow,err)
%              bad=ampcut(dd,amp,cutoff,pow,err,w)
%              bad=ampcut(dd,amp,cutoff,pow,err,w,ax)
%              [bad,cutoff]=ampcut(...)
%              [bad,cutoff,ax]=ampcut(...)
%
%    Description:
%     BAD=AMPCUT(DD,AMP) creates an interactive plot of amplitudes AMP
%     versus degree distance DD where the user is expected to indicate the
%     limits from the best fitting trend for removing outliers.  The
%     default cutoff is 3 standard deviations from the best fitting line
%     and is altered when the user left clicks on the plot.  The user can
%     perform the cut and exit by middle clicking the plot.  The trend by
%     default is linear (see the following usage form to adjust this).  The
%     indices of the cut amplitudes are returned in BAD.
%
%     BAD=AMPCUT(DD,AMP,CUTOFF) adjusts the default cutoff to CUTOFF in
%     standard deviations.  The default is 3 standard deviations from the
%     best fit line.
%
%     BAD=AMPCUT(DD,AMP,CUTOFF,POW) alters the power of the polynomial fit.
%     The default value is 1.
%
%     BAD=AMPCUT(DD,AMP,CUTOFF,POW,ERR) includes the amplitude errors in
%     ERR by drawing whiskers around the points.  This does not effect the
%     fit (use weighting to do that).
%
%     BAD=AMPCUT(DD,AMP,CUTOFF,POW,ERR,W) uses the weights given in W for
%     determining the fit.  W should be the same size as AMP.
%
%     BAD=AMPCUT(DD,AMP,CUTOFF,POW,ERR,W,AX) draws the plot in the axes
%     given by AX.
%
%     [BAD,CUTOFF]=AMPCUT(...) also returns the final outlier cutoff in
%     standard deviations from the best fit.
%
%     [BAD,CUTOFF,AX]=AMPCUT(...) returns the axes of the plot.
%
%    Notes:
%
%    Examples:
%     % Fit to some simple exponentially decaying amplitudes:
%     dd=5*rand(1,100);
%     amp=exp(-dd).*(1+randn(1,100)/10);
%     ampcut(dd,amp,[],1,amp.*rand(1,100)/2);
%
%    See also: ARRCUT, ERRCUT, SNRCUT, POPCUT

%     Version History:
%        Sep. 17, 2010 - initial version
%        Sep. 28, 2010 - natural log not log base 10
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 28, 2010 at 20:00 GMT

% todo:

% check nargin
error(nargchk(2,7,nargin));

% defaults
if(nargin<3 || isempty(cutoff)); cutoff=3; end
if(nargin<4 || isempty(pow)); pow=1; end
if(nargin<5); err=[]; end
if(nargin<6); w=[]; end
if(nargin<7); ax=[]; end

% check inputs
if(~isreal(dd) || ~isreal(amp) || isempty(dd) ...
        || ~isequal(size(dd),size(amp)))
    error('seizmo:ampcut:badInput',...
        'DD & AMP must be equal sized nonempty real valued arrays!');
end
if(~isreal(cutoff) || ~isscalar(cutoff) || cutoff<0)
    error('seizmo:ampcut:badInput',...
        'CUTOFF must be a positive value in standard deviations!');
end
if(~isreal(pow) || ~isscalar(pow) || pow<0 || pow~=fix(pow))
    error('seizmo:ampcut:badInput',...
        'POW must be a positive integer >=0');
end
if(~isempty(err) && (~isreal(err) || ~isequalsizeorscalar(amp,err)))
    error('seizmo:ampcut:badInput',...
        'ERR & AMP must be equal sized real valued arrays!');
end
if(~isempty(w) && (~isreal(w) || ~isequalsizeorscalar(amp,w)))
    error('seizmo:ampcut:badInput',...
        'W & AMP must be equal sized real valued arrays!');
end
if(~isempty(ax) && (~isscalar(ax) || ~isreal(ax)))
    error('seizmo:ampcut:badInput',...
        'AX must be a handle to a single axis!');
end

% expand scalars
if(isscalar(err)); err=expandscalars(err,amp); end
if(isscalar(w)); w=expandscalars(w,amp); end

% get the fit
m=wlinem(dd,log(amp),pow,diag(log(amp)-log(amp+err)),diag(w))';

% residuals
resid=log(amp)-polyval(fliplr(m),dd);
std=sqrt(var(resid));
resid=abs(resid);

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
pamp=polyval(fliplr(m),pdd);
hfit=plot(pdd,pamp,'b','linewidth',2);
hold(ax,'on');
hcut=plot(pdd,[pamp+cutoff*std pamp-cutoff*std],'r--','linewidth',2);

% draw the points (w/ or w/o errorbars)
if(isempty(err))
    hpnts=plot(dd,log(amp),'ko');
else
    h=ploterr(dd,log(amp),[],{log(amp-err) log(amp+err)},'ko');
    hpnts=h(1);
    hy=h(2);
end
hold(ax,'off');

% label plot
set(ax,'fontsize',10,'fontweight','bold');
xlabel('Distance (^o)','fontsize',10,'fontweight','bold');
ylabel('ln(A)','fontsize',10,'fontweight','bold');
title({'Left Click = Change Cutoff';
    'Middle Click = Implement Cut';
    ['Cutoff = ' num2str(cutoff) ' stddev'];
    '';
    ['ln(A) = ' polystr(fliplr(m),'\Delta')]},...
    'fontsize',10,'fontweight','bold');

% let user adjust the limits
unhappy=true;
while(unhappy)
    axis(ax);
    [x,y,button]=ginput(1);
    switch button
        case 1
            % get cutoff
            cutoff=abs(y-polyval(fliplr(m),x))/std;
            
            % adjust limits
            set(hcut(1),'ydata',pamp+cutoff*std);
            set(hcut(2),'ydata',pamp-cutoff*std);
            
            % reset title
            set(get(ax,'title'),'string',...
                {'Left Click = Change Cutoff';
                'Middle Click = Implement Cut';
                ['Cutoff = ' num2str(cutoff) ' stddev'];
                '';
                ['ln(A) = ' polystr(fliplr(m),'\Delta')]});
        case 2
            bad=find(resid>cutoff*std);
            unhappy=false;
    end
end

% output if desired
if(nargout>1); varargout={cutoff ax}; end

end
