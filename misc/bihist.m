function [varargout]=bihist(varargin)
%BIHIST    Plots histograms for 2 datasets with one going up & one down
%
%    Usage:    bihist(x1,x2)
%              bihist(x1,x2,n)
%              bihist(x1,x2,c)
%              bihist(ax,...)
%              [h1,h2]=bihist(...)
%
%    Description:
%     BIHIST(X1,X2) counts the elements of X1 & X2 in 10 equal width bins
%     and plots 2 histograms (vertical bar plots) with the X1 histogram
%     going up and the X2 histogram going down.  If X1 or X2 are matrices,
%     the individual columns are counted and plotted side-by-side in the
%     histogram.
%
%     BIHIST(X1,X2,N) if N is scalar, uses N bins.
%
%     BIHIST(X1,X2,C) uses C as the centers for the bins.  Note the first
%     and last bins extend to -/+ Inf.
%
%     BIHIST(AX,...) plots into the axis given by handle AX instead of in
%     the current axis GCA.
%
%     [H1,H2]=BIHIST(...) returns the handles to the histogram objects.
%
%    Notes:
%     - Both directions are labeled positive but there are some quirks:
%       - using YLIM will not update the tick locations - zoom in & out
%         afterwards to update the ticks.
%       - YLIM requires negative values for the lower histogram
%       - data cursor returns negative values for the lower histogram
%
%    Examples:
%     % Compare a normal and uniform distribution
%     bihist(8*(rand(10000,1)-.5),randn(10000,1),100)
%
%    See also: HIST, HIST3, HISTC, BAR

%     Version History:
%        Mar. 17, 2009 - initial version (FEX submission #23312)
%        Sep. 18, 2009 - minor updates
%        June 18, 2011 - rewrite by gge (histogram handle output, more
%                        flexible inputs, ticks redraw on zoom/pan)
%        Jan. 27, 2014 - use axparse instead of axescheck for octave
%
%     Written by Mauro
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 18:55 GMT

% todo:

% get axes handle if there is one
[ax,varargin]=axparse(varargin{:});
nargs=numel(varargin);
if(isempty(ax)); ax=gca; end

% check nargin
error(nargchk(2,3,nargs));

% default n to 10
if(nargs==2 || isempty(varargin{3})); varargin{3}=10; end

% check inputs
if(~(isnumeric(varargin{1}) || islogical(varargin{1})) ...
        || ~(isnumeric(varargin{2}) || islogical(varargin{2})) ...
        || ~(isnumeric(varargin{3}) || islogical(varargin{3})))
    error('seizmo:bihist:badInputs',...
        'All bihist inputs must be numeric!');
end

% combine to get bin centers
[tmp,bins]=hist([varargin{1}(:); varargin{2}(:)],varargin{3});

% now get the bin heights
h1=hist(varargin{1},bins);
h2=hist(varargin{2},bins);

% flip histogram of dataset 2 to be negative
h2=-h2;

% plot the histograms
top=bar(ax,bins,h1,'style','hist');
hold(ax,'on');
bot=bar(ax,bins,h2,'style','hist');
hold(ax,'off');

% change the color of the histograms for those who read this
set(top,'facecolor','y');
set(bot,'facecolor','g');

% change negative y ticks to positive
set(ax,'yticklabelmode','manual','yticklabel',abs(get(ax,'ytick')));

% set id of axes for zoom/pan callbacks
set(ax,'createfcn','bihist');

% setup zoom/pan override
zh=zoom(get(ax,'parent'));
set(zh,'ActionPostCallback',@myzoompostcallback);
ph=pan(get(ax,'parent'));
set(ph,'ActionPreCallback',@myzoompostcallback);
set(ph,'ActionPostCallback',@myzoompostcallback);

% output
if(nargout); varargout={top bot}; end

end

function myzoompostcallback(obj,evd)
% relabel negative y ticks to positive if we are a bihist plot
% - note we remove out of range ticks because otherwise their labels are
%   used in the incorrect location as the ticks don't update during a drag
if(strcmp('bihist',get(evd.Axes,'createfcn')))
    set(evd.Axes,'ytickmode','auto');
    ytick=get(evd.Axes,'ytick');
    ylimit=get(evd.Axes,'ylim');
    set(evd.Axes,'ytick',ytick(ytick>=ylimit(1) & ytick<=ylimit(2)),...
        'yticklabel',abs(ytick(ytick>=ylimit(1) & ytick<=ylimit(2))));
    %set(evd.Axes,'yticklabel',abs(get(evd.Axes,'ytick')));
end
end
