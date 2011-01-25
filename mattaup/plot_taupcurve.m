function [varargout]=plot_taupcurve(tt,offset,flip,varargin)
%PLOT_TAUPCURVE    Plots taupcurve output
%
%    Usage:    plot_taupcurve(tt)
%              plot_taupcurve(tt,offset)
%              plot_taupcurve(tt,offset,flip)
%              plot_taupcurve(tt,offset,flip,'option',value,...)
%              h=plot_taupcurve(...)
%
%    Description:
%     PLOT_TAUPCURVE(TT) plots the travel time curves in TT on a new plot.
%     TT is the output from TAUPCURVE.  See that function for more details
%     on generating travel time curves.
%
%     PLOT_TAUPCURVE(TT,OFFSET) offsets the travel time curves by OFFSET
%     seconds.  This is quite useful when plotting against records that are
%     not relative to the event time.
%
%     PLOT_TAUPCURVE(TT,OFFSET,FLIP) indicates if the x-axis is time or
%     distance.  The default (FALSE) is to have the x-axis as distance.
%
%     PLOT_TAUPCURVE(TT,OFFSET,FLIP,'OPTION',VALUE,...) passes option/value
%     pairs to PLOT.  This is useful for defining the axis ('parent') or
%     the width of the curves ('linewidth').
%
%     H=PLOT_TAUPCURVE(...) returns the handles to the travel time curves.
%
%    Notes:
%
%    Examples:
%     % Align some records on Pdiff, then plot some travel time curves
%     % underneath the records to enhance analysis:
%     evdp=getheader(data(1),'evdp')/1000;
%     data=timeshift(data,-getarrival(data,'Pdiff'));
%     ax=recordsection(data);
%     tt=taupcurve('dep',evdp);
%     idx=find(strcmp('Pdiff',{tt.phase}));
%     intrcpt=tt(idx).time(1)-tt(idx).distance(1)*tt(idx).rayparameter(1);
%     tt=taupcurve('dep',evdp,...
%                  'reddeg',1/tt(idx).rayparameter(1),'ph','ttall');
%     hold(ax,'on');
%     h=plot_taupcurve(tt,-intrcpt,true,'parent',ax,'linewidth',5);
%     movekids(h,'back');
%     hold(ax,'off');
%
%    See also: TAUPCURVE, RECORDSECTION

%     Version History:
%        Jan. 23, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2011 at 17:15 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>3 && ~mod(nargin,2))
    error('matTaup:plot_taupcurve:badNumOptions',...
        'Unpaired option(s)!');
end

% check taupcurve output
if(~isstruct(tt) || any(~isfield(tt,...
        {'time' 'distance' 'depth' 'phase' 'rayparameter'})))
    error('matTaup:plot_taupcurve:badInput',...
        'TT was not created by TAUPCURVE!');
end

% default/check offset
if(nargin<2 || isempty(offset)); offset=0; end
if(~isscalar(offset) || ~isreal(offset))
    error('matTaup:plot_taupcurve:badInput',...
        'OFFSET must be a real-valued scalar!');
end

% default/check flip
if(nargin<3 || isempty(flip)); flip=false; end
if(~isscalar(flip) || (~isreal(flip) && ~islogical(flip)))
    error('matTaup:plot_taupcurve:badInput',...
        'FLIP must be TRUE or FALSE!');
end

% plot curves
ntt=numel(tt);
colors=hsv(ntt)/3;
h=nan(ntt,1);
first=true;
ax=-1;
for i=1:ntt
    if(isempty(tt(i).time)); continue; end
    if(flip) % time vs dist
        h(i)=plot(tt(i).time+offset,tt(i).distance,...
            'color',colors(i,:),varargin{:});
    else     % dist vs time
        h(i)=plot(tt(i).distance,tt(i).time+offset,...
            'color',colors(i,:),varargin{:});
    end
    set(h(i),'tag',tt(i).phase);
    if(first)
        ax=get(h(i),'parent');
        hold(ax,'on');
        first=false;
    end
end
if(ishandle(ax))
    hold(ax,'off');
else
    error('matTaup:plot_taupcurve:badInput',...
        'TT contains no travel time curves!');
end

% output if desired
if(nargout); varargout{1}=h(~isnan(h)); end

end
