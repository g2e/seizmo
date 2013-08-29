function [varargout]=plot_taupcurve_dt(tt,offset,flip,varargin)
%PLOT_TAUPCURVE_DT    Distance-time plot of taupcurve output
%
%    Usage:    plot_taupcurve_dt(tt)
%              plot_taupcurve_dt(tt,offset)
%              plot_taupcurve_dt(tt,offset,flip)
%              plot_taupcurve_dt(tt,offset,flip,'option',value,...)
%              h=plot_taupcurve_dt(...)
%
%    Description:
%     PLOT_TAUPCURVE_DT(TT) plots the travel time curves in TT in the
%     current axes as distance vs traveltime.  TT is the output from
%     TAUPCURVE.  See that function for more details on generating travel
%     time curves.
%
%     PLOT_TAUPCURVE_DT(TT,OFFSET) offsets the travel time curves by OFFSET
%     seconds.  This is quite useful when overlaying curves against records
%     that are not relative to the origin time.
%
%     PLOT_TAUPCURVE_DT(TT,OFFSET,FLIP) indicates if the x-axis is time or
%     distance.  The default (FALSE) is to have the x-axis as distance.
%
%     PLOT_TAUPCURVE_DT(TT,OFFSET,FLIP,'OPTION',VALUE,...) passes
%     option/value pairs to PLOT.  This is useful for defining the axis
%     ('parent') or the width of the curves ('linewidth').
%
%     H=PLOT_TAUPCURVE_DT(...) returns the handles to the travel time
%     curves.
%
%    Notes:
%
%    Examples:
%     % Align some records on Pdiff, then plot some travel time curves
%     % underneath the records to enhance analysis:
%     evdp=getheader(data(1),'evdp'); % depth in meters!
%     data=timeshift(data,-findpicks(data,'Pdiff',1));
%     ax=recordsection(data);
%     tt=taupcurve('dep',evdp/1000); % depth in kilometers!
%     idx=find(strcmp('Pdiff',{tt.phase}));
%     intrcpt=tt(idx).time(1)-tt(idx).distance(1)*tt(idx).rayparameter(1);
%     tt=taupcurve('dep',evdp/1000,...
%                  'reddeg',1/tt(idx).rayparameter(1),'ph','ttall');
%     hold(ax,'on');
%     h=plot_taupcurve_dt(tt,-intrcpt,true,'parent',ax,'linewidth',5);
%     movekids(h,'back');
%     hold(ax,'off');
%
%    See also: TAUPCURVE, PLOT_TAUPPATH, RECORDSECTION

%     Version History:
%        Jan. 23, 2011 - initial version
%        May  21, 2011 - minor doc update
%        Feb. 24, 2012 - doc update, code update
%        Mar. 15, 2012 - minor fix in example
%        Aug. 29, 2013 - fixed m=>km bug in example (thanks Cheng!)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 29, 2013 at 17:15 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>3 && ~mod(nargin,2))
    error('TauP:plot_taupcurve_dt:badNumOptions',...
        'Unpaired option(s)!');
end

% check taupcurve output
if(~isstruct(tt) || any(~isfield(tt,...
        {'modelname' 'depth' 'distance' 'mindistance' ...
        'phase' 'puristphase' 'time' 'rayparameter'})))
    error('TauP:plot_taupcurve:badInput',...
        'TT was not created by TAUPCURVE!');
elseif(isempty(tt))
    error('TauP:plot_taupcurve:badInput',...
        'TT is empty!');
end

% default/check offset
if(nargin<2 || isempty(offset)); offset=0; end
if(~isscalar(offset) || ~isreal(offset))
    error('TauP:plot_taupcurve_dt:badInput',...
        'OFFSET must be a real-valued scalar!');
end

% default/check flip
if(nargin<3 || isempty(flip)); flip=false; end
if(~isscalar(flip) || (~isreal(flip) && ~islogical(flip)))
    error('TauP:plot_taupcurve_dt:badInput',...
        'FLIP must be TRUE or FALSE!');
end

% plot curves
ntt=numel(tt);
colors=hsv(ntt)/3;
h=nan(ntt,1);
first=true;
ax=-1;
for i=1:ntt
    if(flip) % time vs dist
        h(i)=plot(tt(i).time+offset,tt(i).distance,...
            'color',colors(i,:),'displayname',tt(i).phase,...
            'tag','taupcurve_timecurve',varargin{:});
    else     % dist vs time
        h(i)=plot(tt(i).distance,tt(i).time+offset,...
            'color',colors(i,:),'displayname',tt(i).phase,...
            'tag','taupcurve_timecurve',varargin{:});
    end
    if(first)
        ax=get(h(i),'parent');
        held=ishold(ax);
        if(~held); hold(ax,'on'); end
        first=false;
    end
end
if(~held); hold(ax,'off'); end

% output if desired
if(nargout); varargout{1}=h; end

end
