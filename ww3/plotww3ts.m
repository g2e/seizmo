function [varargout]=plotww3ts(varargin)
%PLOTWW3TS    Plots WaveWatch III hindcast data as a time-series
%
%    Usage:    plotww3ts(s)
%              plotww3ts(s,...)
%              plotww3ts(ax,...)
%              h=plotww3ts(...)
%
%    Description:
%     PLOTWW3TS(S) plots the WaveWatch III hindcast data in the structure S
%     created by WW3STRUCT as a time-series rather than as a map or map
%     movie.  All data types are plotted in in the same axes (wind u/v
%     components are both in the plot).
%
%     PLOTWW3TS(S,...) passes additional inputs onto PLOT.
%
%     PLOTWW3TS(AX,...) sets the axes to draw in.  This is useful for
%     subplots, guis, etc.  The default draws in the current axes.
%
%     H=PLOTWW3TS(...) returns the line handles.
%
%    Notes:
%
%    Examples:
%     % Read in some WW3 primary wave period hindcast data at a specific
%     % location and then plot it up:
%     s=ww3struct('nww3.tp.200608.grb',[],[],[],[0 0],[0 1]);
%     h=plotww3ts(s);
%     datetick(get(h,'parent'),'x');
%
%    See also: LOCATE_STORMS, WW3STRUCT, WW3CAT, PLOTWW3

%     Version History:
%        May  31, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  31, 2012 at 15:05 GMT

% todo:

% check number of inputs
error(nargchk(1,inf,nargin));

% find ww3 struct
if(isww3(varargin{1}))
    ax={};
    s=varargin{1};
    varargin(1)=[];
elseif(isww3(varargin{2}))
    ax=varargin(1);
    s=varargin{2};
    varargin(1:2)=[];
else
    error('seizmo:plotww3ts:badInput',...
        'WW3 struct input required!');
end

% require scalar struct
if(~isscalar(s))
    error('seizmo:plotww3ts:badWW3',...
        'PLOTWW3TS can only handle 1 file!');
end

% plot
h=[];
for i=1:numel(s.data)
    h=cat(2,h,plot(ax{:},s.time(:),s.data{i}(:,:),varargin{:}));
end

% output
if(nargout); varargout={h}; end

end

function [tf]=isww3(s)
valid={'path' 'name' 'description' 'units' 'data' ...
    'lat' 'lon' 'time' 'latstep' 'lonstep' 'timestep'};
tf=all(ismember(valid,fieldnames(s)));
end

