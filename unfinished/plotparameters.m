function [option]=plotparameters(setglobal,varargin)
%PLOTPARAMETERS    Handles SEIZMO plot function options

% wishlist:
%%%%%%%%%%%%%%%%
% global options -- general (plotparameters) and function specific (plotX)
% group options -- target specific parts of a figure
%               -- figure, subplot, plot, title, (xyz)label, legend, marker
% 'all' options -- set unset options in that category (allbgcolor sets
%                  the background color for all groups) (note that this
%                  sets the 'color' field even though it is 'bgcolor')
% - markers to be in the background (plot first or move to back)
% - marker text location settable (top,bottom,alternate)
% - text and spacing is normalized to plot
% - way of displaying info when a trace or subplot is clicked

% - what is difference between groups?:
%   axis (not a group -- part of special)
%   plot (axes properties of the plot)
%   record (properties associated with function drawing the records)
%   special (properties handled by the calling function)

% - how to handle things like
%      axis tight/auto/etc, xlim, ylim, normalization, colormap

% how to set the title,xyzlabel?
% - titlestring or xyzlabelstring

% 

% input check
if(~mod(nargin,2))
    error('seizmo:plotparameters:noValue',...
        'Plotting parameters must be field / value pairs!');
end

% who is calling?
fname=dbstack;
fname=fname.name(min(2,size(fname)));

% group field setup
gf={'all' 'figure' 'subplot' 'plot' 'record' ...
    'label' 'title' 'xlabel' 'ylabel' 'zlabel' ...
    'legend' 'markerlabel' 'amarkerlabel' 'fmarkerlabel' ...
    'omarkerlabel' 'tmarkerlabel' 'marker' 'amarker' ...
    'fmarker' 'omarker' 'tmarker' 'special'};
ugf=upper(gf);
len1=cellfun('prodofsize',gf)+1;

% general defaults
option.ALL={'linewidth' 0.5 'units' 'normalized' ...
    'bgcolor' 'k' 'fgcolor' 'w' 'fontsize' 0.1 'interpreter' 'tex'};
option.FIGURE={'name' [fname ' -- SEIZMO Plotting Utility'] ...
    'numbertitle' 'off' 'pointer' 'crosshair'};
option.SUBPLOT={};
option.PLOT={'box' 'on' 'grid' 'on' 'ticklength' [0 0]};
option.RECORD={'linewidth' 1};
option.LABEL={};
option.TITLE={};
option.XLABEL={};
option.YLABEL={};
option.ZLABEL={};
option.LEGEND={'box' 'on' 'interpreter' 'none' 'location' 'best'};
option.MARKER={'linewidth' 2};
option.AMARKER={'color' 'g'};
option.FMARKER={'color' 'r'};
option.OMARKER={'color' [1 0.5 0]};
option.TMARKER={'color' [0.5 0.5 0.5]};
option.MARKERLABEL={};
option.AMARKERLABEL={'fontcolor' 'g'};
option.FMARKERLABEL={'fontcolor' 'r'};
option.OMARKERLABEL={'fontcolor' [1 0.5 0]};
option.TMARKERLABEL={'fontcolor' [0.5 0.5 0.5]};
option.SPECIAL={'normstyle' 'single' 'norm2yaxis' true 'normmax' 0.1 ...
    'showlegend' false 'showmarkers' true 'showmarkerlabels' true};

% specific function defaults
fname=upper(fname);
switch fname
    case 'PLOT0'
        option.SPECIAL={'namesonyaxis' false};
    case 'PLOT1'
        option.SPECIAL={};
    case 'PLOT2'
        option.SPECIAL={'norm' false};
    case 'RECSEC'
        option.SPECIAL={'yfield' 'gcarc'};
    case 'PLOTDENDRO'
        option.SPECIAL={};
    otherwise
        error('seizmo:plotparameters:badFunctionName',...
        'Unknown calling function: %s !',fname);
end

% get options set by SEIZMO global
% general
global SEIZMO;
if(isfield(SEIZMO,'PLOTPARAMETERS'))
    % loop thru group fields
    for i=ugf
        j=i{:};
        if(isfield(SEIZMO.PLOTPARAMETERS,j));
            option.(j)={option.(j){:} SEIZMO.PLOTPARAMETERS.(j){:}};
        end
    end
end
% calling function specific
if(isfield(SEIZMO,fname))
    % loop thru group fields
    for i=ugf
        j=i{:};
        if(isfield(SEIZMO.(fname),j));
            option.(j)={option.(j){:} SEIZMO.(fname).(j){:}};
        end
    end
end

% get options passed with the function call
fields=strvcat(varargin{1:2:end}); %#ok<VCAT>
% loop thru group fields
for i=1:numel(gf)
    % find matches
    f=strmatch(gf{i},fields);
    v=2*f; v1=v-1;
    
    % skip if none
    if(isempty(f))
        continue;
    end
    
    % trim fields
    varargin(v1)=cellstr(fields(f,len1(i):end));
    
    % save
    option.(ugf{i})={option.(ugf{i}){:} varargin([v1; v])};
    
    % trim fields, varargin
    fields(f,:)=[];
    varargin([v1 v])=[];
end

% handle specific fields
if(~isempty(varargin))
    option.SPECIAL={option.SPECIAL{:} varargin{:}};
end

% handle heirarchy
% push general options on front (so specific options overwrite)
option.HEIRARCHY=struct(...
'bgcolor',{{'figurecolor' 'plotcolor' 'legendcolor'}},...
'fgcolor',{{'allfontcolor' 'plotxaxiscolor' 'plotyaxiscolor'...
    'plotzaxiscolor' 'legendcolor' 'markercolor' 'markerlabelcolor' ...
    'OCOLOR' 'ACOLOR' 'FCOLOR' 'TCOLOR' 'OTHERCOLOR'}},...
'FONTCOLOR',{{'TITLEFONTCOLOR' 'XLABELFONTCOLOR'...
    'YLABELFONTCOLOR' 'LEGENDFONTCOLOR' 'OMARKERFONTCOLOR'...
    'AMARKERFONTCOLOR' 'FMARKERFONTCOLOR' 'TMARKERFONTCOLOR'}},...
'FONTNAME',{{'AXISFONT' 'TITLEFONT' 'XLABELFONT'...
    'YLABELFONT' 'LEGENDFONT' 'MARKERFONT'}},...
'FONTSIZE',{{'AXISFONTSIZE' 'TITLEFONTSIZE' 'XLABELFONTSIZE'...
    'YLABELFONTSIZE' 'LEGENDFONTSIZE' 'MARKERFONTSIZE'}},...
'FONTWEIGHT',{{'AXISFONTWEIGHT' 'TITLEFONTWEIGHT' 'XLABELFONTWEIGHT'...
    'YLABELFONTWEIGHT' 'LEGENDFONTWEIGHT' 'MARKERFONTWEIGHT'}});


option.ALL={'linewidth' 0.5 'units' 'normalized' ...
    'bgcolor' 'k' 'fgcolor' 'w' 'fontsize' 0.1 'interpreter' 'tex'};
% linewidth units interpreter bgcolor fgcolor
% fontname fontweight fontcolor fontunits fontsize



end
