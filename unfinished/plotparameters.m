function [option]=plotparameters(setglobal,varargin)
%PLOTPARAMETERS    Handles SEIZMO plot function options

% I want to use functions called getplotopt/setplotopt to alter the
% plotting parameters.  Maybe it would be better to just have a plotopt
% function to set the parameters and returns the overall settings if
% desired.  This is what I do for my check functions.
%
% Q. How do I specify parameters that affect all plot types vs just a
%    specific type?
% A. First off, parameters passed directly to plotting functions only
%    affect that one plot.  To affect all plots the settings must be
%    changed through plotopt or the underlying function plotparameters.
%    To choose to affect a certain type of plot use 'plottype*option'.  If
%    no plottype is given then we assume it is not specific.
% Q. How do I target my options at specific matlab functions in the plots?
% A. Use 'matlabfunctionoption' where there is no separating character
%    between the two.
% Q. What if I want to target just a specific element?
% A. Many elements follow the same form: 'elementoption'.
%
% 'plottype*elem/func*option',value

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

% so what is the order
% - alldefaults defaults general specific

% input check
if(~mod(nargin,2))
    error('seizmo:plotparameters:noValue',...
        'Plotting parameters must be field / value pairs!');
end

% who is calling?
fname=dbstack;
fname=fname(min(2,size(fname))).name;
if(strcmpi(fname,'plotopt')); fname='plotparameters'; end

% group field setup
gf={'figure' 'axes' 'subplot' 'line' 'text' 'record' ...
    'label' 'title' 'xlabel' 'ylabel' 'zlabel' ...
    'legend' 'markerlabel' 'amarkerlabel' 'fmarkerlabel' ...
    'omarkerlabel' 'tmarkerlabel' 'marker' 'amarker' ...
    'fmarker' 'omarker' 'tmarker'};
ugf=upper(gf);
len1=cellfun('prodofsize',gf)+1;

% general defaults
option.ALL=struct(...
    'LINEWIDTH',0.5,...
    'BGCOLOR','k',...
    'FGCOLOR','w',...
    'FONTNAME',[],...
    'FONTWEIGHT',[],...
    'FONTANGLE',[],...
    'FONTUNITS',[],...
    'FONTSIZE',[],...
    'FONTCOLOR',[],...
    'INTERPRETER',[]);
option.FIGURE={'name' ['SEIZMO -- ' upper(fname)] ...
    'numbertitle' 'off' 'pointer' 'crosshair'};
option.AXES={'box' 'on' 'xgrid' 'on' 'ygrid' 'on' 'ticklength' [0 0]};
option.AXIS={};

%option.CONTOUR={};
%option.STEM={};
%option.BAR={};
%option.STAIR={};
%option.AREA={};
%option.IMAGE={};
%option.SURFACE={};
%option.QUIVER={};
%option.RECTANGLE={};
%option.LIGHT={};
%option.ERRORBAR={};
%option.PATCH={};

option.SUBPLOT={};
option.LINE={};
option.TEXT={};
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

% SPECIAL must be a struct
option.SPECIAL=struct(...
    'HANDLE',[],...
    'COLORMAP','hsv',...
    'NORMSTYLE','single',...
    'NORM2YAXIS',true,...
    'NORMMAX',0.1,...
    'SHOWLEGEND',false,...
    'SHOWMARKERS',true,...
    'SHOWMARKERLABELS',true,...
    'SHOWINFO',false);

% specific function defaults
fname=upper(fname);
switch fname
    case 'PLOTEVEN'
        option.SPECIAL.NAMESONYAXIS=false;
    case 'PLOT1'
        option.NCOLS=[];
    case 'PLOTOVERLAY'
        option.SPECIAL.NORM2YAXIS=false;
    case 'RECORDSECTION'
        option.SPECIAL.YFIELD='gcarc';
    case 'PLOTDENDRO'
        option.SPECIAL.TREELIMIT=0.2;
        option.SPECIAL.TREELIMITCOLOR='r';
        option.SPECIAL.TREECOLORMAP=option.SPECIAL.COLORMAP;
        option.SPECIAL.TREEDEFCOLOR=[0.5 0.5 0.5];
        option.SPECIAL.TREETITLE=[];
        option.SPECIAL.TREEXLABEL=[];
end

% get options set in SEIZMO global
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
%%%%%%%%
%% STILL NEED TO PARSE SPECIFIC
%% STILL NEED TO SET SEIZMO GLOBAL HERE
%%%%%%%%
fields=strvcat(varargin{1:2:end}); %#ok<VCAT>
% loop thru group fields
for i=1:numel(gf)
    % find matches
    f=strmatch(gf{i},fields);
    v=2*f; v1=v-1;
    v0=[v1 v]'; v0=v0(:);
    
    % skip if none
    if(isempty(f)); continue; end
    
    % trim fields
    varargin(v1)=cellstr(fields(f,len1(i):end));
    
    % save
    option.(ugf{i})={option.(ugf{i}){:} varargin(v0)};
    
    % trim fields, varargin
    fields(f,:)=[];
    varargin(v0)=[];
end

% get options passed to axis
%%%%%%%%
%% STILL NEED TO PARSE SPECIFIC
%% STILL NEED TO SET SEIZMO GLOBAL HERE
%%%%%%%%

% set the rest as special
%%%%%%%%
%% STILL NEED TO PARSE SPECIFIC
%% STILL NEED TO SET SEIZMO GLOBAL HERE
%%%%%%%%

% deal with all heirarchy


%{
% handle specific fields
if(~isempty(varargin))
    option.SPECIAL={option.SPECIAL{:} varargin{:}};
end

% handle heirarchy of all
% push general options on front (so specific options overwrite)
option.AXES={'fontname' option.ALL.FONTNAME ...
    'fontunits' option.ALL.FONTUNITS 'fontsize' option.ALL.FONTSIZE ...
    'fontweight' option.ALL.FONTWEIGHT 'fontangle' option.ALL.FONTANGLE ...
    'linewidth' option.ALL.LINEWIDTH 'color' option.ALL.BGCOLOR ...
    'xcolor' option.ALL.FGCOLOR 'ycolor' option.ALL.FGCOLOR ...
    'zcolor' option.ALL.FGCOLOR option.AXES{:}};
option.TEXT={'interpreter'  option.TEXT{:}};

option.HEIRARCHY=struct(...
'LINEWIDTH',{{'AXESLINEWIDTH','LINELINEWIDTH'}},...
'BGCOLOR',{{'FIGURECOLOR' 'AXESCOLOR' 'LEGENDCOLOR'}},...
'FGCOLOR',{{'ALLFONTCOLOR' 'AXESXCOLOR' 'AXESYCOLOR'...
    'AXESZCOLOR' 'LEGENDEDGECOLOR' 'MARKERCOLOR'}},...
'INTERPRETER',{{'LEGENDINTERPRETER' 'TEXTINTERPRETER'}},...
'FONTCOLOR',{{'TEXTCOLOR'}},...
'FONTUNITS',{{'AXESFONTUNITS' 'TEXTFONTUNITS'}},...
'FONTNAME',{{'AXESFONTNAME' 'TEXTFONTNAME'}},...
'FONTSIZE',{{'AXESFONTSIZE' 'TEXTFONTSIZE'}},...
'FONTWEIGHT',{{'AXESFONTWEIGHT' 'TEXTFONTWEIGHT'}},...
'FONTANGLE',{{'AXESFONTANGLE' 'TEXTFONTWEIGHT'}});
%}


end
