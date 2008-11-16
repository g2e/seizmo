function [conf]=pconf()
%PCONF    Returns default configuration structure for SEIZMO plotting
%
%    Description: Returns a structure for controlling plotting parameters
%     in SEIZMO plotting functions.  Lists all parameters currently 
%     supported and their defaults.
%
%    Usage: CONF=pconf()
%
%    Examples:
%     pconf contains the default plotting parameters that can be overridden
%     in plotting functions by the global SEIZMO and also by function call
%     options.  Once all these are done pconffix should be run to fill in 
%     unset options that need to be replaced by general ones.
%
%     Order of operations:
%       1. CONF=pconf()               (defaults)
%       2. SEIZMO --> CONF            (user global settings)
%       3. function options --> CONF  (user function-specific settings)
%       4. CONF=pconffix(CONF)        (sets unset options using heirarchy)
%       5. plotting
%
%    See also: p1, p2, p3, recsec

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 16:20 GMT

% SETTING CONFIGURE STRUCTURE DEFAULTS
conf=struct(...
...% FIGURE PROPERTIES
'NAME','Seismogram Plotting Utility',...    % FIGURE NAME
'MENUBAR','figure',...                      % DISPLAY MENUBAR
'TOOLBAR','auto',...                        % DISPLAY TOOLBAR
'NUMBERTITLE','off',...                     % DISPLAY FIGURE NUMBER
'RENDERER','painters',...                   % PLOT RENDERER
'RENDERERMODE','auto',...                   % PLOT RENDERER MODE
'DOUBLEBUFFER','on',...                     % DOUBLE THE BUFFER
'POINTER','crosshair',...                   % POINTER STYLE
'FIGBGCOLOR',[],...                         % FIGURE COLOR
...
...% GENERAL
'BGCOLOR','k',...                           % BACKGROUND COLOR
'FGCOLOR','w',...                           % FONT AND AXIS COLOR
'FONTSIZE',6,...                            % FONT SIZE
'FONTWEIGHT','light',...                    % FONT WEIGHT
'FONTNAME','arial',...                      % FONT TYPE
'FONTCOLOR',[],...                          % FONT COLOR
...
...% AXIS PROPERTIES
'AXISFONT',[],...                           % AXIS TICKS FONT
'AXISFONTWEIGHT',[],...                     % AXIS TICKS FONT WEIGHT
'AXISFONTSIZE',[],...                       % AXIS TICKS FONT SIZE
'AXISLINEWIDTH',0.5,...                     % AXIS LINE WIDTH
'TICKDIR','out',...                         % TICK DIRECTION
'TICKLEN',[0 0],...                         % TICK LENGTH [2D 3D]
'XMINORTICK','off',...                      % X-AXIS MINOR TICKS
'YMINORTICK','off',...                      % Y-AXIS MINOR TICKS
'AXIS',{{'tight'}},...                      % AXIS PROPERTIES
'BOX','on',...                              % AXIS BOX
'GRID','on',...                             % GRIDDING
'XLIMITS',[],...                            % X-AXIS LIMITS
'XAXISCOLOR',[],...                         % X-AXIS COLOR
'XAXISLOCATION','bottom',...                % X-AXIS LOCATION
'YLIMITS',[],...                            % Y-AXIS LIMITS
'YAXISCOLOR',[],...                         % Y-AXIS COLOR
'YAXISLOCATION','left',...                  % Y-AXIS LOCATION
'PLOTBGCOLOR',[],...                        % BACKGROUND COLOR OF PLOT
...
...% AXIS LABELS
'TITLE',[],...                              % PLOT TITLE & SETTINGS
'TITLEFONT',[],...
'TITLEFONTCOLOR',[],...
'TITLEFONTSIZE',[],...
'TITLEFONTWEIGHT',[],...
'TITLEINTERP','tex',...
'XLABEL',[],...                             % XAXIS LABEL & SETTINGS
'XLABELFONT',[],...
'XLABELFONTCOLOR',[],...
'XLABELFONTSIZE',[],...
'XLABELFONTWEIGHT',[],...
'XLABELINTERP','tex',...
'YLABEL',[],...                             % YAXIS LABEL & SETTINGS
'YLABELFONT',[],...
'YLABELFONTCOLOR',[],...
'YLABELFONTSIZE',[],...
'YLABELFONTWEIGHT',[],...
'YLABELINTERP','tex',...
...
...% RECORDS
'RECWIDTH',1,...                            % RECORD LINE WIDTH
'COLORMAP','hsv',...                        % RECORD COLORING
...
...% NORMALIZATION
'NORMSTYLE','single',...                    % NORMALIZATION STYLE (SINGLE OR GROUP)
'NORM2YRANGE',true,...                      % NORMALIZE TO Y-AXIS RANGE (TRUE) OR UNITS (FALSE)
'NORMMAX',1/10,...                          % NORMALIZER (Y-AXIS FRACTION OR UNITS)
...
...% LEGEND PROPERTIES
'LEGEND',false,...                          % SHOW LEGEND IF TRUE
'LEGENDINTERP','none',...                   % LEGEND TEXT INTERPOLATION
'LEGENDLOCATION','best',...                 % LEGEND LOCATION
'LEGENDBOX','on',...                        % LEGEND BORDER
'LEGENDBOXCOLOR',[],...                     % LEGEND BORDER COLOR
'LEGENDBOXWIDTH',0.5,...                    % LEGEND BORDER WIDTH
'LEGENDBGCOLOR',[],...                      % LEGEND BACKGROUND COLOR
'LEGENDFONT',[],...                         % LEGEND FONT SETTINGS
'LEGENDFONTSIZE',[],...
'LEGENDFONTWEIGHT',[],...
'LEGENDFONTCOLOR',[],...
...
...% MARKERS
'MARKERWIDTH',2,...                         % WIDTH OF MARKER LINES
'MARKERS',true,...                          % PLOT HEADER MARKERS
'MARKERLABELS',true,...                     % PLOT MARKER LABELS
'OCOLOR',[1 0.5 0],...                      % COLOR OF O MARKER & LABEL
'OMARKERFONTCOLOR',[1 0.5 0],...
'ACOLOR','g',...                            % COLOR OF A MARKER & LABEL
'AMARKERFONTCOLOR','g',...
'FCOLOR','r',...                            % COLOR OF F MARKER & LABEL
'FMARKERFONTCOLOR','r',...
'TCOLOR','y',...                            % COLOR OF T MARKERS & LABELS
'TMARKERFONTCOLOR','y',...
'OTHERCOLOR',[0.5 0.5 0.5],...              % COLOR OF ANYTHING ELSE (?)
'MARKXPAD',0.02,...                         % PADDING IN X DIRECTION FOR MARK (NORMED TO RANGE)
'MARKYPAD',0.05,...                         % PADDING IN Y DIRECTION FOR MARK (NORMED TO RANGE)
'MARKERFONT',[],...                         % MARKER FONT SETTINGS
'MARKERFONTSIZE',[],...
'MARKERFONTWEIGHT',[],...
...
...% P1 SPECIFIC
'NCOLS',[],...                              % NUMBER OF SUBPLOT COLUMNS
...
...% P2 SPECIFIC
'P2NORM',false,...                          % NORMALIZE RECORDS IN P2
...
...% P3 SPECIFIC
'NAMESONYAXIS',false,...                    % DISPLAY RECORD NAMES ON Y-AXIS
...
...% RECSEC SPECIFIC
'YFIELD','gcarc',...                        % HEADER FIELD THAT DEFINES RECORDS Y-AXIS POSITION
...
...% PDENDRO SPECIFIC
'TREELIMIT',0.03,...                        % DISSIMILARITY LIMIT FOR COLORING
'TREELIMITCOLOR','r',...                    % COLOR OF DISSIMILARITY LIMIT MARKER
'DISSIM',true,...                           % X-AXIS AS DISSIMILARITY (FALSE = SIMILARITY)
'TREECOLORMAP','hsv',...                    % GROUP COLORING IN DENDROGRAM
'TREEDEFCOLOR',[0.5 0.5 0.5],...            % DEFAULT COLORING OF TREE
'TREETITLE',[],...                          % TITLE FOR DENDROGRAM PLOT
'TREEXLABEL',[],...                         % XLABEL FOR DENDROGRAM PLOT
...
...% FIGURE HANDLES
'FIGHANDLE',[],...                          % FIGURE HANDLE
'SUBHANDLE',[]...                           % SUBPLOT HANDLE
...
);

% GENERAL RELATIONS
conf.GENERAL=struct(...
'BGCOLOR',{{'FIGBGCOLOR' 'PLOTBGCOLOR' 'LEGENDBGCOLOR'}},...
'FGCOLOR',{{'FONTCOLOR' 'XAXISCOLOR' 'YAXISCOLOR' 'LEGENDBOXCOLOR'...
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

end
