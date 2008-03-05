function [conf]=pconf()
%PCONF    Returns default configuration structure for SAClab plotting
%
%    Description: Returns a structure for controlling plotting parameters
%     in SAClab plotting functions.
%
%    Usage: CONF=pconf()
%
%    Examples:
%
%    See also: p1, p2, p3, recsec

% SETTING CONFIGURE STRUCTURE DEFAULTS
conf=struct(...
...% MAIN PLOTTING
'NAME','Seismogram Plotting Utility',...    % FIGURE NAME
'NUMBERTITLE','off',...                     % DISPLAY FIGURE NUMBER
'POINTER','crosshair',...                   % POINTER STYLE
'BGCOLOR','k',...                           % BACKGROUND COLOR
'FGCOLOR','w',...                           % FONT AND AXIS COLOR
'FONTSIZE',6,...                            % FONT SIZE
'FONTWEIGHT','light',...                    % FONT WEIGHT
'FONTNAME','arial',...                      % FONT TYPE
'TRACEWIDTH',1,...                          % TRACE LINE WIDTH
'TICKDIR','out',...                         % TICK DIRECTION
'TICKLEN',[0 0],...                         % TICK LENGTH [2D 3D]
'AXIS','tight',...                          % AXIS SCALING STYLE
'BOX','on',...                              % AXIS BOX
'GRID','on',...                             % GRIDDING
'COLORMAP','hsv',...                        % PLOT COLORING 
...
...% MORE DETAILS
'LINEWIDTH',3,...               % WIDTH OF OTHER LINES
'OCOLOR',[1 0.5 0],...          % COLOR OF O FIELD
'ACOLOR','g',...                % COLOR OF A FIELD
'FCOLOR','r',...                % COLOR OF F FIELD
'TCOLOR','y',...                % COLOR OF T FIELDS
'OTHERCOLOR',[0.5 0.5 0.5],...  % COLOR OF ANYTHING ELSE
'LABELXPAD',0.02,...            % PADDING IN X DIRECTION FOR LABEL
'LABELYPAD',0.05...             % PADDING IN Y DIRECTION FOR LABEL
);

end
