function sliceomatic(p1,p2,xmesh,ymesh,zmesh)
% SLICEOMATIC - Slice and isosurface volume exploration GUI
%
% SLICEOMATIC(DATA) - Use 3D double matrix DATA as a volume data
%
% Example:
%
%       [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
%       v = x .* exp(-x.^2 - y.^2 - z.^2);
%       sliceomatic(v)
%
% Using SLICEOMATIC with no arguments is equivalent to the above
% example.
%
% SLICEOMATIC(DATA, X, Y, Z) - Run sliceomatic using the specified data 
% coordinates for the volume DATA. X, Y, and Z are the vectors over which
% DATA is defined.
%
% ex:
%       x = -2:.2:2; y = -2:.25:2; z = -2:.16:2;
%       [X,Y,Z] = meshgrid(x,y,z);
%       v = X .* exp(-X.^2 - Y.^2 - Z.^2);
%       sliceomatic(v,x,y,z)
%
%
% Using the GUI:
% -------------
%
% Create/Delete slices:
% 
% The white bars on the top, left, and right allow insertion of
% new slices on the X, Y, and Z planes.  Click in an empty area to
% add a new slice or surface.  Right click on a control arrow to
% reconfigure or delete that slice.
%
% Create/Delete isosurfaces:
%
% The colored bar at the bottom is used to place and position an
% isosurface.  The color in the bar indicates a position (as seen
% in the slice) where the isosurface will go.  Right click on a
% control arrow to reconfigure or delete the isosurface.
%
% Orientation of the view:
%
% When the rotate camera button is on, the popup menu will control
% the camera.  Turn off camera rotation in order to get individual
% control over properties of the slices and isosurfaces.
%
% Changing Defaults:
%
% The defaults menu provides default features of newly created
% slices and surfaces.  The AllSlices menu controls properties of
% all the slices and surfaces in the scene.  Use popup menus on the
% objects themselves, or the control arrows to change individual
% properties.
%
% Color & Alpha Maps:
%
% The Colormap popdown controls the currently active colormap.
% This map is used to color the slices.  The Alphamap popdown
% controls the alphamap used on the slices.
%
% Use the color or alpha maps to change how areas of your data are
% highlighted.
%
% Controls Control:
%
% The Controls menu allows you to adjust how the controls look.  An
% important feature here is the "Animate" item.  You can enable or
% disable an animation when some changes are made.  Since this is
% decorative, it may be important to turn this off for large data
% sets.
%
% Doing Cool Stuff:
% ----------------
%
% Exploration:
% You can get a quick feel of the current data set by adding a
% slice using the ColorTexture option.  Such a slice can be dragged
% through the data very quickly.
%
% Highlight an Area:
% If certain values in your data are interesting (very large, very
% small, or very median values) you can use transparency to make
% parts of your slices disappear.  Choose AlphaTexture options from
% the defaults, and sweep some slices across your data.  Use the
% AlphaMap to pick out certain data sets.  The example given here
% looks best with the `vdown' alphamap.
%
% Contours on slices:
% You can add a contour onto a slice to further extract shapes
% from the data you are exploring.  Auto-selecting contour limits
% will choose contours on a per slice basis.  Auto-selecting
% contour lines from a volume arbitrarily specifies 10 levels
% based on the limits of the volume.
%
% Hidden Shapes:
% Use the isosurface control bar to create an isosurface.  Be
% patient with this control.  It often takes a while for the
% surface to be created.  Click and hold the mouse button down
% until the first surface appears.  Drag the surface through the
% values until you get something you like, then let go.  If your
% data set is very large, you will need to wait while the new and
% more accurate isosurface is created.
%
% Volumes:
% You can simulate a volume object by creating lots of stacked
% slices.  Simply use the proper Alphamap and transparent textures
% to highlight the correct data, and a large stack of slices will
% let you simulate a volume object.
%
% Customized Graphics:
% -------------------
%
% To add your own graphics into the sliceomatic display, whatever
% that may be, you can use the following technique:
%
% 1) click on a control arrow
% 2) use gco to get the data for that object
%    slice = getappdata(gco,'arrowslice')
% 3) use GET to get the cdata and position data which you can use
%    to add your own graphics.
%
% Setting Default Values:
% ----------------------
%
% If you want to change some default setup feature of sliceomatic,
% use the "Save Preferences" menu item.  Future sliceomatic
% sessions will then retrieve those settings.
%
%
% BUGS:
% ----
%
% 1) Inaccurate Slice
%    Sliceomatic does not use the `slice' command.  All slices are
%    created by explicitly extracting data from the volume.  As such,
%    only slices at integer values are allowed.
% 
% 2) Crashes MATLAB
%    Sliceomatic uses the default OpenGL setup.  If you encounter
%    frequent crashes you can start by enabling software OpenGL
%    rendering.  This should always fix the problem, and would
%    likely slow things down too.  You should also update your
%    graphics driver for your video card.  On Windows, in
%    particular, drivers are updated frequently.  For detail on how
%    to overcome these problems, visit this web page:
%  http://www.mathworks.com/support/tech-notes/1200/1201.html
%
% See Also: SLICE, ISOSURFACE, ISOCAPS, CONTOURC, COLORMAP, SURFACE
  
% This is version 2.3 of sliceomatic.
%
% Sliceomatic is a tool I wrote for fun.  There are no warrenties
% expressed or implied.

% Written by Eric Ludlam <eludlam@mathworks.com>
% Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008
%           The MathWorks Inc
%
% Modified by Emiliano Spezi <emiliano.spezi@physics.org> on 22 May 2003
% Added capability: axes limits control, slice leveling controls,
% and contour level specification.
%
% Patch from David Schwartz on Sep 4, 2008
% Warning fix, and isosurface color fixes.
%
% Patch from Gerhard Stoeckel on Nov 17, 2008
% Fix colormap setting to rand.
    
  if nargin==0
    [x,y,z] = meshgrid(-2:.2:2, -2:.25:2, -2:.16:2);
    v = x .* exp(-x.^2 - y.^2 - z.^2);
    sliceomatic(v)
    return
  end

  if isnumeric(p1)
% $$$     if nargin==4
% $$$       d.Xv=p1;
% $$$       d.Yv=p2;
% $$$       d.Zv=p3
% $$$       p1=p4;
% $$$     else
% $$$       d.Yv=1:size(p1,1);
% $$$       d.Xv=1:size(p1,2);
% $$$       d.Zv=1:size(p1,3);
% $$$     end
    
    d.data=p1;
    
    if nargin>=4
      
      if nargin==4
        zmesh=ymesh;
        ymesh=xmesh;
        xmesh=p2;
      end
      
      d = sliceomaticfigure(d,xmesh,ymesh,zmesh);
      d = sliceomaticsetdata(d,xmesh,ymesh,zmesh);
    else
      d = sliceomaticfigure(d);
      d = sliceomaticsetdata(d);
    end
    
    setappdata(gcf,'sliceomatic',d);

  elseif isa(p1,'char')
    % Interpret commands
    d=getappdata(gcf,'sliceomatic');
    try
      switch p1
       case 'Xnew'
        if strcmp(get(gcf,'selectiontype'),'normal')
          pt=get(gcbo,'currentpoint');
          axis(gcbo);
          X=pt(1,1);
          newa=arrow(gcbo,'down',[X 0]);
          set(gcf,'currentaxes',d.axmain);
          new=localslice(d.data, X, [], []);
          setappdata(new,'controlarrow',newa);
          setappdata(newa(2),'arrowslice',new);
          set(new,'alphadata',get(new,'cdata'),'alphadatamapping','scaled');
          set(newa,'buttondownfcn','sliceomatic Xmove');
          set([new newa],'uicontextmenu',d.uic);
          % Make sure whatever buttonupfcn on the figure is run now to "turn
          % off" whatever was going on before we got our callback on the
          % arrow.
          buf = get(gcf,'windowbuttonupfcn');
          if ~strcmp(buf,'')
            eval(buf);
          end
          d.draggedarrow=newa(2);
          dragprep(newa(2));
          setpointer(gcf,'SOM leftright');
          set(d.motionmetaslice,'visible','off');
        end
       case 'Ynew'
        if strcmp(get(gcf,'selectiontype'),'normal')
          pt=get(gcbo,'currentpoint');
          Y=pt(1,2);
          newa=arrow(gcbo,'right',[0 Y]);
          set(gcf,'currentaxes',d.axmain);
          new=localslice(d.data, [], Y, []);
          setappdata(new,'controlarrow',newa);
          setappdata(newa(2),'arrowslice',new);
          set(new,'alphadata',get(new,'cdata'),'alphadatamapping','scaled');
          set(newa,'buttondownfcn','sliceomatic Ymove');
          set([new newa],'uicontextmenu',d.uic);
          % Make sure whatever buttonupfcn on the figure is run now to "turn
          % off" whatever was going on before we got our callback on the
          % arrow.
          buf = get(gcf,'windowbuttonupfcn');
          if ~strcmp(buf,'')
            eval(buf);
          end
          d.draggedarrow=newa(2);
          dragprep(newa(2));
          setpointer(gcf,'SOM topbottom');
          set(d.motionmetaslice,'visible','off');
        end % if strcmp(get(gcf,
       case 'Znew'
        if strcmp(get(gcf,'selectiontype'),'normal')
          pt=get(gcbo,'currentpoint');
          Y=pt(1,2);
          newa=arrow(gcbo,'left', [0 Y]);
          set(gcf,'currentaxes',d.axmain);
          new=localslice(d.data, [], [], Y);
          set(new,'alphadata',get(new,'cdata'),'alphadatamapping','scaled');
          setappdata(new,'controlarrow',newa);
          setappdata(newa(2),'arrowslice',new);
          set(newa,'buttondownfcn','sliceomatic Zmove');
          set([new newa],'uicontextmenu',d.uic);
          % Make sure whatever buttonupfcn on the figure is run now to "turn
          % off" whatever was going on before we got our callback on the
          % arrow.
          buf = get(gcf,'windowbuttonupfcn');
          if ~strcmp(buf,'')
            eval(buf);
          end
          d.draggedarrow=newa(2);
          dragprep(newa(2));
          setpointer(gcf,'SOM topbottom');
          set(d.motionmetaslice,'visible','off');
        end % if strcmp(get(gcf,
       case 'ISO'
        if strcmp(get(gcf,'selectiontype'),'normal')
          pt=get(gcbo,'currentpoint');
          V=pt(1,1);
          newa=arrow(gcbo,'up',[V 0]);
          set(gcf,'currentaxes',d.axmain);
          new=localisosurface(d.reducelims,d.reduce,d.reducesmooth,V);
          set([newa new],'uicontextmenu',d.uiciso);
          setappdata(new,'controlarrow',newa);
          setappdata(new,'reduced',1);
          setappdata(newa(2),'arrowiso',new);
          set(newa,'buttondownfcn','sliceomatic ISOmove');
          % Make sure whatever buttonupfcn on the figure is run now to "turn
          % off" whatever was going on before we got our callback on the
          % arrow.
          buf = get(gcf,'windowbuttonupfcn');
          if ~strcmp(buf,'')
            eval(buf);
          end	
          d.draggedarrow=newa(2);
          dragprep(newa(2));
          setpointer(gcf,'SOM leftright');
        end % if strcmp(get(gcf,
       case 'Xmove'
        if strcmp(get(gcf,'selectiontype'),'normal')
          [a s]=getarrowslice;
          d.draggedarrow=a;
          dragprep(a);
        end
       case 'Ymove'
        if strcmp(get(gcf,'selectiontype'),'normal')
          [a s]=getarrowslice;
          d.draggedarrow=a;
          dragprep(a);
        end
       case 'Zmove'
        if strcmp(get(gcf,'selectiontype'),'normal')
          [a s]=getarrowslice;
          d.draggedarrow=a;
          dragprep(a);
        end
       case 'ISOmove'
        if strcmp(get(gcf,'selectiontype'),'normal')
          [a s]=getarrowslice;
          d.draggedarrow=a;
          dragprep(a);
        end      
       case 'up'
        if strcmp(get(gcf,'selectiontype'),'normal')
          dragfinis(d.draggedarrow);
        end
       case 'motion'
        % Make sure our cursor is ok
        a=d.draggedarrow;			% The arrow being dragged
        s=getappdata(a,'arrowslice');	% The slice to 'move'
        if isempty(s)
          s=getappdata(a,'arrowiso');	% or the isosurface
        end
        aa=get(a,'parent');		% arrow's parent axes
        pos=getappdata(a,'arrowcenter');	% the line the arrow points at.
        apos=get(aa,'currentpoint');
        
        % Bind the axes position to the limits of that axes.
        xlimits = get(aa,'xlim');
        ylimits = get(aa,'ylim');

        if apos(1,1) < xlimits(1)
          apos(1,1) = xlimits(1);
        elseif apos(1,1) > xlimits(2)
          apos(1,1) = xlimits(2);
        end
          
        if apos(1,2) < ylimits(1)
          apos(1,2) = ylimits(1);
        elseif apos(1,2) > ylimits(2)
          apos(1,2) = ylimits(2);
        end
          
        if aa==d.axx || aa==d.axiso
          % We are moving an X slice
          xdiff=apos(1,1)-pos(1,1);
          v=get(a,'vertices');
          v(:,1)=v(:,1)+xdiff;
          set([a getappdata(a,'arrowedge')],'vertices',v);
          np=[ apos(1,1) 0 ];
          % This might be a slice, or an isosurface!
          if aa==d.axiso
            new=localisosurface(d.reducelims,d.reduce,d.reducesmooth,...
                                apos(1,1),s);
            setappdata(new,'reduced',1);
            movetipforarrow(a, aa, apos(1,1), [ apos(1,1) 6 ], 'bottom','center')
          else
            %disp([ 'apos = ' num2str(apos(1,1))])
            %disp([ 'pos  = ' num2str(pos(1,1))])
            %disp([ 'change=' num2str(round(apos(1,1))~=round(pos(1,1)))]);
            if round(apos(1,1))~=round(pos(1,1))
              localslice(d.data, apos(1,1), [], [],s);
            end
            movetipforarrow(a, aa, apos(1,1), [ apos(1,1) .5 ],'top','center')
          end
        else
          % We are moving a Y or Z slice
          ydiff=apos(1,2)-pos(1,2);
          v=get(a,'vertices');
          v(:,2)=v(:,2)+ydiff;
          set([a getappdata(a,'arrowedge')],'vertices',v);
          np=[ 0 apos(1,2) ];
          if aa==d.axy
            if round(apos(1,2))~=round(pos(1,2))
              localslice(d.data, [], apos(1,2), [], s);
            end
            movetipforarrow(a, aa, apos(1,2), [ 5.5 apos(1,2) ], 'middle','left');
          else
            if round(apos(1,2))~=round(pos(1,2))
              localslice(d.data, [], [], apos(1,2), s);
            end
            movetipforarrow(a, aa, apos(1,2), [ .5 apos(1,2) ], 'middle','right');
          end
        end
        setappdata(a,'arrowcenter',np);
        % Question: Is anyone dependant on versions of MATLAB
        % that does not support Java Figures?
        %
        %% This improves anaimation speed in Java Figures/R14
        %% The Rule: Java Figures don't want this drawnow.
        %try
        %  if isempty(get(gcf,'javaframe'))
        %    drawnow;
        %  end
        %catch
        %  drawnow;
        %end
        %
        % IsoSurface context menu items
        %
       case 'isotogglevisible'
        [a s]=getarrowslice;
        if propcheck(s,'visible','on')
          set(s,'visible','off');
        else
          set(s,'visible','on');
        end
       case 'isodelete'
        [a s]=getarrowslice;
        if numel(a)==1
          delete(getappdata(a,'arrowedge'));
        end
        cap=getappdata(s,'sliceomaticisocap');
        if ~isempty(cap)
          delete(cap);
        end
        delete(s);
        delete(a);
       case 'isoflatlight'
        [a s]=getarrowslice;
        set(s,'facelighting','flat');
       case 'isosmoothlight'
        [a s]=getarrowslice;
        set(s,'facelighting','phong');
       case 'isocolor'
        [a s]=getarrowslice;
        c=uisetcolor(get(s,'facecolor'));
        slowset(s,'facecolor',c,d.animincrement);
       case 'isoalpha'
        [a s]=getarrowslice;
        if nargin ~= 2
          error('Not enough arguments to sliceomatic.');
        end
        slowset(s,'facealpha',eval(p2),d.animincrement);
       case 'isocaps'
        [a s]=getarrowslice;
        cap=getappdata(s,'isosurfacecap');
        if isempty(cap)
          new=localisocaps(s);
          set(new,'uicontextmenu',d.uiciso);
        else
          delete(cap);
          setappdata(s,'isosurfacecap',[]);
        end
        %
        % Now for slice context menu items
        %
       case 'togglevisible'
        [a s]=getarrowslice;
        switch get(s,'visible')
         case 'on'
          set(s,'visible','off');
          pushset(a,'facealpha',.2);
         case 'off'
          set(s,'visible','on');
          popset(a,'facealpha');
        end
       case 'setfaceted'
        [a s]=getarrowslice;
        set(s,'edgec','k','facec','flat');
        if ischar(get(s,'facea')) && strcmp(get(s,'facea'),'texturemap')
          set(s,'facea','flat');
        end
        textureizeslice(s,'off');
       case 'setflat'
        [a s]=getarrowslice;
        set(s,'edgec','n','facec','flat');
        if ischar(get(s,'facea')) && strcmp(get(s,'facea'),'texturemap')
          set(s,'facea','flat');
        end
        textureizeslice(s,'off');
       case 'setinterp'
        [a s]=getarrowslice;
        set(s,'edgec','n','facec','interp');
        if ischar(get(s,'facea')) && strcmp(get(s,'facea'),'texturemap')
          set(s,'facea','interp');
        end
        textureizeslice(s,'off');
       case 'settexture'
        [a s]=getarrowslice;
        set(s,'facecolor','texture','edgec','none');
        if ischar(get(s,'facea'))
          set(s,'facealpha','texturemap');
        end
        textureizeslice(s,'on');
       case 'setnone'
        [a s]=getarrowslice;
        set(s,'facecolor','none','edgec','none');
        textureizeslice(s,'off');
       case 'setalphanone'
        [a s]=getarrowslice;
        slowset(s,'facealpha',1,d.animincrement);
       case 'setalphapoint5'
        [a s]=getarrowslice;
        slowset(s,'facealpha',.5,d.animincrement);
       case 'setalphaflat'
        [a s]=getarrowslice;
        set(s,'facealpha','flat');
        if ischar(get(s,'facec')) && strcmp(get(s,'facec'),'texturemap')
          set(s,'facecolor','flat');
          textureizeslice(s,'off');
        end
       case 'setalphainterp'
        [a s]=getarrowslice;
        set(s,'facealpha','interp');
        if ischar(get(s,'facec')) && strcmp(get(s,'facec'),'texturemap')
          set(s,'facecolor','interp');
          textureizeslice(s,'off');
        end
       case 'setalphatexture'
        [a s]=getarrowslice;
        set(s,'facealpha','texturemap');
        if ischar(get(s,'facec'))
          set(s,'facecolor','texturemap');
          textureizeslice(s,'on');
        end
       case 'slicecontour'
        [a s]=getarrowslice;
        localcontour(s, getappdata(s,'contour'));
       case 'slicecontourfullauto'
        [a s]=getarrowslice;
        d = getappdata(gcf, 'sliceomatic');
        minmax = get(d.axiso,'clim');
        levels = minmax(1):(minmax(2)-minmax(1))/10:minmax(2);
        setappdata(s, 'contourlevels', levels);
        localcontour(s, getappdata(s,'contour'),levels);        
       case 'slicecontour_setauto'
        [a s]=getarrowslice;
        setappdata(s, 'contourlevels', []);
        localcontour(s, getappdata(s,'contour'));
       case 'slicecontour_setfullauto'
        [a s]=getarrowslice;
        minmax = get(d.axiso,'clim');
        levels = minmax(1):(minmax(2)-minmax(1))/10:minmax(2);
        setappdata(s, 'contourlevels', levels);
        localcontour(s, getappdata(s,'contour'),levels);
       case 'slicecontour_select'
        [a s]=getarrowslice;
        d = getappdata(gcf, 'sliceomatic');
        xl = get(d.axiso,'xlim');
        levels = selectcontourlevels(get(s,'cdata'), xl(1), xl(2));
        setappdata(s, 'contourlevels', levels);
        localcontour(s, getappdata(s,'contour'),levels);
       case 'slicecontour_setlevels'
        [a s]=getarrowslice;
        d = getappdata(gcf, 'sliceomatic');
        xl = get(d.axiso,'xlim');
        levels = selectcontourlevels(get(s,'cdata'), xl(1), xl(2));
        setappdata(s, 'contourlevels', levels);
        localcontour(s, getappdata(s,'contour'),levels);
       case 'deleteslice'
        [a s]=getarrowslice;
        if prod(size(a))==1
          delete(getappdata(a,'arrowedge'));
        end
        if ~isempty(getappdata(s,'contour'))
          delete(getappdata(s,'contour'));
        end
        delete(s);
        delete(a);
       case 'deleteslicecontour'
        [a s]=getarrowslice;
        if ~isempty(getappdata(s,'contour'))
          delete(getappdata(s,'contour'));
        end
        temp=getappdata(s);
        try 
          temp.contourlevels;
          setappdata(s,'contourlevels',[]);
        end
        setappdata(s,'contour',[]);
       case 'slicecontourflat'
        [a s]=getarrowslice;
        c = getappdata(s,'contour');
        if ~isempty(c)
          set(c,'edgecolor','flat');
        end
       case 'slicecontourinterp'
        [a s]=getarrowslice;
        c = getappdata(s,'contour');
        if ~isempty(c)
          set(c,'edgecolor','interp');
        end
       case 'slicecontourblack'
        [a s]=getarrowslice;
        c = getappdata(s,'contour');
        if ~isempty(c)
          set(c,'edgecolor','black');
        end
       case 'slicecontourwhite'
        [a s]=getarrowslice;
        c = getappdata(s,'contour');
        if ~isempty(c)
          set(c,'edgecolor','white');
        end
       case 'slicecontoursmooth'
        [a s]=getarrowslice;
        c = getappdata(s,'contour');
        onoff = get(gcbo,'checked');
        switch onoff
         case 'off'
          set(c,'linesmoothing','on');
         case 'on'
          set(c,'linesmoothing','off');
        end
       case 'slicecontourcolor'
        [a s]=getarrowslice;
        c = getappdata(s,'contour');
        if ~isempty(c)
          inputcolor = get(c,'edgecolor');
          if isa(inputcolor,'char')
            inputcolor=[ 1 1 1 ];
          end
          slowset(c,'edgecolor',uisetcolor(inputcolor),d.animincrement);
        end
       case 'slicecontourlinewidth'
        [a s]=getarrowslice;
        c = getappdata(s,'contour');
        if ~isempty(c)
          if isa(p2,'char')
            slowset(c,'linewidth',str2num(p2),d.animincrement);
          else
            slowset(c,'linewidth',p2,d.animincrement);
          end
        end
        %
        % Menu All Slices
        %
       case 'allfacet'
        s=allSlices;
        set(s,'facec','flat','edgec','k');
        textureizeslice(s,'off');
       case 'allflat'
        s=allSlices;
        set(s,'facec','flat','edgec','none');
        textureizeslice(s,'off');
       case 'allinterp'
        s=allSlices;
        set(s,'facec','interp','edgec','none');
        textureizeslice(s,'off');
       case 'alltex'
        s=allSlices;
        set(s,'facec','texturemap','edgec','none');
        textureizeslice(s,'on');
       case 'allnone'
        s=allSlices;
        set(s,'facec','none','edgec','none');
        textureizeslice(s,'off');
       case 'alltnone'
        s=allSlices;
        set(s,'facea',1);
        textureizeslice(s,'off');
       case 'alltp5'
        s=allSlices;
        set(s,'facea',.5);
        textureizeslice(s,'off');
       case 'alltflat'
        s=allSlices;
        set(s,'facea','flat');
        textureizeslice(s,'off');
       case 'alltinterp'
        s=allSlices;
        set(s,'facea','interp');
        textureizeslice(s,'off');
       case 'allttex'
        s=allSlices;
        set(s,'facea','texturemap');
        textureizeslice(s,'on');
        %
        % Menu Defaults callbacks
        %
       case	'defaultfaceted'
        d.defcolor='faceted';
       case	'defaultflat'
        d.defcolor='flat';
       case	'defaultinterp'
        d.defcolor='interp';
       case	'defaulttexture'
        d.defcolor='texture';
        if strcmp(d.defalpha,'flat') || strcmp(d.defalpha,'interp')
          d.defalpha='texture';
        end
       case	'defaultinterp'
        d.defcolor='none';
       case	'defaulttransnone'
        d.defalpha='none';
       case	'defaulttransflat'
        d.defalpha='flat';
       case	'defaulttransinterp'
        d.defalpha='interp';
       case	'defaulttranstexture'
        d.defalpha='texture';
        d.defcolor='texture';
       case      'defaultlightflat'
        d.deflight='flat';
       case      'defaultlightsmooth'
        d.deflight='smooth';
       case 'defaultcontoursmooth'
        d.defaultcontoursmooth='on';
       case 'defaultcontourflat'
        d.defcontourcolor='flat';
       case 'defaultcontourinterp'
        d.defcontourcolor='interp';
       case 'defaultcontourblack'
        d.defcontourcolor='black';
       case 'defaultcontourwhite'
        d.defcontourcolor='white';
       case 'defaultcontourlinewidth'
        if isa(p2,'char')
          d.defcontourlinewidth=str2num(p2);
        else
          d.defcontourlinewidth=p2;
        end
        %
        % Camera toolbar Toggling
        %
       case 'cameratoolbar'
        cameratoolbar('Toggle');
       case 'annotationtoolbar'
        if propcheck(d.toolbar,'visible','on')
          set(d.toolbar,'vis','off');
        else
          set(d.toolbar,'vis','on');
        end
        %
        % Controller Preferences
        %
       case 'controlalpha'
        val=str2num(p2);
        iso=findobj(d.axiso,'type','image');
        if val == 0
          set([d.pxx d.pxy d.pxz iso],'visible','off');
        else
          set([d.pxx d.pxy d.pxz iso],'visible','on');
          slowset([d.pxx d.pxy d.pxz] , 'facealpha',val,d.animincrement);
          slowset(iso,'alphadata',val,d.animincrement);
        end
       case 'toggleanimation'
        if d.animincrement == 0
          d.animincrement = 10;
        else
          d.animincrement = 0;
        end
       case 'controllabels'
        l = get(d.axx,'xticklabel');
        if isempty(l)
          set([d.axx d.axiso],'xticklabelmode','auto');
          set([d.axy d.axz],'yticklabelmode','auto');
        else
          set([d.axx d.axiso],'xticklabel',[]);
          set([d.axy d.axz],'yticklabel',[]);
        end
       case 'controlvisible'
        objs=findobj([d.axiso d.axx d.axy d.axz]);
        if strcmp(get(d.axx,'visible'),'on')
          set(objs,'visible','off');
          set(d.axmain,'pos',[.1 .1 .9 .8]);
        else
          set(objs,'visible','on');
          set(d.axmain,'pos',[.2  .2 .6 .6]);
        end
        %
        % UICONTROL callbacks
        %
       case 'colormap'
        str=get(gcbo,'string');
        val=str{get(gcbo,'value')};
        size(val);
        if strcmp(val,'custom')
          cmapeditor
        else
          if strcmp(val, 'rand')
              cm = get(gcf, 'colormap');
              val = feval(@rand, size(cm));
          else
              val = feval(val);
          end            
          slowset(gcf,'colormap',val,d.animincrement);
        end
       case 'alphamap'
        str=get(gcbo,'string');
        str=str{get(gcbo,'value')};
        if strcmp(str, 'rand')
            val=rand(size(alphamap));
        else
            val=alphamap(str);
        end
        slowset(gcf,'alphamap',val,d.animincrement);
        %
        % Commands
        %
       case 'copy'
        copyobj(gca,figure);set(gca,'pos',[.1 .1 .9 .8]);
       case 'print'
        newf=figure('visible','off','renderer',get(gcf,'renderer'));
        copyobj(d.axmain,newf);
        set(gca,'pos',[.1 .1 .9 .8])
        printdlg(newf);
        close(newf);
       otherwise
        error('Bad slice-o-matic command.');
      end
    catch
      disp(get(0,'errormessage'));
    end
    setappdata(gcf,'sliceomatic',d);
  else
    disp('Sliceomatic data must be DOUBLE');
  end

function dragprep(arrowtodrag)

  arrows=findall(gcf,'tag','sliceomaticarrow');

  pushset(arrows,'facecolor',[1 0 0]);
  pushset(arrows,'facealpha',.2);

  pushset(arrowtodrag,'facecolor',[0 1 0]);
  pushset(arrowtodrag,'facealpha',.7);

  slices=allSlices;

  for i=1:length(slices)
    fa=get(slices(i),'facea');
    if isa(fa,'double') && fa>.3
      pushset(slices(i),'facealpha',.3);
      pushset(slices(i),'edgecolor','n');
    else
      pushset(slices(i),'facealpha',fa);
      pushset(slices(i),'edgecolor',get(slices(i),'edgec'));
    end
  end

  isosurfs=allIsos;

  for i=1:length(isosurfs)
    fa=get(isosurfs(i),'facea');
    if isa(fa,'double') && fa>.3
      pushset(isosurfs(i),'facealpha',.3);
      pushset(isosurfs(i),'edgecolor','n');
    else
      pushset(isosurfs(i),'facealpha',fa);
      pushset(isosurfs(i),'edgecolor',get(isosurfs(i),'edgec'));
    end
    cap=getappdata(isosurfs(i),'isosurfacecap');
    if ~isempty(cap)
      pushset(cap,'visible','off');
    end
  end

  ss=getappdata(arrowtodrag,'arrowslice');

  if isempty(ss)
    ss=getappdata(arrowtodrag,'arrowiso');
  end

  popset(ss,'facealpha');
  popset(ss,'edgecolor');

  pushset(gcf,'windowbuttonupfcn','sliceomatic up');
  pushset(gcf,'windowbuttonmotionfcn','sliceomatic motion');

  % Doing this makes the tip invisible when visible is on.
  showarrowtip(arrowtodrag);
  
function dragfinis(arrowtodrag)

  arrows=findall(gcf,'tag','sliceomaticarrow');

  popset(arrowtodrag,'facecolor');
  popset(arrowtodrag,'facealpha');

  popset(arrows,'facecolor');
  popset(arrows,'facealpha');

  ss=getappdata(arrowtodrag,'arrowslice');
  if isempty(ss)
    ss=getappdata(arrowtodrag,'arrowiso');
  end

  % These pushes are junk which will be undone when all slices or
  % isosurfs are reset below.
  pushset(ss,'facealpha',1);
  pushset(ss,'edgecolor','k');

  slices=allSlices;

  if ~isempty(slices)
    popset(slices,'facealpha');
    popset(slices,'edgecolor');
  end

  isosurfs=allIsos;

  if ~isempty(isosurfs)
    popset(isosurfs,'facealpha');
    popset(isosurfs,'edgecolor');
  end

  d=getappdata(gcf,'sliceomatic');
  
  if isnan(d.xmesh)==1
    for i=1:length(isosurfs)
      cap=getappdata(isosurfs(i),'isosurfacecap');
      if ~isempty(cap)
        popset(cap,'visible');
        localisocaps(isosurfs(i),cap);
      end
      if getappdata(isosurfs(i), 'reduced')
        setappdata(isosurfs(i),'reduced',0);
        localisosurface({},d.data,d.smooth,...
                        getappdata(isosurfs(i),'isosurfacevalue'),...
                        isosurfs(i));
      end
    end
  else
    for i=1:length(isosurfs)
      cap=getappdata(isosurfs(i),'isosurfacecap');
      if ~isempty(cap)
        popset(cap,'visible');
        localisocaps(isosurfs(i),cap);
      end
      if getappdata(isosurfs(i), 'reduced')
        setappdata(isosurfs(i),'reduced',0);
        realvolume={d.xmesh d.ymesh d.zmesh};
        localisosurface(realvolume,d.data,d.smooth,...
                        getappdata(isosurfs(i),'isosurfacevalue'),...
                        isosurfs(i));
      end
    end
  end

  popset(gcf,'windowbuttonupfcn');
  popset(gcf,'windowbuttonmotionfcn');

  showarrowtip([]);
  
  % Make sure whatever buttonupfcn on the figure is run now to "turn
  % off" whatever was going on before we got our callback on the
  % arrow.

  buf = get(gcf,'windowbuttonupfcn');
  if ~strcmp(buf,'')
    eval(buf);
  end

function movetipforarrow(arrow, ax, value, position, va, ha)
% Setup the current data tip for a slice arrow, and show it's
% control value
  
  tipdata.parentaxes = ax;
  tipdata.value = value;
  tipdata.position = position;
  tipdata.verticalalign = va;
  tipdata.horizontalalign = ha;
  
  setappdata(arrow, 'tipdata', tipdata);
  
  showarrowtip(arrow);
% Put it onto d.axisiso so that
% it always appears on top.
%set(t,'parent',d.axiso);

function p=arrow(parent,dir,pos)

%   21012    21012      12345     12345
% 5  *-*   5   *     2   *     2   *  
% 4  | |   4  / \    1 *-*\    1  /*-*
% 3 ** **  3 ** **   0 |   *   0 *   |
% 2  \ /   2  | |   -1 *-*/   -1  \*-*
% 1   *    1  *-*   -2   *    -2   *  

  switch dir
   case 'down'
    pts=[ 0 1; -2 3; -1 3; -1 5; 1 5; 1 3; 2 3 ];
    mp = 'SOM leftright';
   case 'up'
    pts=[ 0 5; 2 3; 1 3; 1 1; -1 1; -1 3; -2 3; ];
    mp = 'SOM leftright';
   case 'right'
    pts=[ 5 0; 3 -2; 3 -1; 1 -1; 1 1; 3 1; 3 2 ];
    mp = 'SOM topbottom';
   case 'left'
    pts=[ 1 0; 3 2; 3 1; 5 1; 5 -1; 3 -1; 3 -2 ];
    mp = 'SOM topbottom';
  end

  f=[1 2 7; 3 4 5; 3 5 6 ];

  % Modify the arrows to look good no matter what
  % the data aspect ratio may be.
  if pos(1)
    lim=get(parent,'xlim');
    fivep=abs(lim(1)-lim(2))/15/5;
    pts(:,1)=pts(:,1)*fivep+pos(1);
  elseif pos(2)
    lim=get(parent,'ylim');
    fivep=abs(lim(1)-lim(2))/15/5;
    pts(:,2)=pts(:,2)*fivep+pos(2);
  end

  % Create the patches, and add app data to them to remember what
  % They are associated with.
  p(1)=patch('vertices',pts,'faces',1:size(pts,1),'facec','n','edgec','k',...
             'linewidth',2,'hittest','off',...
             'parent',parent);
  p(2)=patch('vertices',pts,'faces',f,'facec','g','facea',.5,'edgec','n',...
             'parent',parent,'tag','sliceomaticarrow');
  setappdata(p(2),'arrowcenter',pos);
  setappdata(p(2),'arrowedge',p(1));
  setappdata(p(2),'motionpointer',mp);

  
function p=localisocaps(isosurface,isocap)
% Isocap management
  
% Get relevant info from the isosurface.
  if nargin<2 || ~strcmp(get(isocap,'visible'),'off')
    d=getappdata(gcf,'sliceomatic');
    data=getappdata(isosurface,'isosurfacedata');
    if isnan(d.xmesh)==1
      caps=isocaps(data,getappdata(isosurface,'isosurfacevalue'));
    else
      caps=isocaps(d.xmesh,d.ymesh,d.zmesh,data,getappdata(isosurface,'isosurfacevalue'));
    end
  end

  if nargin==2
    if ~strcmp(get(isocap,'visible'),'off')
      set(isocap,caps);
    end
    p=isocap;
  else
    p=patch(caps,'edgecolor','none','facecolor','flat',...
            'facelighting','none',...
            'tag','sliceomaticisocap');

    setappdata(p,'isosurface',isosurface);
    setappdata(isosurface,'isosurfacecap',p);
    
    d=getappdata(gcf,'sliceomatic');
    
    switch d.defcolor
     case 'faceted'
      set(p,'facec','flat','edgec','black');
     case 'flat'
      set(p,'facec','flat','edgec','none');
     case 'interp'
      set(p,'facec','interp','edgec','none');
     case 'texture'
      set(p,'facec','flat','edgec','none');
     case 'none'
      set(p,'facec','none','edgec','none');
    end
    switch d.defalpha
     case 'none'
      set(p,'facea',1);
     case 'flat'
      set(p,'facea','flat');
     case 'interp'
      set(p,'facea','interp');
     case 'texture'
      set(p,'facea','flat');
    end    
  end


function p=localisosurface(volume, data, datanormals, value, oldiso)
% Isosurface management
  
  pushset(gcf, 'pointer','watch');

  d=getappdata(gcf,'sliceomatic');
  
  fv = isosurface(volume{:},data, value);
  
  clim=get(gca,'clim');
  cmap=get(gcf,'colormap');
  clen=clim(2)-clim(1);
  idx=fix((value-clim(1))*length(cmap)/clen)+1;

  if nargin==5
    try
      set(oldiso,fv,'facecolor',cmap(idx,:));
    catch
      set(oldiso,fv,'facecolor','none');
    end
    p=oldiso;
    cap=getappdata(p,'isosurfacecap');
    if ~isempty(cap)
      localisocaps(p,cap);
    end
  else
    if isnan(d.xmesh)==1
      p=patch(fv,'edgecolor','none','facecolor',cmap(idx,:),...
              'tag', 'sliceomaticisosurface');
    else
      p=patch(fv,'edgecolor','none','facecolor',cmap(idx,:),...
              'tag', 'sliceomaticisosurface');
    end
    % d=getappdata(gcf,'sliceomatic');
    switch d.deflight
     case 'flat'
      set(p,'facelighting','flat');
     case 'smooth'
      set(p,'facelighting','phong');
    end
    setappdata(p,'isosurfacecap',[]);
  end

  setappdata(p,'isosurfacevalue',value);
  setappdata(p,'isosurfacedata',data);

  reducepatch(p,10000);
  isonormals(volume{:},datanormals,p);
  
  popset(gcf,'pointer');

function s=localslice(data, X, Y, Z, oldslice)
% Slice Management.  Uses specialized slicomatic slices, not slices
% created with the SLICE command.

  s=[];
  d=getappdata(gcf,'sliceomatic');

  ds=size(data);

  if ~isempty(X)
    xi=round(X);
    if isnan(d.xmesh) == 1
      if xi > 0 && xi <= ds(2)
        cdata=reshape(data(:,xi,:),ds(1),ds(3));
        [xdata ydata zdata]=meshgrid(xi,1:ds(1),1:ds(3));
        st = 'X';
      else
        return
      end
    else
      if isequal(d.xdir,'reverse')==1
        locate_xi=histc(xi,flipdim(d.xmesh,2));
        slice_number=find(locate_xi);
        slice_number=length(d.xmesh)-slice_number+1;
      else 
        locate_xi=histc(xi,d.xmesh);
        slice_number=find(locate_xi);
      end
      if ~isempty(slice_number) && slice_number > 0 && slice_number <= ds(2)
        cdata=reshape(data(:,slice_number,:),ds(1),ds(3));
        [xdata ydata zdata]=meshgrid(X,d.ymesh,d.zmesh);
        st = 'X';
      else
        return
      end
    end
    
  elseif ~isempty(Y)
    yi=round(Y);
    if isnan(d.ymesh) == 1
      if yi > 0 && yi <= ds(1)
        cdata=reshape(data(yi,:,:),ds(2),ds(3));
        [xdata ydata zdata]=meshgrid(1:ds(2),yi,1:ds(3));
        st = 'Y';
      else
        return    
      end
    else
      if isequal(d.ydir,'reverse')==1
        locate_yi=histc(yi,flipdim(d.ymesh,2));
        slice_number=find(locate_yi);
        slice_number=length(d.ymesh)-slice_number+1;
      else 
        locate_yi=histc(yi,d.ymesh);
        slice_number=find(locate_yi);
      end
      if ~isempty(slice_number) && slice_number > 0 && slice_number <= ds(1)
        cdata=reshape(data(slice_number,:,:),ds(2),ds(3));
        [xdata ydata zdata]=meshgrid(d.xmesh,Y,d.zmesh);
        st = 'Y';
      else
        return
      end
    end
    
  elseif ~isempty(Z)
    zi=round(Z);
    if isnan(d.zmesh) == 1
      if zi > 0 && zi <= ds(3)
        cdata=reshape(data(:,:,zi),ds(1),ds(2));
        [xdata ydata zdata]=meshgrid(1:ds(2),1:ds(1),zi);
        st = 'Z';
      else
        return
      end
    else
      if isequal(d.zdir,'reverse')==1
        locate_zi=histc(zi,flipdim(d.zmesh,2));
        slice_number=find(locate_zi);
        slice_number=length(d.zmesh)-slice_number+1;
      else 
        locate_zi=histc(zi,d.zmesh);
        slice_number=find(locate_zi);
      end
      if ~isempty(slice_number) && slice_number > 0 && slice_number <= ds(3)
        cdata=reshape(data(:,:,slice_number),ds(1),ds(2));
        [xdata ydata zdata]=meshgrid(d.xmesh,d.ymesh,Z);
        st = 'Z';
      else
        return
      end
    end
  else
    error('Nothing was passed into LOCALSLICE.');
  end

  cdata=squeeze(cdata);
  xdata=squeeze(xdata);
  ydata=squeeze(ydata);
  zdata=squeeze(zdata);

  if nargin == 5
    % Recycle the old slice
    set(oldslice,'cdata',cdata,'alphadata',cdata, 'xdata',xdata, ...
                 'ydata',ydata, 'zdata',zdata);
    s=oldslice;
    %delete(news);
    if propcheck(s,'facec','texturemap')
      textureizeslice(s,'on');
    end
    setappdata(s,'slicetype',st);
  else
    % setup the alphadata
    news=surface('cdata',cdata,'alphadata',cdata, 'xdata',xdata, ...
                 'ydata',ydata, 'zdata',zdata);
    set(news,'alphadata',cdata,'alphadatamapping','scaled','tag','sliceomaticslice',...
             'facelighting','none',...
             'uicontextmenu',d.uic);
    s=news;
    setappdata(s,'slicetype',st);
    switch d.defcolor
     case 'faceted'
      set(s,'facec','flat','edgec','k');
     case 'flat'
      set(s,'facec','flat','edgec','n');
     case 'interp'
      set(s,'facec','interp','edgec','n');
     case 'texture'
      set(s,'facec','texture','edgec','n');
    end
    switch d.defalpha
     case 'none'
      set(s,'facea',1);
     case 'flat'
      set(s,'facea','flat');
     case 'interp'
      set(s,'facea','interp');
     case 'texture'
      set(s,'facea','texture');
    end    
  end
  
  contour = getappdata(s,'contour');
  if ~isempty(contour)
    try 
      levels = getappdata(s, 'contourlevels');
      if isempty(levels)~=1
        localcontour(s,contour,levels);
      else
        localcontour(s, contour);
      end
    catch
      localcontour(s, contour);
    end
  end
  

function textureizeslice(slice,onoff)
% Convert a regular slice into a texture map slice, or a texture
% slice into a regular slice.
  
  for k=1:prod(size(slice))

    d=getappdata(slice(k),'textureoptimizeations');

    switch onoff
     case 'on'
      d.xdata=get(slice(k),'xdata');
      d.ydata=get(slice(k),'ydata');
      d.zdata=get(slice(k),'zdata');
      setappdata(slice(k),'textureoptimizeations',d);
      if max(size(d.xdata)==1)
        nx=[d.xdata(1) d.xdata(end)];
      else
        nx=[d.xdata(1,1)   d.xdata(1,end);
            d.xdata(end,1) d.xdata(end,end)];
      end
      if max(size(d.ydata)==1)
        ny=[d.ydata(1) d.ydata(end)];
      else
        ny=[d.ydata(1,1)   d.ydata(1,end);
            d.ydata(end,1) d.ydata(end,end)];
      end
      if max(size(d.zdata)==1)
        nz=[d.zdata(1) d.zdata(end)];
      else
        nz=[d.zdata(1,1)   d.zdata(1,end);
            d.zdata(end,1) d.zdata(end,end)];
      end
      set(slice(k),'xdata',nx, 'ydata', ny, 'zdata', nz,...
                   'facec','texturemap');
      if ischar(get(slice(k),'facea'))
        set(slice(k),'facea','texturemap');
      end
      if ischar(get(slice(k),'facec'))
        set(slice(k),'facec','texturemap');
      end
     case 'off'
      if ~isempty(d)
        set(slice(k),'xdata',d.xdata,'ydata',d.ydata,'zdata',d.zdata);
        setappdata(slice(k),'textureoptimizeations',[]);
      end
      if ischar(get(slice(k),'facea')) && strcmp(get(slice(k),'facea'),'texturemap')
        set(slice(k),'facea','flat');
      end
      if ischar(get(slice(k),'facec')) && strcmp(get(slice(k),'facec'),'texturemap')
        set(slice(k),'facec','flat');
      end
    end
  end


function localcontour(slice,oldcontour,levels)
% Create a contour on SLICE
% When OLDCONTROUR, recycle that contour patch.
% This does not use the CONTOURSLICE command, but instead uses a
% specialized slice created for sliceomantic.

  d=getappdata(gcf,'sliceomatic');
  
  cdata = get(slice,'cdata');
  st = getappdata(slice,'slicetype');

  % Calculate the new contour for CDATA's values.
  if nargin < 3
    if isnan(d.zmesh)==1
      c = contourc(double(cdata));
    else
      switch st
       case 'X'
        c = contours(d.zmesh,d.ymesh,cdata);
       case 'Y'
        c = contours(d.zmesh,d.xmesh,cdata);
       case'Z'
        c = contours(d.xmesh,d.ymesh,cdata);
      end
    end
  else
    if isnan(d.zmesh)==1
      c = contourc(double(cdata),levels);
    else
      switch st
       case 'X'
        c = contours(d.zmesh,d.ymesh,cdata,levels);
       case 'Y'
        c = contours(d.zmesh,d.xmesh,cdata,levels);
       case 'Z'
        c = contours(d.xmesh,d.ymesh,cdata,levels);
      end
    end
  end

  newvertices = [];
  newfaces = {};
  longest = 1;
  cdata = [];
  
  limit = size(c,2);
  i = 1;
  while(i < limit)
    z_level = c(1,i);
    npoints = c(2,i);
    nexti = i+npoints+1;

    xdata = c(1,i+1:i+npoints);
    ydata = c(2,i+1:i+npoints);

    switch st
     case 'X'
      xv = get(slice,'xdata');
      lzdata = xv(1,1) + 0*xdata;
      vertices = [lzdata.', ydata.', xdata.'];
     case 'Y'
      yv = get(slice,'ydata');
      lzdata = yv(1,1) + 0*xdata;
      vertices = [ydata.', lzdata.', xdata.'];
     case 'Z'
      zv = get(slice,'zdata');
      lzdata = zv(1,1) + 0*xdata;
      vertices = [xdata.', ydata.', lzdata.'];
    end    
    
    faces = 1:length(vertices);
    faces = faces + size(newvertices,1);

    longest=max(longest,size(faces,2));
    
    newvertices = [ newvertices ; vertices ];
    newfaces{end+1} = faces;
    
    tcdata =  (z_level + 0*xdata).';
    
    cdata = [ cdata; tcdata ]; % need to be same size as faces
    
    i = nexti;
  end

  % Glom a NAN on the end for loop-breaking
  newvertices = [ newvertices ; nan nan nan ];
  cdata = [ cdata ; nan ];
  
  vertmax = size(newvertices,1);
  
  % Fix up FACES, which is a cell array.
  faces = [];
  for i = 1:size(newfaces,2)
    faces = [ faces;
              newfaces{i} ones(1,longest-size(newfaces{i},2))*vertmax vertmax ];
  end
  
  if isempty(oldcontour)
    oldcontour = patch('facecolor','none', 'edgecolor',d.defcontourcolor,...
                       'linewidth',d.defcontourlinewidth);
    try
      set(oldcontour,'linesmoothing',d.defcontoursmooth);
    catch
    end
    setappdata(slice,'contour',oldcontour);
  end

  set(oldcontour,'vertices',newvertices,...
                 'faces',faces,...
                 'facevertexcdata',cdata);

function ss=allSlices
  ss=findobj(gcf,'type','surface','tag','sliceomaticslice');

function ss=allIsos
  ss=findobj(gcf,'type','patch','tag','sliceomaticisosurface');

function ss=allCaps
  ss=findobj(gcf,'type','patch','tag','sliceomaticisocap');

function working(onoff)

  ax=getappdata(gcf,'workingaxis');

  if isempty(ax)
    ax=axes('units','norm','pos',[.3 .4 .4 .2],...
            'box','on','ytick',[],'xtick',[],...
            'xlim',[-1 1],'ylim',[-1 1],...
            'color','none','handlevis','off');
    text('parent',ax,'string','Working...','fontsize',64,...
         'pos',[0 0], ...
         'horizontalalignment','center',...
         'verticalalignment','middle',...
         'erasemode','xor');
    setappdata(gcf,'workingaxis',ax);
  end

  disp(['Working...' onoff]);
  set([ax get(ax,'children')],'vis',onoff);

function activelabel(label, string)
% ACTIVELABEL(LABEL, STRING) - Create a label on GCA which is
%     active.  LABEL is the property of GCA whose label you are
%     setting.  STRING is the initial text string for the label.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  l = get(gca,label);
  
  set(l,'string',string);
  set(l,'buttondownfcn',@activelabelbuttondown);
  
function activelabelbuttondown(obj, action)
% Callback when one of our active labels is clicked on.
  
  set(obj,'edit','on');
  
function outd = figmenus(d)
% Set up sliceomatic's gui menus within structure D

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

%%% 2/18/05 RAB part in 3 parts %%% 

% Main Figure Menu
  set(gcf,'menubar','none');
  
  % File menu
  d.filemenu = uimenu(gcf,'label','File');
  d.fcopy = uimenu(d.filemenu, 'label', 'Copy figure','callback', 'sliceomatic copy');
  d.fprint  = uimenu(d.filemenu,'label','Print...','callback','sliceomatic print');

  %%% start patch 1of3 RAB 2/18/05 %%%
  d.fsaveprefs = uimenu(d.filemenu,'label','Save preferences','callback',@SavePrefs);
  %%% end patch 1of3 RAB 2/18/05 %%%

  % How do get these props onto the print figure?
  %d.fprints = uimenu(d.filemenu,'label','Print Setup...','callback','printdlg -setup');
  % ---
  d.fexit = uimenu(d.filemenu, 'label', 'Close','callback','closereq',...
                   'separator','on');

  % Controls Menu
  d.defcontrols = uimenu(gcf,'label','Controls', 'callback',@controlmenu);
  if exist('uitoolfactory') == 2
    d.anntoolbar = uimenu(d.defcontrols,'label','Annotations toolbar','callback', 'sliceomatic annotationtoolbar');
  end
  d.camtoolbar = uimenu(d.defcontrols,'label','Camera toolbar','callback', 'sliceomatic cameratoolbar');
  d.dcalpha = uimenu(d.defcontrols,'label','Controls Transparency');
  d.dcalpha1= uimenu(d.dcalpha,'label','1','callback','sliceomatic controlalpha 1');
  d.dcalpha8= uimenu(d.dcalpha,'label','.8','callback','sliceomatic controlalpha .8');
  d.dcalpha6= uimenu(d.dcalpha,'label','.6','callback','sliceomatic controlalpha .6');
  d.dcalpha5= uimenu(d.dcalpha,'label','.5','callback','sliceomatic controlalpha .5');
  d.dcalpha4= uimenu(d.dcalpha,'label','.4','callback','sliceomatic controlalpha .4');
  d.dcalpha2= uimenu(d.dcalpha,'label','.2','callback','sliceomatic controlalpha .2');
  d.dcalpha0= uimenu(d.dcalpha,'label','0','callback','sliceomatic controlalpha 0');
  d.dcanimstep = uimenu(d.defcontrols,'label','Animation','callback', 'sliceomatic toggleanimation');
  d.dclabels= uimenu(d.defcontrols','label','Tick Labels','callback','sliceomatic controllabels');
  d.dcvis   = uimenu(d.defcontrols','label','Visible','callback','sliceomatic controlvisible');
  % d.dsetrange= uimenu(d.defcontrols','label','Set Range','callback','@setvolumerange');
  % d.dcslice = uimenu(d.defcontrols,'label','Slice Controls','callback','sliceomatic useslicecontrols');
  % d.dciso   = uimenu(d.defcontrols,'label','Iso Surface Control','callback','sliceomatic useisocontrols','separator','on');

  % Remove this once we have more controls to enable and disable.
  %  set(d.defcontrols,'vis','off');
  
  % Default for new slices menu
  d.defmenu = uimenu(gcf,'label','Object_Defaults', 'callback', @defaultmenu);
  d.dfacet  = uimenu(d.defmenu,'label','Slice Color Faceted','callback','sliceomatic defaultfaceted');
  d.dflat   = uimenu(d.defmenu,'label','Slice Color Flat',   'callback','sliceomatic defaultflat');
  d.dinterp = uimenu(d.defmenu,'label','Slice Color Interp', 'callback','sliceomatic defaultinterp');
  d.dtex    = uimenu(d.defmenu,'label','Slice Color Texture','callback','sliceomatic defaulttexture');
  d.dcnone  = uimenu(d.defmenu,'label','Slice Color None','callback','sliceomatic defaultcolornone');
  d.dtnone  = uimenu(d.defmenu,'label','Slice Transparency None','callback','sliceomatic defaulttransnone','separator','on');
  d.dtflat  = uimenu(d.defmenu,'label','Slice Transparency Flat','callback','sliceomatic defaulttransflat');
  d.dtinterp= uimenu(d.defmenu,'label','Slice Transparency Interp','callback','sliceomatic defaulttransinterp');
  d.dttex   = uimenu(d.defmenu,'label','Slice Transparency Texture','callback','sliceomatic defaulttranstexture');
  d.dlflat  = uimenu(d.defmenu,'label','IsoSurface Lighting Flat','callback','sliceomatic defaultlightflat','separator','on');
  d.dlsmooth= uimenu(d.defmenu,'label','IsoSurface Lighting Smooth','callback','sliceomatic defaultlightsmooth');
  %d.dcsmooth= uimenu(d.defmenu,'label','Contour Line Smoothing','callback','sliceomatic defaultcontoursmooth');
  d.dcflat  = uimenu(d.defmenu,'label','Contour Color Flat',   'callback','sliceomatic defaultcontourflat','separator','on');
  d.dcinterp= uimenu(d.defmenu,'label','Contour Color Interp', 'callback','sliceomatic defaultcontourinterp');
  d.dcblack = uimenu(d.defmenu,'label','Contour Color Black',  'callback','sliceomatic defaultcontourblack');
  d.dcwhite = uimenu(d.defmenu,'label','Contour Color White',  'callback','sliceomatic defaultcontourwhite');
  d.dclinew = uimenu(d.defmenu,'label','Contour Line Width');
  d.dcl1    = uimenu(d.dclinew,'label','1','callback','sliceomatic defaultcontourlinewidth 1');
  d.dcl2    = uimenu(d.dclinew,'label','2','callback','sliceomatic defaultcontourlinewidth 2');
  d.dcl3    = uimenu(d.dclinew,'label','3','callback','sliceomatic defaultcontourlinewidth 3');
  d.dcl4    = uimenu(d.dclinew,'label','4','callback','sliceomatic defaultcontourlinewidth 4');
  d.dcl5    = uimenu(d.dclinew,'label','5','callback','sliceomatic defaultcontourlinewidth 5');
  d.dcl6    = uimenu(d.dclinew,'label','6','callback','sliceomatic defaultcontourlinewidth 6');
  
  d.defcolor='texture';
  d.defalpha='texture';
  d.deflight='smooth';
  d.defcontourcolor='black';
  d.defcontourlinewidth=1;
  % This exposes an unpleasant R14 bug
  d.defcontoursmooth='off';

  % investigate hardware opengl.
  inc = 0;
  try
    od = opengl('data');
    if isfield(od,'Software')
      % R14 version of MATLAB
      if ~od.Software
        inc = 10;
      end
    else
      % Older version of MATLAB
      if ~(strcmp(od.Renderer,'Mesa X11') || ...
           strcmp(od.Renderer, 'GDI Generic'))
        inc = 10;
      end
    end
  end
  
  d.animincrement=inc;
  
  %%% start patch 2of3 RAB 2/18/05 %%%
  d = OverrideStickyUserPreferences(d);
  %%% end patch 2of3 RAB 2/18/05 %%%
  
  % Set props for all slices menu
  d.allmenu = uimenu(gcf,'label','AllSlices');
  uimenu(d.allmenu,'label','Color Faceted','callback','sliceomatic allfacet');
  uimenu(d.allmenu,'label','Color Flat','callback','sliceomatic allflat');
  uimenu(d.allmenu,'label','Color Interp','callback','sliceomatic allinterp');
  uimenu(d.allmenu,'label','Color Texture','callback','sliceomatic alltex');
  uimenu(d.allmenu,'label','Color None','callback','sliceomatic allnone');
  uimenu(d.allmenu,'label','Transparency None','callback','sliceomatic alltnone','separator','on');
  uimenu(d.allmenu,'label','Transparency .5','callback','sliceomatic alltp5');
  uimenu(d.allmenu,'label','Transparency Flat','callback','sliceomatic alltflat');
  uimenu(d.allmenu,'label','Transparency Interp','callback','sliceomatic alltinterp');
  uimenu(d.allmenu,'label','Transparency Texture','callback','sliceomatic allttex');

  % Setup Help style options
  d.helpmenu = uimenu(gcf,'label','Help');
  uimenu(d.helpmenu,'label','Help','callback','doc sliceomatic/sliceomatic');
  uimenu(d.helpmenu,'label','Check for Updates','callback','web http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=764&objectType=FILE');
  uimenu(d.helpmenu,'label','About Author','callback','web http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objectId=803709&objectType=author');
  
  % Context Menus
  % Slice Context Menu
  d.uic=uicontextmenu('callback', @slicecontextmenu);
  d.vistog = uimenu(d.uic,'label','Visible','callback','sliceomatic togglevisible');
  d.uicdelete = uimenu(d.uic,'label','Delete','callback','sliceomatic deleteslice');
  d.smcolorm  = uimenu(d.uic,'label','Color','separator','on');
  d.smfacet   = uimenu(d.smcolorm,'label','Color Faceted','callback','sliceomatic setfaceted');
  d.smflat    = uimenu(d.smcolorm,'label','Color Flat','callback','sliceomatic setflat');
  d.sminterp  = uimenu(d.smcolorm,'label','Color Interp','callback','sliceomatic setinterp');
  d.smtex     = uimenu(d.smcolorm,'label','Color Texture','callback','sliceomatic settexture');
  d.smnone    = uimenu(d.smcolorm,'label','Color None','callback','sliceomatic setnone');
  d.smtransm  = uimenu(d.uic,'label','Transparency');
  d.smtnone   = uimenu(d.smtransm,'label','Transparency None','callback','sliceomatic setalphanone');
  d.smtp5     = uimenu(d.smtransm,'label','Transparency .5','callback','sliceomatic setalphapoint5');
  d.smtflat   = uimenu(d.smtransm,'label','Transparency Flat','callback','sliceomatic setalphaflat');
  d.smtinterp = uimenu(d.smtransm,'label','Transparency Interp','callback','sliceomatic setalphainterp');
  d.smttex    = uimenu(d.smtransm,'label','Transparency Texture','callback','sliceomatic setalphatexture');
  d.smcontour = uimenu(d.uic,'label','Add Contour','separator','on');
  d.smcont0   = uimenu(d.smcontour,'label','Auto (Slice)','callback','sliceomatic slicecontour');
  d.smcont0v  = uimenu(d.smcontour,'label','Auto (Volume)','callback','sliceomatic slicecontourfullauto');
  d.smcont1   = uimenu(d.smcontour,'label','Select Levels','callback','sliceomatic slicecontour_select','separator','on');
  d.smcsetauto= uimenu(d.uic,'label','Set Auto Levels (Slice)','callback','sliceomatic slicecontour_setauto');
  d.smcsetav  = uimenu(d.uic,'label','Set Auto Levels (Volume)','callback','sliceomatic slicecontour_setfullauto');
  d.smclevels = uimenu(d.uic,'label','Set Levels','callback','sliceomatic slicecontour_setlevels');
  d.smrcontour= uimenu(d.uic,'label','Remove Contour','callback','sliceomatic deleteslicecontour');
  d.smccm     = uimenu(d.uic,'label','Contour Colors');
  d.smcflat   = uimenu(d.smccm,'label','Contour Flat','callback','sliceomatic slicecontourflat');
  d.smcinterp = uimenu(d.smccm,'label','Contour Interp','callback','sliceomatic slicecontourinterp');
  d.smcblack  = uimenu(d.smccm,'label','Contour Black','callback','sliceomatic slicecontourblack');
  d.smcwhite  = uimenu(d.smccm,'label','Contour White','callback','sliceomatic slicecontourwhite');
  d.smccolor  = uimenu(d.smccm,'label','Contour Color','callback','sliceomatic slicecontourcolor');
  d.smcsmooth = uimenu(d.uic,'visible','off','label','Smooth Contour Lines','callback','sliceomatic slicecontoursmooth');
  d.smclinew  = uimenu(d.uic,'label','Contour Line Width');
  d.smcl1     = uimenu(d.smclinew,'label','1','callback','sliceomatic slicecontourlinewidth 1');
  d.smcl2     = uimenu(d.smclinew,'label','2','callback','sliceomatic slicecontourlinewidth 2');
  d.smcl3     = uimenu(d.smclinew,'label','3','callback','sliceomatic slicecontourlinewidth 3');
  d.smcl4     = uimenu(d.smclinew,'label','4','callback','sliceomatic slicecontourlinewidth 4');
  d.smcl5     = uimenu(d.smclinew,'label','5','callback','sliceomatic slicecontourlinewidth 5');
  d.smcl6     = uimenu(d.smclinew,'label','6','callback','sliceomatic slicecontourlinewidth 6');
  
  % Isosurface Context Menu
  d.uiciso=uicontextmenu('callback',@isocontextmenu);
  d.vistogiso = uimenu(d.uiciso,'label','Visible','callback','sliceomatic isotogglevisible');
  d.isodelete = uimenu(d.uiciso,'label','Delete','callback','sliceomatic isodelete');
  d.isoflatlight=uimenu(d.uiciso,'label','Lighting Flat','callback','sliceomatic isoflatlight','separator','on');
  d.isosmoothlight=uimenu(d.uiciso,'label','Lighting Smooth','callback','sliceomatic isosmoothlight');
  d.isocolor = uimenu(d.uiciso,'label','Change Color','callback','sliceomatic isocolor','separator','on');
  d.isoalpha=uimenu(d.uiciso,'label','Change Transparency');
  uimenu(d.isoalpha,'label','.2','callback','sliceomatic isoalpha .2');
  uimenu(d.isoalpha,'label','.5','callback','sliceomatic isoalpha .5');
  uimenu(d.isoalpha,'label','.8','callback','sliceomatic isoalpha .8');
  uimenu(d.isoalpha,'label','1','callback','sliceomatic isoalpha 1');
  d.isocap=uimenu(d.uiciso,'label','Add IsoCaps','callback','sliceomatic isocaps','separator','on');
  
  outd = d;
  
function controlmenu(fig, action)  
% Handle doing things to the CONTROLS menu
  
  d=getappdata(gcf,'sliceomatic');

  if cameratoolbar('getvisible')
    set(d.camtoolbar,'checked','on');
  else
    set(d.camtoolbar,'checked','off');
  end
  
  if exist('uitoolfactory') == 2
    if propcheck(d.toolbar,'visible','on')
      set(d.anntoolbar,'checked','on');
    else
      set(d.anntoolbar,'checked','off');
    end
  end
  
  set([d.dcalpha1 d.dcalpha8 d.dcalpha6 d.dcalpha5 d.dcalpha6 d.dcalpha2 d.dcalpha0...
       d.dclabels d.dcvis ],...
      'checked','off');

  switch get(d.pxx,'facealpha')
   case 1,  set(d.dcalpha1,'checked','on');
   case .8, set(d.dcalpha8,'checked','on');
   case .6, set(d.dcalpha6,'checked','on');
   case .5, set(d.dcalpha5,'checked','on');
   case .4, set(d.dcalpha4,'checked','on');
   case .2, set(d.dcalpha2,'checked','on');
   case 0,  set(d.dcalpha0,'checked','on');
  end

  if d.animincrement == 0
    set(d.dcanimstep,'checked','off');
  else
    set(d.dcanimstep,'checked','on');
  end
  
  if ~isempty(get(d.axx,'xticklabel'))
    set(d.dclabels,'checked','on');
  end
  
  if strcmp(get(d.axx,'visible'),'on')
    set(d.dcvis,'checked','on');
  end
  
  if 0
    xt = get(get(d.axx,'title'),'string');
    switch xt
     case 'X Slice Controller'
      set(d.dcslice,'checked','on');
    end
    
    xt = get(get(d.axiso,'title'),'string');
    switch xt
     case 'Iso Surface Controller'
      set(d.dciso,'checked','on');
    end
  end
  
function defaultmenu(fig, action)
% Handle toggling bits on the slice defaults menu
  
  d=getappdata(gcf,'sliceomatic');
  
  set([d.dfacet d.dflat d.dinterp d.dtex d.dtnone d.dtflat d.dtinterp ...
       d.dttex d.dcflat d.dcinterp d.dcblack d.dcwhite d.dcnone ...
       d.dlflat d.dlsmooth ...
       d.smcl1 d.smcl2 d.smcl3 d.smcl4 d.smcl5 d.smcl6 ], 'checked','off');
  switch d.defcolor
   case 'faceted'
    set(d.dfacet,'checked','on');
   case 'flat'
    set(d.dflat,'checked','on');
   case 'interp'
    set(d.dinterp,'checked','on');
   case 'texture'
    set(d.dtex,'checked','on');
   case 'none'
    set(d.dcnone,'checked','on');
  end
  switch d.defalpha
   case 'none'
    set(d.dtnone,'checked','on');
   case 'flat'
    set(d.dtflat,'checked','on');
   case 'interp'
    set(d.dtinterp,'checked','on');
   case 'texture'
    set(d.dttex,'checked','on');
  end
  switch d.deflight
   case 'flat'
    set(d.dlflat,'checked','on');
   case 'smooth'
    set(d.dlsmooth,'checked','on');
  end
  switch d.defcontourcolor
   case 'flat'
    set(d.dcflat,'checked','on');
   case 'interp'
    set(d.dcinterp,'checked','on');
   case 'black'
    set(d.dcblack,'checked','on');
   case 'white'
    set(d.dcwhite,'checked','on');
  end
  %set(d.dcsmooth,'checked',d.defcontoursmooth);
  switch d.defcontourlinewidth
   case 1, set(d.dcl1,'checked','on');
   case 2, set(d.dcl2,'checked','on');
   case 3, set(d.dcl3,'checked','on');
   case 4, set(d.dcl4,'checked','on');
   case 5, set(d.dcl5,'checked','on');
   case 6, set(d.dcl6,'checked','on');
  end

function slicecontextmenu(fig,action)
% Context menu state for slices

  d=getappdata(gcf,'sliceomatic');

  [a s]=getarrowslice;
  set([d.smfacet d.smflat d.sminterp d.smtex d.smtnone d.smtp5 ...
       d.smtflat d.smtinterp d.smttex d.smnone d.smcsmooth
      ],'checked','off');
  set(d.vistog,'checked',get(s,'visible'));
  
  if propcheck(s,'edgec',[0 0 0])
    set(d.smfacet,'checked','on');
  elseif propcheck(s,'facec','flat')
    set(d.smflat,'checked','on');
  end
  if propcheck(s,'facec','interp')
    set(d.sminterp,'checked','on');
  end
  if propcheck(s,'facec','texturemap')
    set(d.smtex,'checked','on');
  end
  if propcheck(s,'facec','none')
    set(d.smnone,'checked','on');
  end
  if propcheck(s,'facea',1)
    set(d.smtnone,'checked','on');
  end
  if propcheck(s,'facea',.5)
    set(d.smtp5,'checked','on');
  end
  if propcheck(s,'facea','flat')
    set(d.smtflat,'checked','on');
  end
  if propcheck(s,'facea','interp')
    set(d.smtinterp,'checked','on');
  end
  if propcheck(s,'facea','texturemap')
    set(d.smttex,'checked','on');
  end
  cm = [d.smcflat d.smcinterp d.smcblack d.smcwhite d.smccolor ...
       d.smcl1 d.smcl2 d.smcl3 d.smcl4 d.smcl5 d.smcl6 ];
  set(cm,'checked','off');
  if isempty(getappdata(s,'contour'))
    set(d.smcontour,'enable','on');
    set(d.smcsetauto,'enable','off');
    set(d.smcsetav,'enable','off');
    set(d.smclevels,'enable','off');
    set(d.smrcontour,'enable','off');
    set(d.smcsmooth,'enable','off');
    set(cm,'enable','off');
  else
    set(d.smcontour,'enable','off')
    set(d.smcsetauto,'enable','on');
    set(d.smcsetav,'enable','on');
    set(d.smclevels,'enable','on');
    set(d.smrcontour,'enable','on')
    set(d.smcsmooth,'enable','on');
    set(cm,'enable','on')
    c = getappdata(s,'contour');
    if propcheck(c,'linesmoothing','on')
      set(d.smcsmooth,'checked','on');
    end
    ec = get(c,'edgecolor');
    if isa(ec,'char')
      switch ec
       case 'flat'
        set(d.smcflat,'checked','on');
       case 'interp'
        set(d.smcinterp,'checked','on');
      end
    else
      if ec == [ 1 1 1 ]
        set(d.smcwhite,'checked','on');
      elseif ec == [ 0 0 0 ]
        set(d.smcblack,'checked','on');
      else
        set(d.smccolor,'checked','on');
      end
    end
    clw = get(c,'linewidth');
    switch clw
     case 1, set(d.smcl1,'checked','on');
     case 2, set(d.smcl2,'checked','on');
     case 3, set(d.smcl3,'checked','on');
     case 4, set(d.smcl4,'checked','on');
     case 5, set(d.smcl5,'checked','on');
     case 6, set(d.smcl6,'checked','on');
    end
  end
  
function isocontextmenu(fig,action)
% Context menu state for isosurfaces

  d=getappdata(gcf,'sliceomatic');
  
  [a s]=getarrowslice;
  if propcheck(s,'facelighting','flat')
    set(d.isoflatlight,'checked','on');
    set(d.isosmoothlight,'checked','off');
  else
    set(d.isoflatlight,'checked','off');
    set(d.isosmoothlight,'checked','on');
  end
  set(d.vistogiso,'checked',get(s,'visible'));
  if ~isempty(getappdata(s,'isosurfacecap'))
    set(d.isocap,'checked','on');
  else
    set(d.isocap,'checked','off');
  end


%%% start patch 3of3 RAB 2/18/05 %%%
%----------------------------------------------------------------------
function SavePrefs(obj,event)

%appdata structure knows everything about the implementation
d = getappdata(gcf,'sliceomatic');

%extract only preferences that need to be sticky
prefs.anntoolbar_Checked = get(d.toolbar,'Visible');
prefs.defcolor = d.defcolor;
prefs.defalpha = d.defalpha;
prefs.deflight = d.deflight;
prefs.defcontourcolor = d.defcontourcolor;
prefs.defcontourlinewidth = d.defcontourlinewidth;
prefs.defcontoursmooth = d.defcontoursmooth;

prefs.camtoolbar_checked = cameratoolbar('getvisible');
prefs.ticklabels = get(d.axx,'xticklabelmode');
prefs.animincrement = d.animincrement;
prefs.controlalpha = get(d.pxx,'facealpha');

%store mini structure (locally where Slice-O-Matic installed)
fileName = UserStickyPrefsFileName;
save(fileName,'prefs')
disp([ 'Saved: ' fileName])


%----------------------------------------------------------------------
function dOut = OverrideStickyUserPreferences(d)

%characteristic prefs file (stored locally where Slice-O-Matic installed)
fileName = UserStickyPrefsFileName;

%override particular field values (if file exists)
if exist(fileName,'file')
  load(fileName)
  set(d.toolbar,'visible',prefs.anntoolbar_Checked)
  if prefs.camtoolbar_checked
    cameratoolbar('show');
  else
    cameratoolbar('hide');
  end
  d.defcolor = prefs.defcolor;
  d.defalpha = prefs.defalpha;
  d.deflight = prefs.deflight;
  d.defcontourcolor = prefs.defcontourcolor;
  d.defcontourlinewidth = prefs.defcontourlinewidth;
  d.defcontoursmooth = prefs.defcontoursmooth;
  
  if strcmp('auto',prefs.ticklabels)
    set([d.axx d.axiso],'xticklabelmode','auto');
    set([d.axy d.axz],'yticklabelmode','auto');
  else
    set([d.axx d.axiso],'xticklabel',[]);
    set([d.axy d.axz],'yticklabel',[]);
  end
  d.animincrement = prefs.animincrement;
  set([d.pxx d.pxy d.pxz] , 'facealpha',prefs.controlalpha);
  iso = findobj(d.axiso,'type','image');
  set(iso,'alphadata',prefs.controlalpha);

  disp('Sticky preferences loaded.')
end

%return modified structure
dOut = d;


%----------------------------------------------------------------------
function fileName = UserStickyPrefsFileName
localPath = fileparts(which(mfilename));
fileName = fullfile(localPath,'Sliceomatic.Prefs.mat');

%%% end patch 3of3 RAB 2/18/05 %%%

function outd = figtoolbar(d)
% Set up the toolbar for Sliceomatic within structure D

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  set(gcf,'toolbar','none');
  
  if exist('uitoolfactory') == 2
    
    % Create a toolbar with just the elements useful on sliceomatic
    % on it.

    d.toolbar = uitoolbar('parent',gcf);
    uitoolfactory(d.toolbar, 'Annotation.InsertRectangle');
    uitoolfactory(d.toolbar, 'Annotation.InsertEllipse');
    uitoolfactory(d.toolbar, 'Annotation.InsertTextbox');
    uitoolfactory(d.toolbar, 'Annotation.InsertArrow');
    uitoolfactory(d.toolbar, 'Annotation.InsertLine');
    uitoolfactory(d.toolbar, 'Exploration.ZoomIn');
    uitoolfactory(d.toolbar, 'Exploration.ZoomOut');
    uitoolfactory(d.toolbar, 'Exploration.Pan');
    uitoolfactory(d.toolbar, 'Exploration.Rotate');
    
    cameratoolbar('show');
    cameratoolbar('togglescenelight');
    
  else

    % We are in R13 or earlier
    try
      cameratoolbar('show');
      cameratoolbar('togglescenelight');
      %cameratoolbar('setmode','orbit');
    catch
      disp('Could not display the camera toolbar.');
    end
    
  end
  
  outd = d;
  
  function [a, s]=getarrowslice
% Return the Arrow and Slice based on the GCO

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  if isempty(getappdata(gco,'controlarrow')) && ...
        isempty(getappdata(gco,'isosurface'))
    a = gco;
    s = getappdata(a,'arrowslice');
    if isempty(s)
      s=getappdata(a,'arrowiso');
    end
  else
    s = gco;
    if ~isempty(getappdata(s,'isosurface'))
      s=getappdata(s,'isosurface');
    end
    a = getappdata(s,'controlarrow');
  end

  function isocontrols(fig, onoff)
% Set up FIG to have an ISO surface controller on the bottom.
% ONOFF indicates if the controller is being turned ON or OFF

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

% Check variables
  error(nargchk(2,2,nargin))

  d = getappdata(fig, 'sliceomatic');
  
  if onoff

    lim=[min(min(min(d.data))) max(max(max(d.data)))];
      
    set(d.axiso,'handlevisibility','on');
    set(fig,'currentaxes',d.axiso);
    set(d.axiso, 'xlim',lim,...
                 'ylim',[1 5],...
                 'clim',lim);
    image('parent',d.axiso,'cdata',1:64,'cdatamapping','direct',...
          'xdata',lim,'ydata',[0 5],...
          'alphadata',.6, ...
          'hittest','off');
    activelabel('title','Iso Surface Controller');
    set(d.axiso,'handlevisibility','off');
    
  else
    % Turn off the controller
    
    delete(findobj(d.axis,'type','image'));

  end
  
  function popset(handle,prop)
% POPSET - pop values for a property from a value stack.
%
% POPSET(HANDLE, PROP) will restore a prevously HGPUSHED property value.
%

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

%  nargchk(2,2,'wrong number of arguments.');

  proplist=fieldnames(get(handle(1)));
  prop=proplist{strcmpi(prop,proplist)};

  appstr = [prop '_hgstack'];

  for k=1:prod(size(handle))
    
    olds = getappdata(handle(k),appstr);

    if length(olds) <= 1
      error(['Nothing left to pop for property ' prop '.']);
    end

    set(handle(k),prop,olds{1});
    setappdata(handle(k),appstr,olds{2:end});

  end
  
function tf=propcheck(obj, prop, value)
% Check to see if PROP for OBJ has VALUE

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc
  
  try
    v=get(obj,prop);
  catch
    tf = 0;
    return
  end

  if isa(v,class(value))
    if isa(v,'char')
      tf=strcmp(v,value);
    else
      if v==value
        tf=1;
      else
        tf=0;
      end
    end
  else
    tf=0;
  end
  
function pushset(handle,prop,value)
% PUSHSET - push new properties onto a value stack.
%
% PUSHSET(HANDLE, PROP, VALUE) will take the old value of PROP for
% HANDLE, and save it on a stack associated with HANDLE.  It then
% assigns VALUE as the new value.
%

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  proplist=fieldnames(get(handle(1)));
  prop=proplist{strcmpi(prop,proplist)};
  appstr = [prop '_hgstack'];

  for k=1:prod(size(handle))

    oldv = get(handle(k),prop);
    olds = getappdata(handle(k),appstr);
    set(handle(k),prop,value);
    setappdata(handle(k),appstr,{ oldv olds });
    
  end

function levels = selectcontourlevels(data, min, max)
% L = SELECTCONTOURLEVELS(DATA) - select contour levels for DATA.
%              Returns a vector of levels.
% L = SELECTCONTOURLEVELS(DATA, min, max) - select contour levels for DATA.
%              MINinimum and MAXimum possible contour levels
%              specified.
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  dlg_title = 'Sliceomatic contour select';
  
  if nargin == 1
    min = min(min(data));
    max = max(max(data));
  end

  levels = [];
  
  prompt = ['Enter contour values over '...
            num2str(min) ' to '...
            num2str(max) ':'];
  answer = inputdlg(prompt,dlg_title);
  if isempty(answer)~=1 & isequal(answer{1},'')~=1
    levels = str2num(answer{1});
    if isempty(levels)==1
      warndlg('Contour values must be numeric: data not accepted.');
    end
  else
    warndlg('Empty value: data not accepted.');
  end
  
% End

function setpointer(fig, ptr)
% Set the pointer on the current figure to PTR
% has several specialized SOM (SliceOMatic) pointers
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  switch ptr
   case 'SOM left'
    pd = [ nan nan nan nan 1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan 1   2   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan 1   2   1   1   1   1   1   1   1   1   1   1   1   1
           nan 1   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           nan 1   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           1   2   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           1   2   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           nan 1   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           nan 1   2   2   2   2   2   2   2   2   2   2   2   2   2   1
           nan nan 1   2   1   1   1   1   1   1   1   1   1   1   1   1
           nan nan 1   2   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan nan nan nan nan nan
           nan nan nan nan 1   nan nan nan nan nan nan nan nan nan nan nan ];
    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 8 1 ] , ...
            'pointer','custom');
   case 'SOM right'
    pd = [ nan nan nan nan nan nan nan nan nan nan nan 1  nan nan nan nan
           nan nan nan nan nan nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  2   1   nan nan 
           1   1   1   1   1   1   1   1   1   1   1   1  2   1   nan nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   2   1   
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   2   1   
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   1   1   1   1   1   1   1   1   1   1   1  2   1   nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  2   1   nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan nan nan nan nan nan nan nan 1  nan nan nan nan ];
    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 8 16 ] , ...
            'pointer','custom');
   case 'SOM bottom'
    pd = [ nan nan nan nan 1   1   1   1 1 1   1   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           1   1   1   1   1   2   2   2 2 2   2   1   1   1   1   1  
           nan 1   1   2   2   2   2   2 2 2   2   2   2   1   1   nan
           nan nan nan 1   1   2   2   2 2 2   2   1   1   nan nan nan
           nan nan nan nan nan 1   1   2 2 1   1   nan nan nan nan nan
           nan nan nan nan nan nan nan 1 1 nan nan nan nan nan nan nan ];

    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 16 8 ] , ...
            'pointer','custom');
   
   case 'SOM top'
    pd = [ nan nan nan nan nan nan nan 1 1 nan nan nan nan nan nan nan
           nan nan nan nan nan 1   1   2 2 1   1   nan nan nan nan nan
           nan nan nan 1   1   2   2   2 2 2   2   1   1   nan nan nan
           nan 1   1   2   2   2   2   2 2 2   2   2   2   1   1   nan
           1   1   1   1   1   2   2   2 2 2   2   1   1   1   1   1  
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   1   1   1 1 1   1   1   nan nan nan nan ];

    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 1 8 ] , ...
            'pointer','custom');
   
   case 'SOM leftright'
    pd = [ nan nan nan nan 1   nan nan nan nan nan nan 1  nan nan nan nan
           nan nan nan 1   1   nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan 1   1   nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan 1   2   1   nan nan nan nan nan nan 1  2   1   nan nan 
           nan nan 1   2   1   1   1   1   1   1   1   1  2   1   nan nan 
           nan 1   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           nan 1   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   2   1   
           1   2   2   2   2   2   2   2   2   2   2   2  2   2   2   1   
           nan 1   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           nan 1   2   2   2   2   2   2   2   2   2   2  2   2   1   nan 
           nan nan 1   2   1   1   1   1   1   1   1   1  2   1   nan nan 
           nan nan 1   2   1   nan nan nan nan nan nan 1  2   1   nan nan 
           nan nan nan 1   1   nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan 1   1   nan nan nan nan nan nan 1  1   nan nan nan 
           nan nan nan nan 1   nan nan nan nan nan nan 1  nan nan nan nan ];
    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 8 8 ] , ...
            'pointer','custom');
   
   case 'SOM topbottom'
    pd = [ nan nan nan nan nan nan nan 1 1 nan nan nan nan nan nan nan
           nan nan nan nan nan 1   1   2 2 1   1   nan nan nan nan nan
           nan nan nan 1   1   2   2   2 2 2   2   1   1   nan nan nan
           nan 1   1   2   2   2   2   2 2 2   2   2   2   1   1   nan
           1   1   1   1   1   2   2   2 2 2   2   1   1   1   1   1  
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           nan nan nan nan 1   2   2   2 2 2   2   1   nan nan nan nan
           1   1   1   1   1   2   2   2 2 2   2   1   1   1   1   1  
           nan 1   1   2   2   2   2   2 2 2   2   2   2   1   1   nan
           nan nan nan 1   1   2   2   2 2 2   2   1   1   nan nan nan
           nan nan nan nan nan 1   1   2 2 1   1   nan nan nan nan nan
           nan nan nan nan nan nan nan 1 1 nan nan nan nan nan nan nan ];


    set(fig,'pointershapecdata', pd,...
            'pointershapehotspot', [ 8 8 ] , ...
            'pointer','custom');
    
   otherwise
    % Set it to the string passed in
     set(fig,'pointer', ptr);
  end
  
function setvolumerange
% Query for a new volume range based on the sliceomatic gui
% which should be GCF
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  %d=getappdata(fig,'sliceomatic');
  
  p=get(fig,'position');
  np=[p(1)+20 p(2)+30 400 200];
  figure('position',np);
  
  uicontrol('units','norm','style','text','string','X Range',...
            'position',[0 .6 .3 .3]);
  uicontrol('units','norm','style','text','string','Y Range',...
            'position',[0 .3 .3 .3]);
  uicontrol('units','norm','style','text','string','Z Range',...
            'position',[0 0  .3 .3]);
  
function showarrowtip (arrow)
% Display a tip for ARROW.
% Depends on tipdata being set on the handle to ARROW.
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  d=getappdata(gcf,'sliceomatic');
  
  if isempty(arrow)
    tipdata = [];
  else
    
    ctrlarrow = getappdata(arrow,'controlarrow');
    
    if ~isempty(ctrlarrow)
      % In this case, the slice or isosurface passed in
      % most likely from the motion callback.  Lets redirect
      % our input argument and show the tip anyway.
      
      % See the "arrow" private function for the setting of this value;
      arrow = ctrlarrow(2);
    end
    
    tipdata = getappdata(arrow,'tipdata');
  end
  
  if ~isempty(tipdata)
  
    set(d.tip,'parent',tipdata.parentaxes, ...
              'string',sprintf('Value: %1.3f',tipdata.value),...
              'units','data', ...
              'position', tipdata.position, ...
              'verticalalignment', tipdata.verticalalign,...
              'horizontalalignment', tipdata.horizontalalign);
    set(d.tip,'units','pixels');
    set(d.tip,'visible','on');
    
  else
    
    set(d.tip,'visible','off');
    
  end

function slicecontrols(fig,onoff,xmesh,ymesh,zmesh,xdir,ydir,zdir)
% Convert figure to contain controls for manipulating slices.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

% Check variables
  error(nargchk(2,8,nargin))

  d = getappdata(fig, 'sliceomatic');
  
  if onoff
    
    if nargin ~= 8
      
      % If there is no user supplied mesh, make one up.
      xmesh(1) = 1;
      xmesh(2) = size(d.data,2);
      ymesh(1) = 1;
      ymesh(2) = size(d.data,1);
      zmesh(1) = 1;
      zmesh(2) = size(d.data,3);
      
      xdir = 'normal';
      ydir = 'normal';
      zdir = 'normal';
    end

    set(0,'currentfigure',fig);
    set([d.axx d.axy d.axz] ,'handlevisibility','on');
    
    set(fig,'currentaxes',d.axx);
    set(d.axx, 'xlim',[xmesh(1) xmesh(end)],...
               'ylim',[1 5]);
    set(d.pxx, 'vertices',[ xmesh(1) xmesh(1) -1; xmesh(end) xmesh(1) -1; xmesh(end) 5 -1; xmesh(1) 5 -1],...
               'faces',[ 1 2 3 ; 1 3 4]);
    
    activelabel('title', 'X Slice Controller');
    
    set(fig,'currentaxes',d.axy);
    set(d.axy, 'xlim',[1 5],...
               'ylim',[ymesh(1) ymesh(end)]);
    set(d.pxy, 'vertices',[ ymesh(1) ymesh(1) -1; ymesh(1) ymesh(end) -1; 5 ymesh(end) -1; 5 ymesh(1) -1],...
               'faces',[ 1 2 3 ; 1 3 4]);
    activelabel('title', 'Y Slice');

    set(fig,'currentaxes',d.axz);
    set(d.axz, 'xlim',[1 5],...
               'ylim',[zmesh(1) zmesh(end)]);
    set(d.pxz, 'vertices',[ zmesh(1) zmesh(1) -1; zmesh(1) zmesh(end) -1; 5 zmesh(end) -1; 5 zmesh(1) -1],...
               'faces',[ 1 2 3 ; 1 3 4]);
    activelabel('title', 'Z Slice');
    
    set([d.axx d.axy d.axz] ,'handlevisibility','off');
    
    set(d.axx,'xdir',xdir);
    set(d.axy,'ydir',ydir);
    set(d.axz,'zdir',zdir);

  else
    
    % Disable these controls.  Perhaps hide all slices?
    
  end

function appdata=sliceomaticfigure(d,xmesh,ymesh,zmesh)
% FIG=SLICEOMATICFIGURE (D,XMESH,YMESH,ZMESH) - 
% Create the figure window to be used by the sliceomatic GUI.
% D is the app data to attach to the figure
% The [XYZ]MESH arguments specify a mesh that the data D falls into.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

% Check variables
error(nargchk(1,4,nargin))

% Init sliceomatic
  figure('name','Sliceomatic','toolbar','none');
  lim=[min(min(min(d.data))) max(max(max(d.data)))];
  if nargin==4
      % Reorder vectors: make them horizontal (prepare to flipdim)
      if size(xmesh,1)>size(xmesh,2) 
          xmesh=xmesh';
      end
      if size(ymesh,1)>size(ymesh,2) 
          ymesh=ymesh';
      end
      if size(zmesh,1)>size(zmesh,2) 
          zmesh=zmesh';
      end
      % Set axis orientation
      xdir='normal';
      ydir='normal';
      zdir='normal';
      if issorted(xmesh)~=1
          xmesh=flipdim(xmesh,2);
          xdir='reverse';
      end
      if issorted(ymesh)~=1
          ymesh=flipdim(ymesh,2);
          ydir='reverse';
      end
      % This should not be the case for medical images
      if issorted(zmesh)~=1
          zmesh=flipdim(zmesh,2);
          zdir='reverse';
      end
      % Update data structure
      
      d.axmain = axes('units','normal','pos',[.2  .2 .6 .6],'box','on',...
          'ylim',[ymesh(1) ymesh(end)],...
          'xlim',[xmesh(1) xmesh(end)],...
          'zlim',[zmesh(1) zmesh(end)],...
          'clim',lim,...
          'alim',lim);
      % Set axes direction
      set(gca,'XDir',xdir,'YDir',ydir,'ZDir',zdir);
    else
      d.axmain = axes('units','normal','pos',[.2  .2 .6 .6],'box','on',...
          'ylim',[1 size(d.data,1)],...
          'xlim',[1 size(d.data,2)],...
          'zlim',[1 size(d.data,3)],...
          'clim',lim,...
          'alim',lim);
  end
  
  activelabel('xlabel', 'X');
  activelabel('ylabel', 'Y');
  activelabel('zlabel', 'Z');
  %activelabel('title', 'Data');
  daspect([1 1 1]);
  view(3);
  axis tight vis3d;
  hold on;
  grid on;
  
  % Set up the four controller axes.
  d.axx    = axes('units','normal','pos',[.2  .81 .6 .1],'box','on',...
                  'ytick',[],'xgrid','on','xaxislocation','top',...
                  'zlim',[-2 1 ],...
                  'layer','top',...
                  'color','none');
  d.pxx    = patch('facecolor',[1 1 1],...
                   'facealpha',.6,...
                   'edgecolor','none',...
                   'hittest','off');
  setappdata(d.axx,'motionpointer','SOM bottom');
  d.axy    = axes('units','normal','pos',[.05 .05 .1 .75],'box','on',...
                  'xtick',[],'ygrid','on',...
                  'zlim',[-2 1 ],...
                  'layer','top',...
                  'color','none');
  d.pxy    = patch('facecolor',[1 1 1],...
                   'facealpha',.6,...
                   'edgecolor','none',...
                   'hittest','off');
  setappdata(d.axy,'motionpointer','SOM right');
  d.axz    = axes('units','normal','pos',[.85 .05 .1 .75],'box','on',...
                  'xtick',[],'ygrid','on','yaxislocation','right',...
                  'zlim',[-2 1 ],...
                  'layer','top',...
                  'color','none');
  d.pxz    = patch('facecolor',[1 1 1],...
                   'facealpha',.6,...
                   'edgecolor','none',...
                   'hittest','off');
  setappdata(d.axz,'motionpointer','SOM left');
  d.axiso  = axes('units','normal','pos',[.2 .05 .6 .1],'box','on',...
                  'ytick',[],'xgrid','off','ygrid','off',...
                  'xaxislocation','bottom',...
                  'zlim',[-1 1],...
                  'color','none',...
                  'layer','top');
  setappdata(d.axiso,'motionpointer','SOM top');
  set([d.axx d.axy d.axz d.axiso],'handlevisibility','off');

  setappdata(gcf,'sliceomatic',d);
  
  % Set up the default sliceomatic controllers
  if nargin == 4 
      slicecontrols(gcf,1,xmesh,ymesh,zmesh,xdir,ydir,zdir);
  else
      slicecontrols(gcf,1);
  end
      
  isocontrols(gcf,1);

  % Button Down Functions
  set(d.axx,'buttondownfcn','sliceomatic Xnew');
  set(d.axy,'buttondownfcn','sliceomatic Ynew');
  set(d.axz,'buttondownfcn','sliceomatic Znew');
  set(d.axiso,'buttondownfcn','sliceomatic ISO');

  % Set up our motion function before cameratoolbar is active.
  d.motionmetaslice = [];
  set(gcf,'windowbuttonmotionfcn',@sliceomaticmotion);
  
  % Try setting up the camera toolbar
  d=figtoolbar(d);
  
  d = figmenus(d);
  
  % Color and alph maps
  uicontrol('style','text','string','ColorMap',...
            'units','normal','pos',[0 .9 .19 .1]);
  uicontrol('style','popup','string',...
            {'jet','hsv','cool','hot','pink','bone','copper','flag','prism','rand','custom'},...
            'callback','sliceomatic colormap',...
            'units','normal','pos',[0 .85 .19 .1]);

  uicontrol('style','text','string','AlphaMap',...
            'units','normal','pos',[.81 .9 .19 .1]);
  uicontrol('style','popup','string',{'rampup','rampdown','vup','vdown','rand'},...
            'callback','sliceomatic alphamap',...
            'units','normal','pos',[.81 .85 .19 .1]);

  % Data tip thingydoo
  d.tip = text('visible','off','fontname','helvetica','fontsize',10,'color','black');
  try
    % Try R13 new feature
    set(d.tip,'backgroundcolor',[1 1 .8],'edgecolor',[.5 .5 .5],'margin',5);
  end
  
  appdata = d;

  set(gcf,'nextplot','new');
  set(gca,'nextplot','new');

function sliceomaticmotion(fig,action)
% Handle generic motion events for the figure window.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  obj = hittest(fig);

  % Some objects get a special pointer when the mouse waves
  % over them.  Get it from appdata.
  if ~isempty(obj)
    t = getappdata(obj,'motionpointer');
    cc = get(fig,'pointer');
  
    if t
      newc = t;
    else
      newc = get(0,'defaultfigurepointer');
    end
  
    if isa(newc,'char') && isa(cc,'char') && ~strcmp(newc,cc)
      setpointer(fig, newc);
    end
  end
  
  d = getappdata(fig,'sliceomatic');

  % The motion meta slice is a line that is managed here.
  % There is only one.
  if isempty(d.motionmetaslice)
    d.motionmetaslice = line('parent',d.axmain,...
                             'vis','off',...
                             'linestyle','--',...
                             'marker','none',...
                             'linewidth',2,...
                             'erasemode','xor','clipping','off');
    setappdata(fig,'sliceomatic',d);
  end

  showarrowtip(obj);
  
  if isempty(obj) || (obj ~= d.axx && obj ~= d.axy && obj ~= d.axz)
    set(d.motionmetaslice,'visible','off');
    
    return
  end
  
  % OBJ can only be an Axes because of the previous IF statement.  
  aa = obj;
  apos=get(aa,'currentpoint');

  xl = d.xlim;
  yl = d.ylim;
  zl = d.zlim;
  
  if aa==d.axx || aa==d.axiso
    if aa==d.axiso
      % eh?
    else
      xdata = [ apos(1,1) apos(1,1) apos(1,1) apos(1,1) apos(1,1) ];
      ydata = [ yl(1) yl(2) yl(2) yl(1) yl(1) ];
      zdata = [ zl(2) zl(2) zl(1) zl(1) zl(2) ];
    end
  else
    % We are moving a Y or Z slice
    if aa==d.axy
      ydata = [ apos(1,2) apos(1,2) apos(1,2) apos(1,2) apos(1,2) ];
      xdata = [ xl(1) xl(2) xl(2) xl(1) xl(1) ];
      zdata = [ zl(2) zl(2) zl(1) zl(1) zl(2) ];
    else
      zdata = [ apos(1,2) apos(1,2) apos(1,2) apos(1,2) apos(1,2) ];
      ydata = [ yl(1) yl(2) yl(2) yl(1) yl(1) ];
      xdata = [ xl(2) xl(2) xl(1) xl(1) xl(2) ];
    end
  end

  set(d.motionmetaslice,'visible','on',...
                    'xdata',xdata,'ydata',ydata,'zdata',zdata);

function appdata = sliceomaticsetdata(d,xmesh,ymesh,zmesh)
% SLICEOMATICSETDATA(rawdata) - Create the data used for
% sliceomatic in the appdata D.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

% Check variables
  error(nargchk(1,4,nargin))

  % Simplify the isonormals
  %as smooth3 is taking the most time for large volumes
  if size(d.data,1)*size(d.data,2)*size(d.data,3)>10^6
      disp('SKIPPING: Smoothing for IsoNormals...');
      disp('Do not use ISO surface controller for this large volume...');
      d.smooth=d.data;
  else
      disp('Smoothing for IsoNormals...');
      d.smooth=smooth3(d.data);  % ,'box',5);
  end
  d.reducenumbers=[floor(size(d.data,2)/20)...
                   floor(size(d.data,1)/20)...
                   floor(size(d.data,3)/20) ];
  d.reducenumbers(d.reducenumbers==0)=1;


  if nargin == 4
    % Reorder vectors: make them horizontal (prepare to flipdim)
    if size(xmesh,1)>size(xmesh,2) 
      xmesh=xmesh';
    end
    if size(ymesh,1)>size(ymesh,2) 
      ymesh=ymesh';
    end
    if size(zmesh,1)>size(zmesh,2) 
      zmesh=zmesh';
    end
    % Set axis orientation
    xdir='normal';
    ydir='normal';
    zdir='normal';
    if issorted(xmesh)~=1
      xmesh=flipdim(xmesh,2);
      xdir='reverse';
      d.xlim = [xmesh(1) xmesh(end)];
      xmesh=flipdim(xmesh,2);
    else
      d.xlim = [xmesh(1) xmesh(end)];
    end
    if issorted(ymesh)~=1
      ymesh=flipdim(ymesh,2);
      ydir='reverse';
      d.ylim = [ymesh(1) ymesh(end)];
      ymesh=flipdim(ymesh,2);
    else
      d.ylim = [ymesh(1) ymesh(end)];
    end
    % This should not be the case for medical images
    if issorted(zmesh)~=1
      zmesh=flipdim(zmesh,2);
      zdir='reverse';
      d.zlim = [zmesh(1) zmesh(end)];
      zmesh=flipdim(zmesh,2);
    else
      d.zlim = [zmesh(1) zmesh(end)];
    end
    
    % Vol vis suite takes numbers in X/Y form.
    ly = 1:d.reducenumbers(1):size(d.data,2);
    lx = 1:d.reducenumbers(2):size(d.data,1);
    lz = 1:d.reducenumbers(3):size(d.data,3);
    
    for i = 1:length(ly)
      ly(i) = xmesh(ly(i));
    end
    for i = 1:length(lx)
      lx(i) = ymesh(lx(i));
    end
    for i = 1:length(lz)
      lz(i) = zmesh(lz(i));
    end
    
    d.reducelims={ ly lx lz };
    disp('Generating reduction volume...');
    d.reduce= reducevolume(d.data,d.reducenumbers);
    d.reducesmooth=smooth3(d.reduce,'box',5);
    % Set axis
    %d.xlim = [xmesh(1) xmesh(end)];
    %d.ylim = [ymesh(1) ymesh(end)];
    %d.zlim = [zmesh(1) zmesh(end)];
    d.xmesh = xmesh;
    d.ymesh = ymesh;
    d.zmesh = zmesh;
    d.xdir = xdir;
    d.ydir = ydir;
    d.zdir = zdir;
    
  else
    % Vol vis suite takes numbers in X/Y form.
    ly = 1:d.reducenumbers(1):size(d.data,2);
    lx = 1:d.reducenumbers(2):size(d.data,1);
    lz = 1:d.reducenumbers(3):size(d.data,3);
    
    d.reducelims={ ly lx lz };
    disp('Generating reduction volume...');
    d.reduce= reducevolume(d.data,d.reducenumbers);
    d.reducesmooth=smooth3(d.reduce,'box',5);

    d.xlim = [1 size(d.data,2)];
    d.ylim = [1 size(d.data,1)];
    d.zlim = [1 size(d.data,3)];
    d.xmesh = nan;
    d.ymesh = nan;
    d.zmesh = nan;
    d.xdir = 'normal';
    d.ydir = 'normal';
    d.zdir = 'normal';
  end
  
  appdata = d;
  
function slowset(handle, prop, value, increment)
% SLOWSET(H, 'PROPERTY', VALUE)
%
% Like SET, except that the property is set against the original
% value over several steps such that the value morphs from the
% starting value to the end value.
%
% H can be a vector of handles, but they must all accept PROPERTY.
%
% Only one property is supported at this time.
%
% Optional fourth argument INCREMENT specifies how many steps of
% animation to use.
  
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

  global INCREMENT;
  
  if nargin == 4
    INCREMENT = increment;
    if INCREMENT==0
      INCREMENT=1;
    end
  else
    INCREMENT=10;
  end
  
%  nargchk(3,3,'wrong number of arguments.');
  
  proplist=fieldnames(get(handle(1)));
  tprop={ proplist{strncmpi(prop,proplist,length(prop))} };
  prop=tprop{1};
  
  hp = [];
  
  for i = 1:length(handle)
    hp(i).handle = handle(i);
    hp(i).start = get(hp(i).handle,prop);
    hp(i).end = value;
    if isnumeric(hp(i).end) && isnumeric(hp(i).start)
      hp(i).values = VectorCalc(hp(i));
    else
      set(hp(i).handle,prop,value);
      hp(i).values = [];
    end
  end
  
  for inc = 1:INCREMENT
    for i = 1:length(handle)
      if ~isempty(hp(i).values)
        newval = reshape(hp(i).values(inc,:,:,:),...
                         size(hp(i).start,1),...
                         size(hp(i).start,2));
        
        set(hp(i).handle,prop,newval);
      end
    end
    pause(.05)
  end
  
function values =  VectorCalc(hp)
% Do nothing but go to end value.

  global INCREMENT;
  
  s = prod(size(hp.end));
  
  values = ones(INCREMENT,size(hp.end,1), size(hp.end,2),size(hp.end,3));
  
  for c = 1:s
    newval = linspace(hp.start(c),hp.end(c),INCREMENT);
    values(:,c) = newval';
  end
  
  values = reshape(values, INCREMENT, size(hp.end,1), size(hp.end,2), ...
                   size(hp.end,3));