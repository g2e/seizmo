function [varargout]=plotmt(x,y,mt,varargin)
%PLOTMT    Plot moment tensor(s) as beach ball(s)
%
%    Usage:    plotmt(x,y,mt)
%              plotmt(x,y,mt,'param1',val1,...,'paramN',valN)
%              h=plotmt(...)
%
%    Description:
%     PLOTMT(X,Y,MT) plots the moment tensor(s) given by MT at the
%     location(s) specified with X & Y in the current plot axes as 2D beach
%     balls (view showing downward radiation pattern).  X & Y are assumed
%     to be East & North making them roughly equivalent to longitude and
%     latitude.  MT must be in Harvard convention with one of the following
%     formats:
%       - Nx6    [Mrr Mtt Mpp Mrt Mrp Mtp]
%       - 3x3xN  [Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp]
%       - scalar struct with fields .mrr .mtt .mpp .mrt .mrp .mtp
%     where N is the number of moment tensors to plot.  Note that X & Y
%     must be scalar or have N elements.
%
%     PLOTMT(X,Y,MT,'PARAM1',VAL1,...'PARAMN',VALN) sets parameters used in
%     plotting the beach balls.  Valid parameters are:
%
%        NAME         DESCRIPTION
%      ---------+----------------------------------
%      'wave'   | wave type used in calculating the beach balls (p, sh, sv)
%      'pov'    | point-of-view from up/North (where are you looking from?)
%      'roll'   | rotate drawn beach ball clockwise by this many degrees
%      'face'   | shown hemisphere of the beach balls (front or back)
%      'tcolor' | color of tension region (contains pressure axis)
%      'pcolor' | color of compression region (contains tension axis)
%      'radius' | radius of beach balls
%      'parent' | handle of axes to plot in
%      'lcolor' | outline color
%      'lwidth' | outline line width
%      'sample' | degree step size in calculating values in beach balls
%
%        NAME         INPUT      MULTI      DEFAULT          ISSUES
%      ---------+--------------+-------+---------------+----------------
%      'wave'   |   'p/sh/sv'  |   Y   | 'p'           | 
%      'pov'    |   [inc az]   |   Y   | [0 180]       | see Notes
%      'roll'   |    angle     |   Y   | 0             | see Notes
%      'face'   | 'front/back' |   Y   | 'back'        | see Notes
%      'tcolor' |   [r g b]    |   Y   | [1 1 1]       | 
%      'pcolor' |   [r g b]    |   Y   | [0 0 0]       | 
%      'radius' |   radius     |   Y   | .5            | 
%      'parent' |   axhandle   |   N   | gca           | 
%      'lcolor' |   [r g b]    |   Y   | [0 0 0]       | 
%      'lwidth' |    width     |   Y   | .5            | 
%      'sample' |   degrees    |   Y   | 5             | 
%
%     H=PLOTMT(...) returns the graphics handles of the patches and lines
%     making up the beach balls in the cell array H.  Patches are in H{:,1}
%     and outlines in H{:,2}.
%
%    Notes:
%     - The point-of-view and face options are defaulted to a view that is
%       a southish perspective (as if your feet are pointed south as you
%       are looking from above at the backside of the solution -- see the
%       'face' option to see the frontside).  Inclination is the angle from
%       up (THIS IS OPPOSITE TO RADPAT!) while azimuth is the clockwise
%       angle from North.  For example, a point-of-view of [90 0] sets the
%       view as from a northern perspective while standing upright on the
%       horizontal and [90 90] gives a similar view from the East.
%     - The roll option is useful for plotting moment tensors on maps and
%       mantle profiles where north/up are not necessarily up in your
%       figure.  MAPCMTS uses this option to orient the moment tensors
%       correctly in the map projection.  Rotates clockwise.
%
%    Examples:
%     % Plot a 10x10 grid of the first 100 cmts in the globalcmt database:
%     y=1:10; y=y(ones(10,1),:); x=y';
%     plotmt(x(:),y(:),findcmt('n',100),'pcolor',gmt_rainbow(100));
%     axis equal tight off;
%
%     % Show off the ability of PLOTMT to "rotate" a moment tensor:
%     mt=findcmt;
%     for i=0:10:360
%         for j=0:10:360
%             plotmt(0,0,mt,'pov',[i+j/36 j]);
%             drawnow;
%             pause(.01);
%         end
%     end
%
%    See also: PLOTMT3, RADPAT, MAPCMTS, FINDCMT, FINDCMTS

%     Version History:
%        Apr. 28, 2011 - initial version
%        May  15, 2011 - implimented point-of-view & roll options
%        May  31, 2011 - plotmt now uses contourf
%        June  1, 2011 - fix for isotropic source
%        June  6, 2011 - fixed bug for 3x3 checking case
%        June 14, 2011 - added handles to userdata, allow character input
%                        for colors
%        Jan. 11, 2012 - improve the notes section
%        Feb.  7, 2012 - use 3x3xN over Nx6 to allow single couples
%        Feb. 22, 2012 - minor fix for example
%        Mar. 23, 2013 - doc update, minor code fixes/refactoring
%        Mar. 25, 2013 - update for mt_change/mt_check
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 23:55 GMT

% todo:
% - tpb points
% - nodal lines

% check args
error(nargchk(3,inf,nargin));

% check optional args
if(~mod(nargin,2)) % must be 3,5,7,... inputs
    error('seizmo:plotmt:unpairedParam',...
        'PLOTMT call has an unpaired pamameter/value!');
elseif(nargin>3 && ~iscellstr(varargin(1:2:end)))
    error('seizmo:plotmt:badParam',...
        'PLOTMT parameters must be specified as strings!');
end

% check x/y/mt
if(~isnumeric(x) || ~isreal(x) || ~isnumeric(y) || ~isreal(y))
    error('seizmo:plotmt:badInput',...
        'X/Y must be real-valued numeric arrays!');
end
error(mt_check(mt));
mt=mt_change('v',mt);
n=size(mt,1);
if(~all([numel(x) numel(y)]==1 | [numel(x) numel(y)]==n))
    error('seizmo:plotmt:badInput',...
        'X/Y must be scalar or arrays with N elements!');
end

% expand x/y to Nx1
if(isscalar(x)); x(1:n)=x; end
if(isscalar(y)); y(1:n)=y; end
[x,y]=deal(x(:),y(:));

% check parameters
pv=check_plotmt_parameters(n,varargin{:});

% set handle to current axis if none give
if(isempty(pv.parent)); pv.parent=gca; end

% get hold state
held=ishold(pv.parent);

% turn off v6 contourf warning (annoying)
warnstate=warning('off','MATLAB:contourf:DeprecatedV6Argument');

% plotting each moment tensor individually
h=cell(n,2);
for i=1:n
    % deviatoric ray directions
    % note: azimuth starts at -90 as a workaround to strange contourf
    %       behavior at 0 (a line is drawn from origin to top)
    inc0=0:pv.sample(i):90; ninc=numel(inc0);
    az0=-90:pv.sample(i):270; naz=numel(az0);
    [az,inc]=meshgrid(az0,inc0); % ninc x naz
    
    % get actual ray directions
    [inc,az]=rotate_incaz(inc,az,pv,i);
    
    % calc radiation pattern in given directions
    rp=reshape(radpat(mt(i,:),inc(:),az(:),pv.wave{i}),[ninc naz]);
    
    % color matrix (for surf)
    %c=permute(pv.tcolor(i,:),[1 3 2]);
    %c=c(ones(ninc,1),ones(1,naz),:);
    %c(rp>=0)=pv.pcolor(i,1);
    %c(find(rp>=0)+ninc*naz)=pv.pcolor(i,2);
    %c(find(rp>=0)+2*ninc*naz)=pv.pcolor(i,3);
    
    % plot moment tensor (for surf)
    %h{i,1}=surf(pv.parent,...
    %    sqrt(2)*sind(inc0'/2)...
    %    *sind(pv.roll(i)+az0)*pv.radius(i)+x(i),...
    %    sqrt(2)*sind(inc0'/2)...
    %    *cosd(pv.roll(i)+az0)*pv.radius(i)+y(i),0*az,c);
    
    % plot moment tensor (for contour)
    % notes:
    % (1) v6 needed to properly draw contours (matlab regression)
    % (2) the -2 forces coloring of the entire radiation pattern
    [h{i,1},h{i,1}]=contourf(...
        'v6',sqrt(2)*sind(inc0'/2)...
        *sind(pv.roll(i)+az0)*pv.radius(i)+x(i),...
        sqrt(2)*sind(inc0'/2)...
        *cosd(pv.roll(i)+az0)*pv.radius(i)+y(i),...
        rp,[-2 0],'parent',pv.parent);
    
    % clean up patches (for contour)
    set(h{i,1},'edgecolor','none','tag','plotmt_patch');
    if(numel(h{i,1})>1)
        cdata=cell2mat(get(h{i,1},'cdata'));
    else
        cdata=get(h{i,1},'cdata');
    end
    set(h{i,1}(cdata>=0),'facecolor',pv.pcolor(i,:));
    set(h{i,1}(cdata>-2 & cdata<0),'facecolor',pv.tcolor(i,:));
    
    % turn hold on after 1st pass if not already
    if(i==1 && ~held); hold(pv.parent,'on'); end
    
    % outline
    lc=pv.lcolor(i,:); vis='on';
    if(isnan(sum(lc))); lc=[0 0 0]; vis='off'; end
    h{i,2}=plot(pv.parent,x(i)+pv.radius(i)*sind(az0),...
        y(i)+pv.radius(i)*cosd(az0),'-','color',lc,...
        'visible',vis,'linewidth',pv.lwidth(i),'tag','plotmt_outline');
    
    % save info to userdata
    tmp.mt=mt(i,:);
    tmp.wave=pv.wave{i};
    tmp.pov=pv.pov(i,:);
    tmp.roll=pv.roll(i);
    tmp.face=pv.face{i};
    tmp.tcolor=pv.tcolor(i,:);
    tmp.pcolor=pv.pcolor(i,:);
    tmp.outline=h(i,2);
    tmp.radius=pv.radius(i);
    tmp.sample=pv.sample(i);
    tmp.patch_handle=h{i,1};
    tmp.outline_handle=h{i,2};
    set(h{i,1},'userdata',tmp);
end

% turn v6 warnings back on
warning(warnstate.state,'MATLAB:contourf:DeprecatedV6Argument');

% adjust view (for surf)
%colormap(pv.parent,[pv.tcolor(1,:); pv.pcolor(1,:)]);
%view(pv.parent,0,90);
%shading(pv.parent,'interp');

% restore hold
if(~held); hold(pv.parent,'off'); end

% output if desired
if(nargout); varargout{1}=h; end

end


function [inc,az]=rotate_incaz(inc,az,pv,i)
% rotates from diavatoric units (downward hemisphere) to true position

% rad/deg
r2d=180/pi;

% flip to front if desired
switch lower(pv.face{i})
    case 'front'
        inc=180-inc;
end

% precompute
% - NOTE 180 rotation for azimuth to have
%   North perspective correspond to 0
cosinc=cosd(pv.pov(i,1));
sininc=sind(pv.pov(i,1));
cosaz=cosd(180+pv.pov(i,2));
sinaz=sind(180+pv.pov(i,2));

% translate to sph
az=(90-az)/r2d;
inc=(inc-90)/r2d;

% translate to cart
[x,y,z]=sph2cart(az,inc,1);

% rotation from viewer specific to east/north/up coordinates
% direction cosines are as follows:
% [e    [ cos(az) sin(az)*cos(inc) -sin(az)*sin(inc)    [x'
%  n  =  -sin(az) cos(az)*cos(inc) -cos(az)*sin(inc)  *  y'
%  u]           0         sin(inc)          cos(inc)]    z']
x1= cosaz*x + sinaz*cosinc*y - sinaz*sininc*z;
y1=-sinaz*x + cosaz*cosinc*y - cosaz*sininc*z;
z1=                 sininc*y +       cosinc*z;

% translate to sph
[az,inc]=cart2sph(x1,y1,z1);

% translate to incaz
az=90-az*r2d;
inc=90+inc*r2d;

end


function [pv,varargin]=check_plotmt_parameters(n,varargin)
% pv pair defaults
%pv.wave='p';       % p/sh/sv    (wavetype - rayleigh/love would be nice)
%pv.pov=[0 180];    % [inc az]   (from north=[90 0], from east=[90 90])
%pv.roll=0;         % angle      (post-rotate beach ball by this amount)
%pv.face='back';    % front/back (which side do you want to see?)
%pv.tcolor=[1 1 1]; % [r g b]    (color of tension)
%pv.pcolor=[0 0 0]; % [r g b]    (color of compression)
%pv.lcolor=[0 0 0]; % [r g b] or [nan nan nan] / [] (color of outline)
%pv.lwidth=0.5;     % width      (line width)
%pv.radius=0.5;     % radius     (radius of moment tensor)
%pv.sample=5;       % degrees    (sample interval for moment tensor values)
%pv.parent=[];      % axishandle (handle of axis to plot in)

% above is defunct (defaults are now parsed/checked/expanded like inputs)
varargin=[{'wave' 'p' 'pov' [0 180] 'roll' 0 'face' 'back' ...
    'tcolor' [1 1 1] 'pcolor' [0 0 0] 'lcolor' [0 0 0] 'lwidth' .5 ...
    'radius' .5 'sample' 5 'parent' []} varargin];

% check/parse given pv pairs
for i=1:2:numel(varargin)
    p=varargin{i};
    v=varargin{i+1};
    switch lower(p)
        case {'w' 'wave' 'type' 'wavetype'}
            if(isempty(v)); continue; end
            if(ischar(v)); v=cellstr(v); end
            if(~iscellstr(v) || ~any(numel(v)==[1 n]) ...
                    || any(~ismember(lower(v),{'p' 'sh' 'sv'})))
                error('seizmo:plotmt:badInput',...
                    'WAVE must be ''P'', ''SV'' or ''SH'' !');
            end
            if(isscalar(v)); v(1:n)=v; end
            pv.wave=v;
        case {'o' 'orient' 'view' 'pov' 'pointofview' 'point-of-view'}
            if(isempty(v)); continue; end
            if(~isnumeric(v) || ~isreal(v) || size(v,2)~=2 ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt:badInput',...
                    'POV must be given as [INCLINATION AZIMUTH]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            pv.pov=v;
        case {'ro' 'roll'}
            if(isempty(v)); continue; end
            if(~isnumeric(v) || ~isreal(v) || ndims(v)~=2 ...
                    || ~any(numel(v)==[1 n]))
                error('seizmo:plotmt:badInput',...
                    'ROLL must be a real valued number!');
            end
            if(numel(v)==1); v=v(ones(1,n),1); end
            pv.roll=v;
        case {'f' 'face'}
            if(isempty(v)); continue; end
            if(ischar(v)); v=cellstr(v); end
            if(~iscellstr(v) || ~any(numel(v)==[1 n]) ...
                    || any(~ismember(lower(v),{'front' 'back'})))
                error('seizmo:plotmt:badInput',...
                    'FACE must be ''FRONT'' or ''BACK'' !');
            end
            if(isscalar(v)); v(1:n)=v; end
            pv.face=v;
        case {'t' 'tc' 'tcolor'}
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt:badInput',...
                    'TCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.tcolor=v;
        case {'c' 'pc' 'pcolor' 'color'}
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt:badInput',...
                    'PCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.pcolor=v;
        case {'lc' 'lcolor' 'outline'}
            if(isempty(v)); pv.lcolor=nan(n,3); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt:badInput',...
                    'LCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.lcolor=v;
        case {'lw' 'lwidth' 'linewidth'}
            if(isempty(v)); continue; end
            if(~isnumeric(v) || ~isreal(v) || ndims(v)~=2 ...
                    || ~any(numel(v)==[1 n]) || any(v<=0))
                error('seizmo:plotmt:badInput',...
                    'LWIDTH must be a positive real valued number!');
            end
            if(numel(v)==1); v=v(ones(1,n),1); end
            pv.lwidth=v;
        case {'r' 'rad' 'radius'}
            if(isempty(v)); continue; end
            if(~isnumeric(v) || ~isreal(v) || ndims(v)~=2 ...
                    || ~any(numel(v)==[1 n]))
                error('seizmo:plotmt:badInput',...
                    'RADIUS must be a real valued number!');
            end
            if(numel(v)==1); v=v(ones(1,n),1); end
            pv.radius=v;
        case {'s' 'sample' 'd' 'delta' 'i' 'interval'}
            if(isempty(v)); continue; end
            if(~isnumeric(v) || ~isreal(v) || ndims(v)~=2 ...
                    || ~any(numel(v)==[1 n]))
                error('seizmo:plotmt:badInput',...
                    'SAMPLE must be a real valued number!');
            end
            if(numel(v)==1); v=v(ones(1,n),1); end
            pv.sample=v;
        case {'p' 'pa' 'parent' 'a' 'ax' 'axis' 'axes'}
            if(isempty(v)); pv.parent=[]; continue; end
            if(~isnumeric(v) || ~isreal(v) || ~isscalar(v) ...
                    || any(~ishandle(v)))
                error('seizmo:plotmt:badInput',...
                    'PARENT must be a valid axis handle!');
            end
            pv.parent=v;
    end
end

end
