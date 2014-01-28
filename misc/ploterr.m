function [varargout]=ploterr(varargin)
%PLOTERR    General errorbar plot
%
%    Usage:    ploterr(y)  &  ploterr(x,y)
%              ploterr(x,y,xerr,yerr)
%              ploterr(x,y,{lx ux},{ly uy})
%              ploterr(x,y,[],yerr)  &  ploterr(x,y,xerr,[])
%              ploterr(x,y,xerr,{ly uy})  &  ploterr(x,y,{lx ux},yerr)
%              ploterr(x,y,xerr,yerr,linestyle)
%              ploterr(...,'property1',value1,'property2',value2,...)
%              ploterr(ax,...)
%              h=ploterr(...)
%              [h,hyerr,hxerr]=ploterr(...)
%
%    Description:
%     PLOTERR(Y) & PLOTERR(X,Y) are the same as PLOT(Y) & PLOT(X,Y).
%
%     PLOTERR(X,Y,XERR,YERR) plots X vs. Y with x errorbars [X-XERR X+XERR]
%     and y errorbars [Y-YERR Y+YERR].
%
%     PLOTERR(X,Y,{LX UX},{LY UY}) plots X vs. Y with errorbars specified
%     by LX and UX in horizontal direction and with errorbars specified by
%     LY and UY in vertical direction.  L and U contain the lower and upper
%     error ranges for each point in X resp. Y (L = lower, U = upper).
%     Each error bar ranges from L(i) to U(i). X, LX and UX sizes may vary
%     in singleton dimensions only, the same accounts for Y, LY and UY. If
%     any of X,Y,LX,UX,LY,UY is a matrix then each column produces a
%     separate line.
%
%     PLOTERR(X,Y,[],YERR) & PLOTERR(X,Y,XERR,[]) will skip plotting
%     the errorbars with no values (empty matrix).
%
%     PLOTERR(X,Y,XERR,{LY UY}) & PLOTERR(X,Y,{LX UX},YERR) allows mixed
%     types of errorbar specifications.
%
%     PLOTERR(X,Y,XERR,YERR,LINESPEC) uses the color and linestyle
%     specified by the string LINESPEC. See PLOT for possibilities.
%     Defaults to a solid line connecting the points (X,Y).  Errorbars only
%     utilize the color value in LINESPEC (forced to be solid lines).  Note
%     that X, Y, XERR, & YERR must all be input before LINESPEC and cannot
%     be omitted (use [] as placeholders if necessary).
%
%     PLOTERR(...,'PROPERTY1',VALUE1,'PROPERTY2',VALUE2,...) allows
%     changing specific properties of PLOTERR.  Property Value pairs must
%     be passed after LINESPEC, however, LINESPEC does not need to be
%     passed.  The Following properties are available:
%      'logx', 'logy', 'logxy'       - toggles for logarithmic scaling
%      'hhx', 'hhy', 'hhxy'          - relative size of handlebar height
%      'abshhx', 'abshhy', 'abshhxy' - absolute size of handlebar height
%
%     The 'log' properties do not require an accompanying value.  On their
%     own they are an "on" switch.  The values accompanying the 'log'
%     properties must evaluate to either TRUE or FALSE where TRUE turns the
%     log scale on and FALSE turns it off.
%
%     The default for 'hhx' and 'hhy' is 1/50, indicating a total width of
%     the x-handlebars is 1/50 of the distance range of datapoints in y and
%     vice-versa for the y-handlebars.  For logarithmic plots that is the
%     distance range on a logarithmic scale.
%
%     Absolute size of errorbar handles for log scale is in units of power
%     of 10 i.e. the command PLOTERR(...,'logy','abshhx',1.0) produces x
%     handlebars where each one spans a decade on the y axis.
%
%     PLOTERR(AX,...) plots into the axes with handle AX.
%
%     H=PLOTERR(...) returns a vector of line handles in the order:
%      H(1) = handle to datapoints
%      H(2) = handle to errorbar y OR errorbar x if error y not specified
%      H(3) = handle to errorbar x if error y specified
%     If more than one line is plotted, the ordering is the following:
%      H(1:n) = handle to lines with datapoints
%      H(n+1:2*n) = handle to y errorbars
%      H(2*n+1:3*n) = handle to x errorbars
%
%     [H,HYERR,HXERR]=PLOTERR(...) returns the line handles separately.
%
%    Notes:
%     - X, Y, xerr and yerr must be of the same size as long as dimensions
%       are not equal to 1, otherwise singletons are expanded. You can pass
%       everything you like, e.g. xerr={0, 20:29} and x=[(1:10)' (2:11)']
%       to indicate a lower bound of 0 for all values and an upper bound of
%       20:29 for both columns of x, which will show up as two separate
%       lines.
%     - If both abshh* and hh* are set, the latter of both in the argument
%       list is used. The same is true for any of the other properties
%       passed twice.
%     - Properties ending on 'xy' may also end on 'yx' or you can leave the
%       ending away: logxy = logyx = log
%
%    Examples:
%     % Draws symmetric errorbars of unit
%     % standard deviation along a sine wave:
%     x=2:11;
%     y=sin(x);
%     e=std(y)*ones(size(x));
%     ploterr(x,y,e);
%   
%     % Draws samples of a noisy exponential function on a logarithmic
%     % y-scale with constant relative errors in x and variable absolute
%     % errors in y with slim errorbar handles. The lineseries objects are
%     % used to set the color of the error bars and to display a legend.
%     x=0:15;
%     y=exp(-x).*(rand(1,16)*0.9+0.1);
%     h=ploterr(x,y,0.3,{exp(-x)*0.1 exp(-x)},'r.','logy','hhxy',1e-2);
%     set(h(2),'Color','b'), set(h(3),'Color','b');
%     legend(h,{'data' 'error x' 'error y'});
%
%     % Plot normal x errorbars and tiny y errorbars
%     % on a logarithmic y-scale and a linear x scale.
%     ploterr(x,y,ex,ey,'logy','hhy',1e-3);
%
%     % Show off ability to pass on prop/val pairs:
%     ploterr(rand(50,1),rand(50,1),rand(50,1)/10,rand(50,1)/10,'sk',...
%         'markerfacecolor','g','markersize',5,'linewidth',2);
%
%    See also: ERRORBAR, PLOT, LOGLOG, SEMILOGX, SEMILOGY

%     Version History:
%        Oct. 23, 2006 - initial version (errorbar_x)
%        Oct. --, 2008 - modification of errorbar_x by Goetz Huesken for
%                        plotting horizontal and/or vertical errorbars and
%                        to support logarithmic scaling.
%        Dec. --, 2008 - changed the user interface. Handle sizes and
%                        logarithmic scaling can now be set via properties.
%                        LineSpec is not compulsory anymore when setting
%                        logscale or handle sizes.
%        Dec. --, 2008 - bugfixes of previous version
%        Jan. --, 2009 - added abshh properties to set handle sizes as
%                        absolute values
%        Feb.  1, 2011 - revamped arg parsing to allow axes input at start
%                        or as a prop/val pair.  Allow passing prop/val
%                        pairs to PLOT.  Organized code & docs.  Output
%                        with handles split up possible now.  Also allow
%                        xy/yx to be at start or end of option.
%        Feb. 10, 2011 - make hhx/hhy/hhxy based on range not mean
%        Jan. 27, 2014 - use axparse instead of axescheck for octave
%
%     Written by Goetz Huesken (goetz.huesken(at)gmx.de)
%                Felix Z�rgiebel (felix_z -> web.de)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 13:35 GMT

% todo:

% parse inputs
[ax,x,y,xerr,yerr,sym,logx,logy,hhx,hhy,varargin]=parseinputs(varargin{:});

% use indices for x if no values given
if(isempty(x)); x=(1:size(y,1))'; end

% check & expand scalar dimensions of x/y
[x,y]=expandarrays(x,y,'X','Y');

% check if xerr is relative or absolute
relxerr=false;
if(~iscell(xerr))
    if(~isempty(xerr))
        relxerr=true;
        xerr={-xerr +xerr};
    else
        xerr={[] []};
    end
elseif(numel(xerr)~=2)
    error('fex:ploterr:badInput',...
        'XERR cell input must have two entries (lower & upper bounds)!');
end

% make xerr, x, y arrays have the same size
[xl,xh]=expandarrays(xerr{1},xerr{2},'LX','UX');
[y,xl,ty,txl]=expandarrays(y,xl,'Y','LX');
if(ty); y=y'; xl=xl'; txl=~txl; end % make sure x and y still match
if(txl); xh=xh'; end % make sure xl and xh still match
[y,xh]=expandarrays(y,xh,'Y','UX'); % propagate expansion to UX
[x,y]=expandarrays(x,y,'X','Y'); % propagate expansion to X
if(relxerr); xl=x+xl; xh=x+xh; end % now error is absolute

% check if yerr is relative or absolute
relyerr=false;
if(~iscell(yerr))
    if(~isempty(yerr))
        relyerr=true;
        yerr={-yerr +yerr};
    else
        yerr={[] []};
    end
elseif(numel(yerr)~=2)
    error('fex:ploterr:badInput',...
        'YERR cell input must have two entries (lower & upper bounds)!');
end

% make everyone have the same size
[yl,yh]=expandarrays(yerr{1},yerr{2},'LY','UY');
[x,yl,tx,tyl]=expandarrays(x,yl,'X','LY');
if(tx); x=x'; yl=yl'; tyl=~tyl; end % make sure x and y still match
if(tyl); yh=yh'; end % make sure yl and yh still match
[x,yh]=expandarrays(x,yh,'X','UY'); % propagate expansion to UY
[x,y]=expandarrays(x,y,'X','Y'); % propagate expansion to Y
[x,xl]=expandarrays(x,xl,'X','LX'); % propagate expansion to LX
[x,xh]=expandarrays(x,xh,'X','UX'); % propagate expansion to UX
if(relyerr); yl=y+yl; yh=y+yh; end % now error is absolute

% choose the appropriate function for the plot
if(     logx &&  logy); plotfct=@loglog;
elseif( logx && ~logy); plotfct=@semilogx;
elseif(~logx &&  logy); plotfct=@semilogy;
else                    plotfct=@plot;
end

% LineSpec setup
[ls,col,mark,msg]=colstyle(sym);
if(~isempty(msg)); error('fex:ploterr:badInput',msg); end
sym=[ls mark col]; % Use marker only on data part
esym=['-' col]; % force errorbars to be solid lines

% what is the hold state?
hold_state=false;
if(~isempty(ax)); hold_state=ishold(ax{:}); end

% default output handles
h=[]; xeh=[]; yeh=[];

% plot specified data
if(~isempty(xl)) % x errorbars
    [bary,barx]=barline(y,xl,xh,logy,hhx);
    xeh=plotfct(ax{:},barx,bary,esym,varargin{:});
    ax={get(xeh(1),'parent')};
    hold(ax{:},'on');
end
if(~isempty(yl)) % y errorbars
    [barx,bary]=barline(x,yl,yh,logx,hhy);
    yeh=plotfct(ax{:},barx,bary,esym,varargin{:});
    ax={get(yeh(1),'parent')};
    hold(ax{:},'on');
end
if(~isempty(y)) % function values
    h=plotfct(ax{:},x,y,sym,varargin{:});
    ax={get(h(1),'parent')};
end

% restore hold state
if(~hold_state); hold(ax{:},'off'); end

% output if desired
if(nargout==1)
    varargout={[h; yeh; xeh]};
elseif(nargout)
    varargout={h yeh xeh};
end

end


function [perp,para]=barline(v,l,h,uselog,handleheight)
%BARLINE    Creates errorbar line input for PLOT

% v: value "perpendicular"
% l: lower bound "parallel"
% h: upper bound "parallel"

% npts & nlines
[npt,n]=size(l);

% math depends on scaling
if(uselog)
    % basic operations for logarithmic spacing
    dist=@rdivide;
    invdist=@times;
    scale=@power;
else
    % basic operations for linear spacing
    dist=@minus;
    invdist=@plus;
    scale=@times;
end

% calculate height of errorbar delimiters (handlebars)
if(handleheight>0) % means handleheight was passed as a relative value
    % Set width of ends of bars to handleheight times distance range of
    % points.  If number of points is under 15, space as if 15 points
    % were there.
    if(dist(max(v(:)),min(v(:)))==0)
        dv=scale(abs(v),1/50)+(abs(v)==0);
    else
        dv=scale(dist(max(v(:)),min(v(:))),handleheight/2);
    end
else % handleheight<=0 means handleheight was passed as an absolute value
    dv=handleheight/2;
    if(uselog); dv=10^dv; end
end

% get end positions of handlebars
vh=invdist(v,dv);
vl=dist(v,dv);

% build up nan-separated vector for bars
para=zeros(npt*9,n);
para(1:9:end,:)=h;
para(2:9:end,:)=l;
para(3:9:end,:)=NaN;
para(4:9:end,:)=h;
para(5:9:end,:)=h;
para(6:9:end,:)=NaN;
para(7:9:end,:)=l;
para(8:9:end,:)=l;
para(9:9:end,:)=NaN;

perp=zeros(npt*9,n);
perp(1:9:end,:)=v;
perp(2:9:end,:)=v;
perp(3:9:end,:)=NaN;
perp(4:9:end,:)=vh;
perp(5:9:end,:)=vl;
perp(6:9:end,:)=NaN;
perp(7:9:end,:)=vh;
perp(8:9:end,:)=vl;
perp(9:9:end,:)=NaN;

end


function [A,B,tA,tB]=expandarrays(A,B,sA,sB)
%EXPANDARRAYS    Reorients & expands inputs to match dimensions

% tA,tB: indicate if A,B have been transposed

% sizes and transpose states
sizA=size(A); tA=false;
sizB=size(B); tB=false;

% do not process empty arrays
if(isempty(A) || isempty(B)); return; end

% require 2D arrays
if(ndims(A)~=2 || ndims(B)~=2)
    error('fex:ploterr:badInput',...
        '%s and %s must not have more than two dimensions!',sA,sB);
elseif(~isreal(A) || ~isreal(B))
    error('fex:ploterr:badInput',...
        '%s and %s must be real-valued!',sA,sB);
end

% make vectors column vectors
if(sizA(1)==1); A=A(:); tA=~tA; sizA=sizA([2 1]); end
if(sizB(1)==1); B=B(:); tB=~tB; sizB=sizB([2 1]); end

% transpose to fit column, if necessary
% - note B is transposed to match A first
if(sizA(2)==1 && sizB(2)~=1 && sizB(2)==sizA(1) && sizB(1)~=sizA(1))
    B=B'; tB=~tB; sizB=sizB([2 1]);
end
if(sizB(2)==1 && sizA(2)~=1 && sizA(2)==sizB(1) && sizB(1)~=sizB(1))
    A=A'; tA=~tA; sizA=sizA([2 1]);
end

% require dimensions are equal or scalar
if(all(sizA==sizB | sizA==1 | sizB==1))
    % expand A
    singletonsA=find(sizA==1 & sizB~=1);
    repA={':' ':'};
    for d=singletonsA; repA{d}=ones(sizB(d),1); end
    A=A(repA{:});
    
    % expand B
    singletonsB=find(sizB==1 & sizA~=1);
    repB={':' ':'};
    for d=singletonsB; repB{d}=ones(sizA(d),1); end
    B=B(repB{:});
else % otherwise return error
    error('fex:ploterr:badInput',...
        '%s & %s must be equally sized for non-singleton dimensions!',...
        sA,sB);
end

end


function [ax,x,y,xerr,yerr,sym,lx,ly,hx,hy,varargin]=parseinputs(varargin)
%PARSEINPUTS    Parses inputs passed to PLOTERR

% find axes handle if given
[ax,varargin]=axparse(varargin{:});
nargin=numel(varargin);
if(isempty(ax)); ax={}; else ax={ax}; end

% check nargin
error(nargchk(1,inf,nargin));

% extract inputs and strip from options
x=varargin{1}; varargin(1)=[];
if(nargin<2); y=x; x=[]; else y=varargin{1}; varargin(1)=[]; end
if(nargin<3); xerr=[]; else xerr=varargin{1}; varargin(1)=[]; end
if(nargin<4); yerr=[]; else yerr=varargin{1}; varargin(1)=[]; end

% option defaults
sym='-';
lx=0;
ly=0;
hx=1/50;
hy=1/50;

% skip rest if nothing left
if(nargin<=4); return; end

% extract leading linestyle arg
[ls,col,mark,errmsg]=colstyle(varargin{1});
if(isempty(errmsg))
    sym=[ls col mark];        % first entry is a linestyle
    varargin=varargin(2:end); % skip symbol for properties
end

% loop over remaining options
nvarg=numel(varargin); idx=1; keep=true(nvarg,1);
while(idx<=nvarg)
    % get current prop/value pair
    prop=varargin{idx};
    val=varargin{idx+1};
    
    % require properties to be strings
    if(~ischar(prop))
        error('fex:ploterr:badInput',...
            'Input PROPERTIES must be strings!');
    end
    
    % look for ploterr properties
    switch lower(prop)
        case {'xlog' 'logx'}
            if(islogical(val)); val=double(val); end
            if(isnumeric(val)); lx=val; keep(idx:idx+1)=false;
            else lx=1; keep(idx)=false; idx=idx-1;
            end
        case {'ylog' 'logy'}
            if(islogical(val)); val=double(val); end
            if(isnumeric(val)); ly=val; keep(idx:idx+1)=false;
            else ly=1; keep(idx)=false; idx=idx-1;
            end
        case {'xylog' 'yxlog' 'log' 'logxy' 'logyx'}
            if(islogical(val)); val=double(val); end
            if(isnumeric(val)); ly=val; lx=val; keep(idx:idx+1)=false;
            else ly=1; lx=1; keep(idx)=false; idx=idx-1;
            end
        case {'xhh' 'hhx'}
            if(isnumeric(val))
                hx=abs(val); keep(idx:idx+1)=false;
            else
                error('fex:ploterr:badInput',...
                    'HHX must be followed by a numerical value!');
            end
        case {'yhh' 'hhy'}
            if(isnumeric(val))
                hy=abs(val); keep(idx:idx+1)=false;
            else
                error('fex:ploterr:badInput',...
                    'HHY must be followed by a numerical value!');
            end
        case {'xyhh' 'yxhh' 'hh' 'hhxy' 'hhyx'}
            if(isnumeric(val))
                hy=abs(val); hx=abs(val); keep(idx:idx+1)=false;
            else
                error('fex:ploterr:badInput',...
                    'HH must be followed by a numerical value!');
            end
        case {'xabshh' 'xhhabs' 'abshhx' 'hhabsx'}
            if(isnumeric(val))
                hx=-abs(val); keep(idx:idx+1)=false;
            else
                error('fex:ploterr:badInput',...
                    'ABSHHX must be followed by a numerical value!');
            end
        case {'yabshh' 'yhhabs' 'abshhy' 'hhabsy'}
            if(isnumeric(val))
                hy=-abs(val); keep(idx:idx+1)=false;
            else
                error('fex:ploterr:badInput',...
                    'ABSHHY must be followed by a numerical value!');
            end
        case {'hhabs' 'abshh' 'abshhxy' 'abshhyx' 'hhabsxy' 'hhabsyx' ...
                'xyabshh' 'yxabshh' 'xyhhabs' 'yxhhabs'}
            if(isnumeric(val))
                hy=-abs(val); hx=-abs(val); keep(idx:idx+1)=false;
            else
                error('fex:ploterr:badInput',...
                    'ABSHH must be followed by a numerical value!');
            end
    end
    idx=idx+2;
end

% remove ploterr properties
varargin(~keep)=[];

end

