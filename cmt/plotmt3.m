function [varargout]=plotmt3(x,y,z,mt,varargin)
%PLOTMT3    Plot moment tensor(s) as 3D beach ball(s)
%
%    Usage:    plotmt3(x,y,z,mt)
%              plotmt3(x,y,z,mt,'param1',val1,...,'paramN',valN)
%              h=plotmt3(...)
%
%    Description:
%     PLOTMT3(X,Y,Z,MT) plots the moment tensor(s) given by MT at the
%     location(s) specified with X, Y, & Z in the current plot axes as 3D
%     beach balls.  X, Y, & Z are assumed to be East, North, & Up making
%     them equivalent to roughly longitude, latitude, and elevation.  MT
%     must be in Harvard convention with one of the following formats:
%       - Nx6    [Mrr Mtt Mpp Mrt Mrp Mtp]
%       - 3x3xN  [Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp]
%       - scalar struct with fields .mrr .mtt .mpp .mrt .mrp .mtp
%     where N is the number of moment tensors to plot.  Note that X, Y, & Z
%     must be scalar or have N elements.
%
%     PLOTMT3(X,Y,Z,MT,'PARAM1',VAL1,...,'PARAMN',VALN) sets parameters
%     used in plotting the beach balls.  Valid parameters are:
%
%        NAME         DESCRIPTION
%      -----------+----------------------------------------
%      'draw'     | moment tensor related objects to draw
%      'tcolor'   | color of tension region (contains pressure axis)
%      'pcolor'   | color of compression region (contains tension axis)
%      'radius'   | radius of beach balls
%      'parent'   | handle of axes to plot in
%      'lcolor'   | color of nodal lines
%      'lwidth'   | width of nodal lines
%      'scolor'   | color of s-wave vector lines
%      'ppcolor'  | color of positive p-wave vector lines
%      'npcolor'  | color of negative p-wave vector lines
%      'enucolor' | color of coordinate axes
%      'tpbcolor' | color of principal axes
%      'font...'  | adjust font properties for 'enu' & 'tpb' axes labels
%                 | (they are passed directly to the function TEXT)
%      'sample'   | degree step size in calculating values in beach balls
%      'usample'  | degree step size in calculating displacement vectors
%
%         NAME         INPUT      MULTI     DEFAULT           ISSUES
%      -----------+--------------+-----+---------------+------------------
%      'draw'     | 'mt' 'nodal' |  N  | 'mt'          | Put mult. entries
%                 | 'enu' 'tpb'  |     |               | into a cell array
%                 |   'p' 's'    |     |               | 
%      'tcolor'   |   [r g b]    |  N  | [1 1 1]       | 
%      'pcolor'   |   [r g b]    |  N  | [0 0 0]       | 
%      'radius'   |   radius     |  Y  | .5            | 
%      'parent'   |   axhandle   |  N  | gca           | 
%      'lcolor'   |   [r g b]    |  Y  | [.5 .5 .5]    | 
%      'lwidth'   |    width     |  Y  | 3             | 
%      'scolor'   |   [r g b]    |  Y  | [1 0 1]       | 
%      'ppcolor'  |   [r g b]    |  Y  | [1 0 0]       | 
%      'npcolor'  |   [r g b]    |  Y  | [0 0 1]       | 
%      'enucolor' |   [r g b]    |  Y  | [0 1 0]       | 
%      'tpbcolor' |   [r g b]    |  Y  | [0 1 1]       | 
%      'font...'  |     N/A      |  N  | none          | See text props...
%      'sample'   |   degrees    |  Y  | 5             | 
%      'usample'  |   degrees    |  Y  | 20            | 
%
%     H=PLOTMT3(...) returns the graphics handles of all the objects making
%     up the beach balls in the cell array H.  All objects for the 1st
%     beach ball are in H{1}, the 2nd in H{2}, and so on.
%
%    Notes:
%
%    Examples:
%     % Plot a 3x3 grid of moment tensors with the s-wave radiation:
%     [x,y]=meshgrid(1:3,1:3);
%     cmts=findcmt('n',9);
%     plotmt3(x(:),y(:),0,cmts,'draw',{'mt' 's'});
%
%     % Plot a moment tensor with everything + the font set to size 16:
%     cmt=findcmt('n',1);
%     plotmt3(0,0,0,cmt,'draw',{'mt' 'nodal' 'enu' 'tpb' 'p' 's'},...
%         'fontsize',16);
%
%     % Some midwest quakes:
%     cmts=findcmts('latrange',[30 50],'lonrange',[-105 -85]);
%     plotmt3(cmts.centroidlon,cmts.centroidlat,-cmts.centroiddep,cmts);
%
%    See also: PLOTMT, RADPAT, FINDCMT, FINDCMTS

%     Version History:
%        Mar. 22, 2013 - initial version
%        Mar. 25, 2013 - update for mt_check/mt_change
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:
% - individual color schemes for mt

% check nargin
error(nargchk(4,inf,nargin));
if(mod(nargin,2))
    error('seizmo:plotmt3:badInput',...
        'PLOTMT3 call has an unpaired parameter/value!');
elseif(nargin>4 && ~iscellstr(varargin(1:2:end)))
    error('seizmo:plotmt3:badParam',...
        'PLOTMT3 parameters must be specified as strings!');
end

% check x/y/z/mt
if(~isnumeric(x) || ~isreal(x) || ~isnumeric(y) || ~isreal(y) ...
        || ~isnumeric(z) || ~isreal(z))
    error('seizmo:plotmt3:badInput',...
        'X/Y/Z must be real-valued numeric arrays!');
end
error(mt_check(mt));
mt=mt_change('v',mt);
n=size(mt,1);
if(~all([numel(x) numel(y) numel(z)]==1 | [numel(x) numel(y) numel(z)]==n))
    error('seizmo:plotmt3:badInput',...
        'X/Y/Z must be scalar or arrays with N elements!');
end

% expand x/y/z to Nx1
if(isscalar(x)); x(1:n)=x; end
if(isscalar(y)); y(1:n)=y; end
if(isscalar(z)); z(1:n)=z; end
[x,y,z]=deal(x(:),y(:),z(:));

% constant
D2R=pi/180;

% check parameters
[pv,varargin]=check_plotmt3_parameters(n,varargin{:});

% set handle to current axis if none give
if(isempty(pv.parent)); pv.parent=gca; end

% get hold state
held=ishold(pv.parent);
hflag=held;

% loop over each moment tensor
h=cell(n,1);
for i=1:n
    % draw!
    if(ismember('mt',lower(pv.draw)))
        [xyz,tri]=sph_tri_auto(pv.sample(i));
        xyz=xyz*pv.radius(i);
        [x0,y0,z0]=deal(xyz(:,1),xyz(:,2),xyz(:,3));
        [theta,phi]=cart2sph(x0,y0,z0);
        [theta,phi]=deal(pi/2-theta,pi/2+phi); % az,inc
        u_p=radpat(mt(i,:),phi/D2R,theta/D2R,'p');
        h{i}=cat(1,h{i},trisurf(tri,x(i)+x0,y(i)+y0,z(i)+z0,sign(u_p),...
            'parent',pv.parent,'tag','plotmt3_mt')); % this uses colormap!!
        if(i==1 && ~hflag); hold(pv.parent,'on'); end
    end
    if(ismember('enu',lower(pv.draw)))
        [u,v,w]=deal([2 0 0]*pv.radius(i),...
            [0 2 0]*pv.radius(i),[0 0 2]*pv.radius(i));
        tmp=quiver3(pv.parent,x([i i i],1),y([i i i],1),z([i i i],1),...
            u(:),v(:),w(:),0);
        set(tmp,'color',pv.enucolor(i,:),'tag','plotmt3_enu');
        h{i}=cat(1,h{i},tmp);
        if(i==1 && ~hflag); hold(pv.parent,'on'); end
        h{i}=cat(1,h{i},text(x(i)+u(1),y(i)+v(1),z(i)+w(1),'E',...
            'parent',pv.parent,'tag','plotmt3_enu_label',varargin{:}));
        h{i}=cat(1,h{i},text(x(i)+u(2),y(i)+v(2),z(i)+w(2),'N',...
            'parent',pv.parent,'tag','plotmt3_enu_label',varargin{:}));
        h{i}=cat(1,h{i},text(x(i)+u(3),y(i)+v(3),z(i)+w(3),'U',...
            'parent',pv.parent,'tag','plotmt3_enu_label',varargin{:}));
    end
    if(ismember('tpb',lower(pv.draw)))
        [t,p,b]=mt2tpb(mt(i,:));
        [u,v,w]=sph2cart((90-[t(3) p(3) b(3)])*D2R,...
            -[t(2) p(2) b(2)]*D2R,2*pv.radius(i));
        tmp=quiver3(pv.parent,x([i i i],1),y([i i i],1),z([i i i],1),...
            u(:),v(:),w(:),0);
        set(tmp,'color',pv.tpbcolor(i,:),'tag','plotmt3_tpb');
        h{i}=cat(1,h{i},tmp);
        if(i==1 && ~hflag); hold(pv.parent,'on'); end
        h{i}=cat(1,h{i},text(x(i)+u(1),y(i)+v(1),z(i)+w(1),'T',...
            'parent',pv.parent,'tag','plotmt3_tpb_label',varargin{:}));
        h{i}=cat(1,h{i},text(x(i)+u(2),y(i)+v(2),z(i)+w(2),'P',...
            'parent',pv.parent,'tag','plotmt3_tpb_label',varargin{:}));
        h{i}=cat(1,h{i},text(x(i)+u(3),y(i)+v(3),z(i)+w(3),'B',...
            'parent',pv.parent,'tag','plotmt3_tpb_label',varargin{:}));
    end
    if(ismember('nodal',lower(pv.draw)))
        [neu1,neu2]=nodallines(mt2sdr(mt(i,:)));
        neu1=neu1*pv.radius(i);
        neu2=neu2*pv.radius(i);
        h{i}=cat(1,h{i},plot3(pv.parent,x(i)+neu1(:,:,2),...
            y(i)+neu1(:,:,1),z(i)+neu1(:,:,3),'color',pv.lcolor(i,:),...
            'linewidth',pv.lwidth(i),'tag','plotmt3_nodal'));
        if(i==1 && ~hflag); hold(pv.parent,'on'); end
        h{i}=cat(1,h{i},plot3(pv.parent,x(i)+neu2(:,:,2),...
            y(i)+neu2(:,:,1),z(i)+neu2(:,:,3),'color',pv.lcolor(i,:),...
            'linewidth',pv.lwidth(i),'tag','plotmt3_nodal'));
    end
    if(ismember('p',lower(pv.draw)))
        xyz=sph_tri_auto(pv.usample(i))*pv.radius(i);
        [x1,y1,z1]=deal(xyz(:,1),xyz(:,2),xyz(:,3));
        [theta,phi]=cart2sph(x1,y1,z1);
        u_p=radpat(mt(i,:),(pi/2+phi)/D2R,(pi/2-theta)/D2R,'p');
        u_p=reshape(u_p,size(theta));
        n=u_p<0; p=~n;
        u_p=abs(u_p);
        [u1,v1,w1]=deal(x1.*u_p,y1.*u_p,z1.*u_p);
        [x1(n),y1(n),z1(n),u1(n),v1(n),w1(n)]=... % flip negative vectors
            deal(x1(n)+u1(n),y1(n)+v1(n),z1(n)+w1(n),-u1(n),-v1(n),-w1(n));
        tmp=quiver3(pv.parent,x(i)+x1(p),y(i)+y1(p),z(i)+z1(p),...
            u1(p),v1(p),w1(p),0);
        set(tmp,'color',pv.ppcolor(i,:),'tag','plotmt3_pp');
        h{i}=cat(1,h{i},tmp);
        if(i==1 && ~hflag); hold(pv.parent,'on'); end
        tmp=quiver3(pv.parent,x(i)+x1(n),y(i)+y1(n),z(i)+z1(n),...
            u1(n),v1(n),w1(n),0);
        set(tmp,'color',pv.npcolor(i,:),'tag','plotmt3_np');
        h{i}=cat(1,h{i},tmp);
    end
    if(ismember('s',lower(pv.draw)))
        xyz=sph_tri_auto(pv.usample(i))*pv.radius(i);
        [x1,y1,z1]=deal(xyz(:,1),xyz(:,2),xyz(:,3));
        [theta,phi]=cart2sph(x1,y1,z1);
        [theta,phi]=deal(pi/2-theta,pi/2+phi); % az,inc
        u=radpat(mt(i,:),phi/D2R,theta/D2R);
        u=reshape(u,[size(theta) 3])*2*pi*pv.radius(i)*pv.usample(i)/360;
        [svx,svy,svz]=deal(cos(phi).*sin(theta),...
            cos(phi).*cos(theta),sin(phi));
        [shx,shy,shz]=deal(cos(theta),-sin(theta),0*theta);
        [u1,v1,w1]=deal(svx.*u(:,:,2)+shx.*u(:,:,3),...
            svy.*u(:,:,2)+shy.*u(:,:,3),...
            svz.*u(:,:,2)+shz.*u(:,:,3));
        tmp=quiver3(pv.parent,x(i)+x1,y(i)+y1,z(i)+z1,u1,v1,w1,0);
        set(tmp,'color',pv.scolor(i,:),'tag','plotmt3_s');
        h{i}=cat(1,h{i},tmp);
        if(i==1 && ~hflag); hold(pv.parent,'on'); end
    end
end

% moment tensor surface adjustment
shading(pv.parent,'interp');
colormap(pv.parent,[pv.tcolor(1,:); pv.pcolor(1,:)]); % TEMP FIX!!!
axis(pv.parent,'equal','tight','vis3d');
rotate3d(pv.parent);

% restore hold
if(~held); hold(pv.parent,'off'); end

% output if desired
if(nargout); varargout{1}=h; end

end


function [pv,varargin]=check_plotmt3_parameters(n,varargin)
%CHECK_PLOTMT3_PARAMETERS    Default/Parse/Check PLOTMT3 parameters
%
%         NAME         INPUT      MULTI     DEFAULT           ISSUES
%      -----------+--------------+-----+---------------+------------------
%      'draw'     | 'mt' 'nodal' |  N  | 'mt'          | Put mult. entries
%                 | 'enu' 'tpb'  |     |               | into a cell array
%                 |   'p' 's'    |     |               | 
%      'tcolor'   |   [r g b]    |  Y  | [1 1 1]       | 
%      'pcolor'   |   [r g b]    |  Y  | [0 0 0]       | 
%      'radius'   |   radius     |  Y  | .5            | 
%      'parent'   |   axhandle   |  N  | gca           | 
%      'lcolor'   |   [r g b]    |  Y  | [.5 .5 .5]    | 
%      'lwidth'   |    width     |  Y  | 3             | 
%      'scolor'   |   [r g b]    |  Y  | [1 0 1]       | 
%      'ppcolor'  |   [r g b]    |  Y  | [1 0 0]       | 
%      'npcolor'  |   [r g b]    |  Y  | [0 0 1]       | 
%      'enucolor' |   [r g b]    |  Y  | [0 1 0]       | 
%      'tpbcolor' |   [r g b]    |  Y  | [0 1 1]       | 
%      'font...'  |     N/A      |  N  | none          | See text props...
%      'sample'   |   degrees    |  Y  | 5             | 
%      'usample'  |   degrees    |  Y  | 20            | 

% defaults
varargin=[{'draw' 'mt' 'tcolor' [1 1 1] 'pcolor' [0 0 0] 'radius' .5 ...
    'parent' [] 'lcolor' [.5 .5 .5] 'lwidth' 3 'scolor' [1 0 1] ...
    'ppcolor' [1 0 0] 'npcolor' [0 0 1] 'enucolor' [0 1 0] ...
    'tpbcolor' [0 1 1] 'sample' 5 'usample' 20} varargin];

% valid draw values
valid={'mt' 'nodal' 'enu' 'tpb' 'p' 's'};

% check/parse given pv pairs
keep=true(numel(varargin),1);
for i=1:2:numel(varargin)
    p=varargin{i};
    v=varargin{i+1};
    switch lower(p)
        case {'draw' 'dr' 'd'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if(ischar(v)); v=cellstr(v); end
            if(~iscellstr(v) || any(~ismember(lower(v),valid)))
                error('seizmo:plotmt3:badInput',...
                    ['DRAW must be one or more of the following:\n' ...
                    sprintf('%s ',valid{:})]);
            end
            pv.draw=v;
        case {'t' 'tc' 'tcolor'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'TCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.tcolor=v;
        case {'c' 'pc' 'pcolor' 'color'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'PCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.pcolor=v;
        case {'lc' 'lcolor' 'nodalcolor'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'LCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.lcolor=v;
        case {'sc' 'scolor' 'swavecolor'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'SCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.scolor=v;
        case {'ppc' 'pos' 'ppcolor' 'poscolor'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'PPCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.ppcolor=v;
        case {'npc' 'neg' 'npcolor' 'negcolor'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'NPCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.npcolor=v;
        case {'enuc' 'enu' 'enucolor'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'ENUCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.enucolor=v;
        case {'tpbc' 'tpb' 'tpbcolor'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if((~ischar(v) ...
                    && (~isnumeric(v) || ~isreal(v) || size(v,2)~=3)) ...
                    || ndims(v)~=2 || ~any(size(v,1)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'TPBCOLOR must be given as [RED GREEN BLUE]!');
            end
            if(size(v,1)==1); v=v(ones(1,n),:); end
            if(ischar(v)); v=name2rgb(v); end
            pv.tpbcolor=v;
        case {'lw' 'lwidth' 'linewidth'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if(~isnumeric(v) || ~isreal(v) || ndims(v)~=2 ...
                    || ~any(numel(v)==[1 n]) || any(v<=0))
                error('seizmo:plotmt3:badInput',...
                    'LWIDTH must be a positive real valued number!');
            end
            if(numel(v)==1); v=v(ones(1,n),1); end
            pv.lwidth=v;
        case {'r' 'rad' 'radius'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if(~isnumeric(v) || ~isreal(v) || ndims(v)~=2 ...
                    || ~any(numel(v)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'RADIUS must be a real valued number!');
            end
            if(numel(v)==1); v=v(ones(1,n),1); end
            pv.radius=v;
        case {'s' 'sample' 'i' 'interval'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if(~isnumeric(v) || ~isreal(v) || ndims(v)~=2 ...
                    || ~any(numel(v)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'SAMPLE must be a real valued number!');
            end
            if(numel(v)==1); v=v(ones(1,n),1); end
            pv.sample=v;
        case {'u' 'us' 'usample' 'ui' 'uinterval'}
            keep(i:i+1)=false;
            if(isempty(v)); continue; end
            if(~isnumeric(v) || ~isreal(v) || ndims(v)~=2 ...
                    || ~any(numel(v)==[1 n]))
                error('seizmo:plotmt3:badInput',...
                    'USAMPLE must be a real valued number!');
            end
            if(numel(v)==1); v=v(ones(1,n),1); end
            pv.usample=v;
        case {'p' 'pa' 'parent' 'a' 'ax' 'axis' 'axes'}
            keep(i:i+1)=false;
            if(isempty(v)); pv.parent=[]; continue; end
            if(~isnumeric(v) || ~isreal(v) || ~isscalar(v) ...
                    || any(~ishandle(v)))
                error('seizmo:plotmt3:badInput',...
                    'PARENT must be a valid axis handle!');
            end
            pv.parent=v;
    end
end

% clear out parsed pv pairs
varargin=varargin(keep);

end

