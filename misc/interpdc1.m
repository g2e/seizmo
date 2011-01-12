function [varargout]=interpdc1(x,y,varargin)
%INTERPDC1    1D interpolation (table lookup) with discontinuity support
%
%    Usage:    [yi_hi,yi_lo]=interpdc1(x,y,xi)
%              [yi_hi,yi_lo]=interpdc1(x,y,xi,method)
%              [yi_hi,yi_lo]=interpdc1(x,y,xi,method,'extrap')
%              [yi_hi,yi_lo]=interpdc1(x,y,xi,method,extrapval)
%              [yi_hi,yi_lo]=interpdc1(x,y,xi,method,extrapval,clamp)
%              pp=interpdc1(x,y,method,'pp')
%              pp=interpdc1(x,y,method,'pp',clamp)
%
%    Description:
%     [YI_HI,YI_LO]=INTERPDC1(X,Y,XI) returns values linearly interpolated
%     at positions XI utilizing the data Y at positions X.  X must be a
%     vector and may include repeat points for defining discontinuities in
%     Y (note that 3 or more points in X can not share the same position
%     but there may be multiple pairs defining multiple discontinuities at
%     different positions).  Y may be a vector (equal length to X) or an
%     array (same number of rows as elements in X).  Interpolation is done
%     across the rows (ie along dimension 1) in Y if it is an array.  YI_HI
%     & YI_LO both give the values at each point in XI and so in general
%     are equal except for points in XI that are on a discontinuity.  For
%     those points YI_HI contains the values on the "upside" of the
%     discontinuity and YI_LO contains the values on the "downside" of the
%     discontinuity.  If Y is a vector, XI, YI_HI & YI_LO will be the same
%     size.  Otherwise, if XI is an array YI_HI & YI_LO have dimensions
%     [SIZE(XI) SIZE_Y(2:END)] or if XI is a vector YI_HI & YI_LO have
%     dimensions [LENGTH(XI) SIZE_Y(2:END)].
%
%     [YI_HI,YI_LO]=INTERPDC1(X,Y,XI,METHOD) sets the interpolation method.
%     May be any of the following:
%      'nearest' - nearest neighbor interp
%      'linear'  - linear interp (THE DEFAULT)
%      'spline'  - piecewise cubic spline interp
%      'pchip'   - shape-preserving piecewise cubic interp
%      'cubic'   - same as 'pchip'
%      'v5cubic' - the cubic interpolation from MATLAB 5, which does not
%                  extrapolate (by default) and uses 'spline' if X is not
%                  equally spaced.
%
%     [YI_HI,YI_LO]=INTERPDC1(X,Y,XI,METHOD,'EXTRAP') extrapolates using
%     METHOD to define values outside the range given by X.  This is the
%     default for METHODs 'spline', 'pchip' & 'cubic'.
%
%     [YI_HI,YI_LO]=INTERPDC1(X,Y,XI,METHOD,EXTRAPVAL) replaces values
%     outside the range spanned by X with EXTRAPVAL.  EXTRAPVAL=NaN is the
%     default for METHODs 'nearest', 'linear' & 'v5cubic'.
%
%     [YI_HI,YI_LO]=INTERPDC1(X,Y,XI,METHOD,EXTRAPVAL,CLAMP) forces the
%     edges of individual sections to match the slopes given in CLAMP.
%     This option is only implemented when METHOD is 'spline'.  CLAMP must
%     be a cell array of size NREGx1, where NREG is the number of
%     continuous sections in X.  Each cell in CLAMP should match dimensions
%     2+ of Y, so it may be concatenated to it along dimension 1.  CLAMP is
%     expected to have 2 rows, the first gives the slope at the lower end
%     of the section and the second row gives the slope at the upper end of
%     the section. CLAMP may be a scalar cell to force equal slope boundary
%     conditions to all continuous sections.
%
%     PP=INTERPDC1(X,Y,METHOD,'PP') uses the specified METHOD to create
%     a piecewise polynomial struct PP.  PPVAL(PP,XI) is equivalent to
%     PPDCVAL(PP,XI) & INTERPDC1(X,Y,XI,METHOD,'EXTRAP') first output.
%     PPDCVAL(PP,XI,FALSE) is equal to the second output of
%     INTERPDC1(X,Y,XI,METHOD,'EXTRAP').  This is useful when you want to
%     call INTERPDC1 on the same X & Y for several XI sets.  Note that the
%     piecewise polynomial spans all discontinuities.
%
%     PP=INTERPDC1(X,Y,METHOD,'PP',CLAMP) returns a piecewise polynomial
%     struct that has clamped edges (see description above for details on
%     CLAMP).
%
%    Notes:
%
%    Examples:
%     Create a dataset with discontinuities:
%      x=[1:10 10:15 15:27 27:40]; % discon at 10, 15, 27
%      y=[rand(1,10) rand(1,6)+2 rand(1,13)-5 rand(1,14)];
%      xi=0:0.1:41;
%      [yihi,yilo]=interpdc1(x,y,xi,'spline');
%      xi=[xi; xi]; yi=[yilo; yihi];
%      figure; plot(xi(:),yi(:));
%
%     Clamping in action (all sections must have slopes of 0 at the edges):
%      x=[1:10 10:15 15:27 27:40]'; % discon at 10, 15, 27
%      y=[rand(1,10) rand(1,6)+2 rand(1,13)-5 rand(1,14)];
%      xi=0:0.1:41;
%      [yihi,yilo]=interpdc1(x,y,xi,'spline','',{[0; 0]});
%      xi=[xi; xi]; yi=[yilo; yihi];
%      figure; plot(xi(:),yi(:),'g');
%      [yihi,yilo]=interpdc1(x,y,xi,'spline');
%      xi=[xi; xi]; yi=[yilo; yihi];
%      hold on; plot(xi(:),yi(:),'r');
%
%    See also: PPDCVAL, INTERP1, INTERP1Q, MKPP, PPVAL

%     Version History:
%        May  18, 2010 - initial version
%        May  19, 2010 - properly handle non-increasing case
%        June  1, 2010 - finally true discontinuity support
%        Aug.  8, 2010 - doc update
%        Jan. 11, 2011 - do not compute negative side if < 2 outputs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 11, 2011 at 03:20 GMT

% todo:

% check nargin
error(nargchk(3,6,nargin));

% check x/y
% - basic shape/type req
if(~isnumeric(x) || ~isnumeric(y))
    error('seizmo:interpdc1:badXY',...
        'X & Y must be numeric arrays!');
elseif(~isvector(x) || ~isreal(x))
    error('seizmo:interpdc1:badXY',...
        'X must be a real-valued vector!');
end
% - get y oriented
sy=size(y);
if(isvector(y))
    % row to column vector
    y=y(:);
    ny=prod(sy);
    yv=true;
    sy2=1;
    psy2=1;
else % array
    % flatten ND arrays
    y=y(:,:);
    ny=sy(1);
    yv=false;
    sy2=sy(2:end);
    psy2=prod(sy2);
end
% - require 2+ points
if(ny<2)
    error('seizmo:interpdc1:badXY',...
        'INTERPDC1 requires 2+ data points!');
end
% - x/y compat
if(numel(x)~=ny)
    if(yv)
        error('seizmo:interpdc1:badXY',...
            'X & Y must be the same length!');
    else
        error('seizmo:interpdc1:badXY',...
            'LENGTH(X) & SIZE(Y,1) must be the same!');
    end
end
x=x(:);
% - nan checks
if(any(isnan(x)))
    % no x nans allowed
    error('seizmo:interpdc1:nanX',...
        'X can not have NaNs!');
elseif(any(isnan(y(:))))
    % y nans are not interpolated over
    warning('seizmo:interpdc1:nanY',...
        'NaN Y values found! NaNs are not removed/interpolated!');
end

% breaking into regions
% - first require monotonicity of x
dx=diff(x);
if(all(dx<=0))
    % flip x/y so monotonically non-decreasing
    x=flipud(x);
    y=flipud(y);
    dx=-flipud(dx);
elseif(~all(dx>=0))
    % x must be non-decreasing
    error('seizmo:interpdc1:nonmonoX',...
        'X must be monotonically non-decreasing!');
end
% - find discontinuities
dctop=find(dx==0);
nreg=numel(dctop)+1;
% - next disallow 3+ equal x
if(nreg>2 && any(diff(dctop)==1))
    error('seizmo:interpdc1:badDiscon',...
        'X discontinuity must only have 2 points!');
end
% - indexing to break x/y into regions
ridx=[[1; dctop+1] [dctop; ny]];
rn=diff(ridx,1,2)+1;
if(any(rn<2))
    error('seizmo:interpdc1:needMoreData',...
        'Each continuous section of X/Y must have 2+ points!');
end

% pp or yi?
% - parsing inputs here b/c of differences
method='linear';
clamp=cell(nreg,1);
%top=true;
if(ischar(varargin{1}))
    % outputing pp
    method=varargin{1};
    ppout=true;
    if(nargin==3 || ~strcmpi('pp',varargin{2}))
        error('seizmo:interpdc1:undefinedCall',...
            'When METHOD is 3rd argument ''pp'' option must follow!');
    end
    % check for clamp
    if(nargin==5 && ~isempty(varargin{3})); clamp=varargin{3}; end
else
    % outputing yi
    ppout=false;
    xi=varargin{1};
    if(yv)
        syi=size(xi);
    else
        if(isvector(xi))
            syi=[numel(xi) sy2];
        else
            syi=[size(xi) sy2];
        end
    end
    xi=xi(:);
    % check xi
    if(~isreal(xi))
        error('seizmo:interpdc1:badXI',...
            'XI must be a real-valued array!');
    end
    % get remaining options if available
    if(nargin>3 && ~isempty(varargin{2})); method=varargin{2}; end
    if(nargin>4 && ~isempty(varargin{3})); extrap=varargin{3}; end
    if(nargin>5 && ~isempty(varargin{4})); clamp=varargin{4}; end
end

% check method
valid.methods={'nearest' 'linear' 'spline' 'pchip' 'cubic' 'v5cubic'};
if(~ischar(method) || isempty(strmatch(method,valid.methods)))
    error('seizmo:interpdc1:badMethod',...
        ['METHOD must be one of the following:\n' ...
        sprintf('''%s'' ',valid.methods)]);
end
methidx=strmatch(method,valid.methods);

% check clamp
% - note that clamp is ignored for all but spline
if(~iscell(clamp) || (~isscalar(clamp) && numel(clamp)~=nreg))
    error('seizmo:interpdc1:badCLAMP',...
        'CLAMP must be a cell array with one cell per region!');
end
if(isscalar(clamp)); clamp(1:nreg,1)=clamp; end

% create pp
% - do last first to preallocate struct
pp(nreg)=ppmake(x(ridx(nreg,1):ridx(nreg,2)),...
    y(ridx(nreg,1):ridx(nreg,2),:),...
    methidx,rn(nreg),sy2,psy2,clamp{nreg});
for i=1:nreg-1
    pp(i)=ppmake(x(ridx(i,1):ridx(i,2)),y(ridx(i,1):ridx(i,2),:),...
        methidx,rn(i),sy2,psy2,clamp{i});
end

% combine pp
pp=ppcombine(pp);

% exit if pp
if(ppout); varargout{1}=pp; return; end

% get values
varargout{1}=ppdcval(pp,xi);
if(nargout>1); varargout{2}=ppdcval(pp,xi,false); end

% undo extrap if wanted
if(~exist('extrap','var'))
    switch methidx
        case {3 4 5} % spline/pchip/cubic
            extrap='extrap';
        otherwise
            extrap=nan;
    end
end
if(~strcmpi(extrap,'extrap'))
    if(~isnumeric(extrap))
        error('seizmo:interpdc1:badEXTRAP',...
            'Invalid EXTRAP option!');
    elseif(~isscalar(extrap))
        error('seizmo:interpdc1:badEXTRAP',...
            'EXTRAP must be scalar!');
    end
    varargout{1}(xi<min(x) | xi>max(x),:)=extrap;
    if(nargout>1); varargout{2}(xi<min(x) | xi>max(x),:)=extrap; end
end

% reshape to match xi
varargout{1}=reshape(varargout{1},syi);
if(nargout>1); varargout{2}=reshape(varargout{2},syi); end

end

function [pp]=ppmake(x,y,m,ny,sy2,psy2,clamp)
% gets pp for this section

% check clamp
if(~isempty(clamp))
    sc=size(clamp);
    if(~isnumeric(clamp) || ~isequal(sc,[2 sy2]))
        error('seizmo:interpdc1:badCLAMP',...
            ['CLAMP must be numeric & have dimension of [2 ' ...
            num2str(sy2) ']!']);
    end
end

% much of this is ripped from interp1's ppinterp
switch m
    case 1 % nearest
        pp=mkpp([x(1); (x(1:end-1)+x(2:end))/2; x(end)].',y.',sy2);
    case 2 % linear
        page1=(diff(y)./repmat(diff(x),[1 psy2])).';
        page2=(reshape(y(1:end-1,:),[ny-1 psy2])).';
        pp=mkpp(x.',cat(3,page1,page2),sy2);
    case 3 % spline
        if(isempty(clamp))
            pp=spline(x.',reshape(y.',[sy2 ny]));
        else
            pp=spline(x.',...
                reshape([clamp(1,:); y; clamp(2,:)].',[sy2 ny+2]));
        end
        % spline will return short orders for short sections
        % - increase to order 4 to make it common
        if(pp.order~=4)
            pp.coefs=[zeros(size(pp.coefs,1),4-pp.order) pp.coefs];
            pp.order=4;
        end
    case {4 5} % pchip/cubic
        pp=pchip(x.',reshape(y.',[sy2 ny]));
    case 6 % v5cubic
        if(ny==2)
            % this avoids breakage for 2pt
            % - just linear interp
            page1=(diff(y)./repmat(diff(x),[1 psy2])).';
            page2=(reshape(y(1:end-1,:),[ny-1 psy2])).';
            pp=mkpp(x.',cat(3,page1,page2),sy2);
            % increase to order 4 to make it common
            if(pp.order~=4)
                pp.coefs=[zeros(size(pp.coefs,1),4-pp.order) pp.coefs];
                pp.order=4;
            end
        else
            b=diff(x);
            if norm(diff(b),Inf)<=eps(norm(x,Inf))
                % equally spaced
                a=repmat(b,[1 psy2]).';
                % this is broken for 2 point sections
                yReorg=[3*y(1,:)-3*y(2,:)+y(3,:); ...
                    y; ...
                    3*y(ny,:)-3*y(ny-1,:)+y(ny-2,:)];
                y1=yReorg(1:end-3,:).';
                y2=yReorg(2:end-2,:).';
                y3=yReorg(3:end-1,:).';
                y4=yReorg(4:end,:).';
                page1=(-y1+3*y2-3*y3+y4)./(2*a.^3);
                page2=(2*y1-5*y2+4*y3-y4)./(2*a.^2);
                page3=(-y1+y3)./(2*a);
                page4=y2;
                coefs=cat(3,page1,page2,page3,page4);
                pp=mkpp(x.',coefs,sy2);
            else
                % not equally spaced
                pp=spline(x.',reshape(y.',[sy2 ny]));
                % spline will return short orders for short sections
                % - increase to order 4 to make it common
                if(pp.order~=4)
                    pp.coefs=[zeros(size(pp.coefs,1),4-pp.order) pp.coefs];
                    pp.order=4;
                end
            end
        end
end
pp.orient='first';
end

function [ppnew]=ppcombine(pp)
% combine pp(s) from each section into 1 pp for entire dataset

% no discontinuities, exit quickly
if(isscalar(pp)); ppnew=pp; return; end

% reorder if regions were decreasing
if(pp(1).breaks(1)>pp(end).breaks(1)); pp=pp(end:-1:1); end

% make new pp struct from pieces
% - this assumes sections are monotonically increasing
% - also assumes arrays are concatenatable (see spline hack in ppmake)
ppnew.form=pp(1).form;
ppnew.breaks=unique([pp.breaks]);
ppnew.coefs=cell2mat({pp.coefs}.');
ppnew.pieces=sum([pp.pieces]);
ppnew.order=max([pp.order]);
ppnew.dim=pp(1).dim;
ppnew.orient=pp(1).orient;

end
