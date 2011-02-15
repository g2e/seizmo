function [varargout]=plot_cmb_measurements(pf,field,varargin)
%PLOT_CMB_MEASUREMENTS    Plots core-diffracted analysis measurements
%
%    Usage:    plot_cmb_measurements(pf,field)
%              plot_cmb_measurements(pf,field,'prop1',val1,...)
%
%    Description:
%     PLOT_CMB_MEASUREMENTS(PF,FIELD) plots the measurements given by FIELD
%     stored in cmb profile struct PF with errorbars in both the x
%     (frequency) and y (slowness or decay constant) directions.  Filter
%     errorbars are based on the corners while measurement errorbars are 1
%     standard deviation.
%
%     PLOT_CMB_MEASUREMENTS(PF,FIELD,'PROP1',VAL1,...) passes
%     property/value pairs on to PLOTERR.  May NOT lead with a LINESPEC
%     string.  See PLOTERR for details.
%
%    Notes:
%
%    Examples:
%     % Plot corrected vs uncorrected of both slowness & decay constant:
%     fh=figure;
%     ax=makesubplots(2,1,1:2,'align','parent',fh);
%     plot_cmb_measurements(pf,'slow','parent',ax(1));
%     plot_cmb_measurements(pf,'cslow','parent',ax(1));
%     plot_cmb_measurements(pf,'decay','parent',ax(2));
%     plot_cmb_measurements(pf,'cdecay','parent',ax(2));
%
%    See also: CMB_PDF_MTX, MAP_CMB_PROFILES, SLOWDECAYPROFILES,
%              SLOWDECAYPAIRS, PLOTERR

%     Version History:
%        Feb.  1, 2011 - initial version
%        Feb. 10, 2011 - reusing axes works, better labeling, doc update
%        Feb. 12, 2011 - altered line style code
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 12, 2011 at 13:35 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% check profile struct
reqfields={'gcdist','azwidth','slow','slowerr','decay','decayerr',...
    'cslow','cslowerr','cdecay','cdecayerr','cluster','kname','st','ev',...
    'delaz','synthetics','earthmodel','corrections','corrcoef','freq',...
    'phase','runname','dirname','time'};
if(~isstruct(pf) || any(~isfield(pf,reqfields)))
    error('seizmo:cmb_pdf_mtx:badInput',...
        ['PF must be a struct with the fields:\n' ...
        sprintf('''%s'' ',reqfields{:}) '!']);
end

% check field
if(~ischar(field) || ...
        ~any(strcmpi(field,{'slow' 'cslow' 'decay' 'cdecay'})))
    error('seizmo:cmb_pdf_mtx:badInput',...
        'FIELD must be ''SLOW'', ''CSLOW'', ''DECAY'', or ''CDECAY''!');
end
field=lower(field);

% get periods
p0=1./cell2mat({pf.freq}');
lp=p0(:,2); up=p0(:,1); p0=mean(p0,2);

% get values & errors
x0=[pf.(field)].';
e0=[pf.([field 'err'])].';

% get earthmodel
emod={pf.earthmodel}';
[umod,idx,idx]=unique(emod);
nmod=max(idx);

% new plot or use existing
ax=axescheck(varargin{:});
if(isempty(ax))
    % new plot
    fh=figure();
    ax=axes('parent',fh);
    oldpts=0;
else
    % bring to front
    axes(ax);
    oldpts=numel(findobj(ax,'tag','points'));
end

% get/set hold state
held=ishold(ax);
hold(ax,'on');

% loop over models
h=nan(nmod,1); hyerr=h; hxerr=h;
for i=1:nmod
    % measurements from this model
    a=find(idx==i);
    
    % get random linespec (is a subfunction)
    [l,m,c,mc]=randlinespec(5+i+oldpts);
    
    % plot measurements
    [h(i),hyerr(i),hxerr(i)]=ploterr(ax,p0(a),x0(a),{lp(a) up(a)},e0(a),...
        [l m],'color',c,'markerfacecolor',mc,'markersize',10,varargin{:});
end

% change names if corrected
switch field
    case {'cslow' 'cdecay'}
        datas=strcmpi('data',umod);
        umod(datas)={'CORRECTED'};
        umod(~datas)=strcat(umod(~datas),{' (3D)'});
end

% tag handles
set(h,'tag','points',{'displayname'},umod);
set(hyerr,'tag','yerrorbars',{'displayname'},umod);
set(hxerr,'tag','xerrorbars',{'displayname'},umod);

% move errorbars to rear (they can clutter)
movekids(findobj(ax,'tag','points'),'front');

% labeling (if none exists)
if(isempty(get(get(ax,'xlabel'),'string')))
    xlabel(ax,'Period (s)');
end
if(isempty(get(get(ax,'ylabel'),'string')))
    switch field
        case {'slow' 'cslow'}
            ylabel(ax,'Ray Parameter (s/^o)');
        case {'decay' 'cdecay'}
            ylabel(ax,'Decay Constant');
    end
end
if(isempty(get(get(ax,'title'),'string')))
    switch field
        case 'slow'
            title(ax,'Ray Parameter Dispersion');
        case 'cslow'
            title(ax,['Ray Parameter Dispersion Corrected ' ...
                'for 3D Heterogeniety']);
        case 'decay'
            title(ax,'Decay Constant Dispersion');
        case 'cdecay'
            title(ax,['Decay Constant Dispersion Corrected ' ...
                     'for Geometrical Spreading']);
    end
end

% legend (regrab points objects to get those from previous plotting)
drawnow; % legends can have issues if plot not drawn yet
lh=legend(findobj(ax,'tag','points'));
set(lh,'interpreter','none');

% reset hold state
if(~held); hold(ax,'off'); end

% output
if(nargout>1)
    varargout={h hyerr hxerr};
elseif(nargout)
    varargout{1}=[h; hyerr; hxerr];
end

end


function [line,marker,color,mcolor]=randlinespec(n)
%RANDLINESPEC    Returns randomly chosen line specifiers

% possibilities
% - must have a line
% - line is not "light colored" (marker can be)
% - color set is limited to a somewhat extended colorset
linestyles={'-' '--' '-.' ':'};
markers={'+' '*' '.' 'x' 'o' 's' 'd' '^' 'v' '>' '<' 'p' 'h'};
colors={'r' 'o' 'b' 'v' 'm' 'p' 'k'};
mcolors={'r' 'o' 'y' 'l' 'g' 'a' 'c' 's' 'b' 'v' 'm' 'p' 'k'};

% random linespec or not
if(nargin)
    % pre-established but few repeats
    line=linestyles{mod(n-1,4)+1};
    tmp=mod(22-n,22)+1;
    if(tmp>13); tmp=tmp-9; end
    marker=markers{tmp};
    color=name2rgb(colors{mod(n-1,7)+1});
    if(mod(n,2))
        mcolor=name2rgb(mcolors{mod(9-n,13)+1});
    else
        mcolor=name2rgb(mcolors{mod(13-n,13)+1});
    end
else % random
    r=randperm(4);
    line=linestyles{r(1)};
    r=randperm(13);
    marker=markers{r(1)};
    r=randperm(7);
    color=name2rgb(colors{r(1)});
    r=randperm(13);
    mcolor=name2rgb(mcolors{r(1)});
end

end

