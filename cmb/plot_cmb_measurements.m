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
%     property/value pairs on to PLOTERR.  May be lead with a LINESPEC
%     string.  See PLOTERR for details.
%
%    Notes:
%
%    Examples:
%     % Plot corrected vs uncorrected of both slowness & decay constant:
%     fh=figure;
%     ax=makesubplots(2,1,1:2,'align','parent',fh);
%     hold(ax(1),'on');
%     hold(ax(2),'on');
%     h(:,1)=plot_cmb_measurements(pf,'slow','sk',...
%         'markerfacecolor','b','markersize',5,'parent',ax(1));
%     h(:,2)=plot_cmb_measurements(pf,'cslow','ok',...
%         'markerfacecolor','g','markersize',5,'parent',ax(1));
%     h(:,3)=plot_cmb_measurements(pf,'decay','sk',...
%         'markerfacecolor','b','markersize',5,'parent',ax(2));
%     h(:,4)=plot_cmb_measurements(pf,'cdecay','ok',...
%         'markerfacecolor','g','markersize',5,'parent',ax(2));
%     movekids(h(2:3,1:2),'back')
%     movekids(h(2:3,3:4),'back')
%     drawnow;
%     legend(h(1,1:2),'Raw','Corrected');
%     legend(h(1,3:4),'Raw','Corrected');
%     xlabel(ax(2),'Freq (Hz)');
%     ylabel(ax(1),'Ray Parameter (s/^o)');
%     ylabel(ax(2),'Decay Constant');
%     title(ax(1),'Ray Parameter Dispersion');
%     title(ax(2),'Decay Constant Dispersion');
%     linkaxes(ax,'x');
%     hold(ax(1),'off');
%     hold(ax(2),'off');
%
%    See also: CMB_PDF_MTX, MAP_CMB_PROFILES, SLOWDECAYPROFILES,
%              SLOWDECAYPAIRS, PLOTERR

%     Version History:
%        Feb.  1, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  1, 2011 at 13:35 GMT

% todo:
% - better labeling
%   - plot each set of points from different sources separately
%     x need a list of linestyle/color/marker combinations
% - how to add to same plot
%   - need to be able to show raw vs corrected easily
%     - detect if axes is from plot_cmb_measurements
%       - tag axes
%     - extract legend info if is
%       - [1,2,3,4]=legend(ax) -- 3 & 4 are what we need
%       - append to it
%     - do not allow mixed slow/decay
%     - do allow mixed phase (may want to show decay together)
%       - need to differentiate these
%       - title loses phase info

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

% get frequencies
f0=cell2mat({pf.freq}');
lf=f0(:,1); uf=f0(:,2); f0=mean(f0,2);

% get values & errors
x0=[pf.(field)].';
e0=[pf.([field 'err'])].';

% get earthmodel
emod={pf.earthmodel}';
[umod,idx,idx]=unique(emod);
nmod=max(idx);

% loop over models
h=nan(nmod,1); hyerr=h; hxerr=h;
for i=1:max(idx)
    % get random linespec
    [l,m,c,mc]=randlinespec(i);
    
    % plot measurements
    [h(i),hyerr(i),hxerr(i)]=ploterr(f0,x0,{lf uf},e0,[l m],...
        'color',c,'markerfacecolor',mc,varargin{:});
end

% move errorbars to rear
movekids(h,'front');

% legend
lh=legend(h,umod);
set(lh,'interpreter','none');

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
    if(tmp>13); tmp=tmp-13; end
    marker=markers{tmp};
    color=name2rgb(colors{mod(n-1,7)+1});
    mcolor=name2rgb(mcolors{mod(13-n,13)+1});
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

