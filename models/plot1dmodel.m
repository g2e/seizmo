function [varargout]=plot1dmodel(model,fields,drng,cmap,fgcolor,bgcolor,ax)
%PLOT1DMODEL    Plots 1D model properties
%
%    Usage:    plot1dmodel(model)
%              plot1dmodel(model,fields)
%              plot1dmodel(model,fields,drng)
%              plot1dmodel(model,fields,drng,cmap)
%              plot1dmodel(model,fields,drng,cmap,fgcolor,bgcolor)
%              plot1dmodel(model,fields,drng,cmap,fgcolor,bgcolor,ax)
%              ax=plot1dmodel(...)
%
%    Description:
%     PLOT1DMODEL(MODEL) makes separate subplots for properties in the 1D
%     models of MODEL for quick inspection & comparison.  The depths extend
%     from the surface (0km) to the center of the Earth (6371km).
%
%     PLOT1DMODEL(MODEL,FIELDS) specifies the fields to plot.  FIELDS must
%     be a char or cellstr array of property fields in MODEL.  The default
%     is all of the depth dependent properties.
%
%     PLOT1DMODEL(MODEL,FIELDS,DRNG) alters the depth range plotted.  This
%     is useful for focusing on a certain region (like the D" region or the
%     transition zone).  The default is [0 6371].
%
%     PLOT1DMODEL(MODEL,FIELDS,DRNG,CMAP) allows setting the line coloring
%     in the plots.  CMAP may be a NMODx3 RGB triplet or a string/function
%     handle to a colormap function.  NMOD is the number of models in
%     MODEL.  The default is 'hsvspin'.
%
%     PLOT1DMODEL(MODEL,FIELDS,DRNG,CMAP,FGCOLOR,BGCOLOR) specifies the
%     foreground & background colors of the plots.  The default values are
%     FGCOLOR = 'w' & BGCOLOR = 'k' (white on black).  FGCOLOR/BGCOLOR may
%     be a color name or a 1x3 RGB triplet.
%
%     PLOT1DMODEL(MODEL,FIELDS,DRNG,CMAP,FGCOLOR,BGCOLOR,AX) indicates the
%     axes that the plotting should be done in.  There should be as many
%     axes as there are fields in FIELDS.  The axes will be replaced.
%
%     AX=PLOT1DMODEL(...) outputs the handles to all the axes.
%
%    Notes:
%
%    Examples:
%     % Perturb PREM and compare:
%     newmod=perturb_1dmodel(prem,'prem+ddp+ulvz',...
%                            'vs',[2600    0 1 1 10;
%                                  2700   -1 4 0  0;
%                                  2700    2 4 2 10;
%                                  2891 -0.2 1 0  0],...
%                            'vs',[2891-25 -10 4 1 5;
%                                  2891    -10 4 0 0]);
%     h=plot1dmodel([prem newmod],[],[2500 3000]);
%     set(findobj(h,'tag','vs'),'xlim',[6 8]); % resetting due to core
%
%    See also: PREM, AK135, IASP91, PREM_PERFECT, PREM2_PERFECT,
%              PERTURB_1DMODEL, AVAILABLE_1DMODELS, CHK1DMODEL,
%              FLATTEN_1DMODEL

%     Version History:
%        May  27, 2010 - initial version
%        May  29, 2010 - turn off tex interpretation of model names
%        Aug. 17, 2010 - fancy titles in LaTeX
%        Feb. 28, 2010 - fix default color code bug
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 15:50 GMT

% todo:

% check nargin
error(nargchk(1,7,nargin));

% check 1d model
error(chk1dmodel(model));

% number of models
nmod=numel(model);

% model names
names={model.name}';

% what fields are available
reqfields={'name' 'ocean' 'crust' 'isotropic' 'refperiod' 'flattened' ...
    'depth'};
ofields=setdiff(fieldnames(model),reqfields);

% default options
if(nargin<2 || isempty(fields)); fields=ofields; end
if(nargin<3 || isempty(drng)); drng=[0 6371]; end
if(nargin<4 || isempty(cmap)); cmap='hsv'; end

% check options
if(ischar(fields)); fields=cellstr(fields); end
if(~iscellstr(fields) || any(~ismember(fields,ofields)))
    error('seizmo:plot1dmodel:badFIELDS',...
        'FIELDS must be a char/cellstr array of fields in MODEL!');
elseif(~isreal(drng) || ~isequal(size(drng),[1 2]) ...
        || any(drng<0 | drng>6371))
    error('seizmo:plot1dmodel:badDEPRNG',...
        'DEPRNG must be a real valued 1x2 vector of depths in [0 6371]!');
elseif(~ischar(cmap) && ~isa(cmap,'function_handle') ...
        && (isreal(cmap) && ~isequal(size(cmap),[nmod 3])))
    error('seizmo:plot1dmodel:badCMAP',...
        'CMAP must be a Nx3 RGB array or a handle to a colormap!');
end
drng=sort(drng);
if(ischar(cmap)); cmap=str2func(cmap); end

% number of fields
nf=numel(fields);

% nrows/ncols
ncols=fix(5*sqrt(nf/5));
if(ncols>nf); ncols=nf; end
nrows=ceil(nf/ncols);

% get line colors
if(isa(cmap,'function_handle'))
    cmap=cmap(nmod);
    if(isreal(cmap) && ~isequal(size(cmap),[nmod 3]))
        error('seizmo:plot1dmodel:badCMAP',...
            'CMAP must be a Nx3 RGB array or a handle to a colormap!');
    end
end

% default colors
if(nargin<5)
    fgcolor='w'; bgcolor='k';
elseif(nargin<6)
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
else
    if(isempty(fgcolor))
        if(isempty(bgcolor))
            fgcolor='w'; bgcolor='k';
        else
            fgcolor=invertcolor(bgcolor,true);
        end
    elseif(isempty(bgcolor))
        bgcolor=invertcolor(fgcolor,true);
    end
end

% change char to rgb
if(ischar(fgcolor)); fgcolor=name2rgb(fgcolor); end
if(ischar(bgcolor)); bgcolor=name2rgb(bgcolor); end

% check colors
if(~isreal(fgcolor) || ~isequal(size(fgcolor),[1 3]) ...
        || any(fgcolor<0 | fgcolor>1))
    error('seizmo:plot1dmodel:badFGCOLOR',...
        'FGCOLOR must be a RGB triplet or color name!');
elseif(~isreal(bgcolor) || ~isequal(size(bgcolor),[1 3]) ...
        || any(bgcolor<0 | bgcolor>1))
    error('seizmo:plot1dmodel:badBGCOLOR',...
        'BGCOLOR must be a RGB triplet or color name!');
end

% check handle
if(nargin<7 || isempty(ax) || numel(ax)~=numel(fields) || ~isreal(ax) ...
        || any(~ishandle(ax)) || any(~strcmp('axes',get(ax,'type'))))
    fh=figure('color',bgcolor,'name','1D EARTH MODEL PROPERTIES');
    ax=makesubplots(nrows,ncols,1:nf,'parent',fh);
end

% get logical array of depths in range
didx=cell(nmod,1);
for i=1:nmod
    didx{i}=model(i).depth>=drng(1) | model(i).depth<=drng(2);
end

% field to title
fancyfield={'rho' 'vp' 'vs' 'vb' 'qk' 'qu' 'shear' 'bulk' 'youngs' ...
    'lambda' 'poisson' 'm' 'g' 'vpv' 'vph' 'vsv' 'vsh' 'eta' 'vbv' 'vbh'...
    'poissonv' 'poissonh' 'shearv' 'shearh' 'bulkv' 'bulkh' 'youngsv' ...
    'youngsh' 'lambdav' 'lambdah' 'drho_dr' 'dvp_dr' 'dvs_dr' 'dvb_dr' ...
    'deta_dr' 'dvpv_dr' 'dvsv_dr' 'dvbv_dr' 'dvph_dr' 'dvsh_dr' ...
    'dvbh_dr' 'd2rho_dr2' 'd2vp_dr2' 'd2vs_dr2' 'd2vb_dr2' 'd2vpv_dr2' ...
    'd2vph_dr2' 'd2vsv_dr2' 'd2vsh_dr2' 'd2vbv_dr2' 'd2vbh_dr2' ...
    'd2eta_dr2' 'radius'};
fancyxlabel={'g/cm3' 'km/s' 'km/s' 'km/s' '' '' ...
    'Pa' 'Pa' 'Pa' 'Pa' '' 'kg' 'm/s2' 'km/s' 'km/s' 'km/s' 'km/s' '' ...
    'km/s' 'km/s' '' '' 'Pa' 'Pa' 'Pa' 'Pa' 'Pa' 'Pa' 'Pa' 'Pa' ...
    'g/cm3/Re' 'km/s/Re' 'km/s/Re' 'km/s/Re' ...
    '1/Re' 'km/s/Re' 'km/s/Re' 'km/s/Re' 'km/s/Re' 'km/s/Re' ...
    'km/s/Re' 'g/cm3/Re' 'km/s/Re2' 'km/s/Re2' 'km/s/Re2' 'km/s/Re2' ...
    'km/s/Re2' 'km/s/Re2' 'km/s/Re2' 'km/s/Re2' 'km/s/Re2' ...
    '1/Re2' 'km'};
fancytitle={'\rho' '\alpha' '\beta' '\Phi' 'Q_\kappa' 'Q_\mu' ...
    '\mu' '\kappa' 'E' '\lambda' '\gamma' 'm' 'g' '\alpha_V' '\alpha_H' ...
    '\beta_V' '\beta_H' '\eta' '\Phi_V' '\Phi_H' '\gamma_V' '\gamma_H' ...
    '\mu_V' '\mu_H' '\kappa_V' '\kappa_H' 'E_V' 'E_H' '\lambda_V' ...
    '\lambda_H' '\displaystyle\frac{\delta\rho}{\delta{r}}' ...
    '\displaystyle\frac{\delta\alpha}{\delta{r}}' ...
    '\displaystyle\frac{\delta\beta}{\delta{r}}' ...
    '\displaystyle\frac{\delta\Phi}{\delta{r}}' ...
    '\displaystyle\frac{\delta\eta}{\delta{r}}' ...
    '\displaystyle\frac{\delta\alpha_V}{\delta{r}}' ...
    '\displaystyle\frac{\delta\beta_V}{\delta{r}}' ...
    '\displaystyle\frac{\delta\Phi_V}{\delta{r}}' ...
    '\displaystyle\frac{\delta\alpha_H}{\delta{r}}' ...
    '\displaystyle\frac{\delta\beta_H}{\delta{r}}' ...
    '\displaystyle\frac{\delta\Phi_H}{\delta{r}}' ...
    '\displaystyle\frac{\delta^{2}\rho}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\alpha}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\beta}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\Phi}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\alpha_V}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\alpha_H}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\beta_V}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\beta_H}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\Phi_V}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\Phi_H}{\delta{r}^2}' ...
    '\displaystyle\frac{\delta^{2}\eta}{\delta{r}^2}' 'r'};

% each property is another subplot
% only far left shows depth values
% make sure depths are synced (linked)
for i=1:nf
    % adjust current axis
    cla(ax(i),'reset');
    
    % make plot
    hold(ax(i),'on');
    for j=1:nmod
        plot(ax(i),model(j).(fields{i})(didx{j}),...
            model(j).depth(didx{j}),...
            'color',cmap(j,:));
    end
    hold(ax(i),'off');
    
    % clean up axes
    set(ax(i),'color',bgcolor,'xcolor',fgcolor,'ycolor',fgcolor);
    set(ax(i),'ydir','reverse','tag',fields{i});
    ylim(ax(i),drng);
    box(ax(i),'on');
    
    % titles
    if(ismember(fields{i},fancyfield))
        title(ax(i),...
            ['$' fancytitle{ismember(fancyfield,fields{i})} '$'],...
            'interpreter','latex','color',fgcolor,'fontsize',16);
    else
        title(fields{i});
    end
    
    % xlabels
    if(ismember(fields{i},fancyfield))
        xlabel(ax(i),...
            fancyxlabel{ismember(fancyfield,fields{i})},...
            'color',fgcolor);
    end
    
    % y tick labels
    if(mod(i-1,ncols))
        % only the first plot has depths labeled
        set(ax(i),'yticklabel',[]);
    end
    if(i==nf)
        lh=legend(ax(i),names,'textcolor',fgcolor,'edgecolor',fgcolor);
        set(lh,'interpreter','none')
    end
end

% link y axes
linkaxes(ax,'y');

% super ylabel
superylabel(ax,'Depth (km)','color',fgcolor,'fontsize',16);

% output
if(nargout); varargout{1}=ax; end

end
