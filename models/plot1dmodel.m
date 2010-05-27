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
%    Description: PLOT1DMODEL(MODEL) makes separate subplots for properties
%     in the 1D models of MODEL for quick inspection & comparison.  The
%     depths extend from the surface (0km) to the center of the Earth
%     (6371km).
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
%     axes that the plotting should be done in.  The axis will be replaced.
%
%     AX=PLOT1DMODEL(...) outputs the handles to all the axes.
%
%    Notes:
%
%    Examples:
%     Perturb PREM and compare:
%      newmod=perturb_1dmodel(prem,'prem+ddp+ulvz',...
%                             'vs',[2600    0 1 1 10;
%                                   2700   -1 4 0  0;
%                                   2700    2 4 2 10;
%                                   2891 -0.2 1 0  0],...
%                             'vs',[2891-25 -10 4 1 5;
%                                   2891    -10 4 0 0]);
%      h=plot1dmodel([prem newmod],[],[2500 3000]);
%      set(findobj(h,'tag','vs'),'xlim',[6 8]); % resetting due to core
%
%    See also: PREM, AK135, IASP91, PERTURB_1DMODEL, AVAILABLE_1DMODELS,
%              CHK1DMODEL, FLATTEN_1DMODEL

%     Version History:
%        May  27, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  27, 2010 at 17:00 GMT

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
if(nargin<4 || isempty(cmap)); cmap=@hsvspin; end

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

% get line colors
if(isa(cmap,'function_handle'))
    cmap=cmap(nmod);
    if(isreal(cmap) && ~isequal(size(cmap),[nmod 3]))
        error('seizmo:plot1dmodel:badCMAP',...
            'CMAP must be a Nx3 RGB array or a handle to a colormap!');
    end
end

% default colors
if(nargin<5); fgcolor='w'; bgcolor='k'; end
if(nargin<6)
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
end
if(nargin<7)
    if(isempty(fgcolor))
        if(isempty(bgcolor))
            fgcolor='w'; bgcolor='k';
        else
            fgcolor=invertcolor(bgcolor,true);
        end
    elseif(isempty(bgcolor))
        if(isempty(fgcolor))
            fgcolor='w'; bgcolor='k';
        else
            bgcolor=invertcolor(fgcolor,true);
        end
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
if(nargin<7 || isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    figure('color',bgcolor,'name','1D EARTH MODEL PROPERTIES');
    ax=gca;
else
    axes(ax);
end

% get axis perimeter & figure, then delete
op=get(ax,'outerposition');
fh=get(ax,'parent');
delete(ax);

% get outer position of each subplot
tpl=min([0.0875 0.2*op(3) 0.2*op(4)]);
jt=2*tpl;
width=(op(3)-tpl)/nf;
height=op(4);
left=[op(1) op(1)+width+tpl:width:op(1)+op(3)-width];
bottom=op(2)*ones(nf,1);
ops=[left(:) bottom(:) ...
    [width+tpl; width*ones(nf-1,1)] height*ones(nf,1)];

% get inner position of each subplot
ipad=0.1;
ips=[ops(:,1)+ipad*width+[tpl; zeros(nf-1,1)] ops(:,2)+jt ...
    width*ones(nf,1)-2*ipad*width ops(:,4)-2*jt];

% get logical array of depths in range
didx=cell(nmod,1);
for i=1:nmod
    didx{i}=model(i).depth>=drng(1) | model(i).depth<=drng(2);
end

% each property is another column subplot
% only far left shows depth values
% make sure depths are synced
h=nan(nf,1);
for i=1:nf
    % make plot
    figure(fh);
    h(i)=axes('outerposition',ops(i,:),'position',ips(i,:));
    plot(model(1).(fields{i})(didx{1}),model(1).depth(didx{1}),...
        'color',cmap(1,:));
    set(h(i),'ydir','reverse','tag',fields{i});
    hold on
    for j=2:nmod
        plot(model(j).(fields{i})(didx{j}),model(j).depth(didx{j}),...
            'color',cmap(j,:));
    end
    hold off
    ylim(drng);
    
    % color & label (or delabel)
    fncyttl=texlabel(fields{i});
    fncyttl(2)=upper(fncyttl(2));
    title(fncyttl);
    set(h(i),'color',bgcolor,'xcolor',fgcolor,'ycolor',fgcolor);
    if(i==1)
        ylabel('Depth (km)');
    else
        % only the first plot has depths labeled
        set(h(i),'yticklabel',[]);
    end
    if(i==nf)
        legend(names,'textcolor',fgcolor,'edgecolor',fgcolor);
    end
    set(get(h(i),'YLabel'),'color',fgcolor);
    set(get(h(i),'title'),'color',fgcolor);
end

% output
if(nargout); varargout{1}=h; end

end
