function [opt]=parse_seizmo_plot_options(varargin)
%PARSE_SEIZMO_PLOT_OPTIONS    Unified SEIZMO plotting option parser

% todo:
% - value checking
% - documentation
% - support more plotting options and functions

% defaults
pf=star69;
switch pf
    case {'p0' 'plot0'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' [] ...
            'ylabel' [] 'title' [] 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'abs' false 'dateformat' [] ...
            'namesonyaxis' false  'fontsize' 6 ...
            'normstyle' 'group' 'normmax' 1/3 'norm2yaxis' true ...
            'xdir' 'normal' 'ydir' 'reverse' 'xscale' 'linear' ...
            'yscale' 'linear' 'ampscale' 'linear'} varargin];
    case {'p1' 'plot1'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' [] ...
            'ylabel' [] 'title' [] 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'ncols' [] 'abs' false 'dateformat' [] ...
            'xdir' 'normal' 'ydir' 'normal' 'fontsize' 6 'align' false ...
            'xscale' 'linear' 'yscale' 'linear'} ...
            varargin];
    case {'p2' 'plot2'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' [] ...
            'ylabel' [] 'title' [] 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'ncols' [] 'abs' false 'dateformat' [] ...
            'xdir' 'normal' 'ydir' 'normal' 'fontsize' 6 'align' false ...
            'xscale' 'linear' 'yscale' 'linear'} ...
            varargin];
    case {'plotdendro'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' [] ...
            'ylabel' [] 'title' [] 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'abs' false 'dateformat' [] ...
            'namesonyaxis' false  'fontsize' 6 ...
            'normstyle' 'group' 'normmax' 1/3 'norm2yaxis' true ...
            'xdir' 'normal' 'ydir' 'reverse' 'xscale' 'linear' ...
            'yscale' 'linear' 'ampscale' 'linear' 'clustercutoff' 0.2 ...
            'othercolor' [.5 .5 .5] 'clustercutoffcolor' 'r'} varargin];
    case {'psp0' 'plotsp0' 'plotspectra0'}
        varargin=[{'fg' [] 'bg' [] 'unwrap' false 'cmp' 'a'} varargin];
    case {'psp1' 'plotsp1' 'plotspectra1'}
        varargin=[{'fg' [] 'bg' [] 'unwrap' false 'cmp' 'a'} varargin];
    case {'psp2' 'plotsp2' 'plotspectra2'}
        varargin=[{'fg' [] 'bg' [] 'unwrap' false 'cmp' 'a'} varargin];
    case {'recsec' 'recordsection'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' [] ...
            'ylabel' [] 'title' [] 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'abs' false 'dateformat' [] ...
            'namesonyaxis' false  'fontsize' 6 'yfield' 'gcarc' ...
            'normstyle' 'group' 'normmax' 1/3 'norm2yaxis' true ...
            'xdir' 'normal' 'ydir' 'normal' 'xscale' 'linear' ...
            'yscale' 'linear' 'ampscale' 'linear'} varargin];
    case {'spsec' 'spectrasection'}
        varargin=[{'fg' [] 'bg' [] 'unwrap' false 'cmp' 'a'} varargin];
    otherwise
        error('seizmo:parse_seizmo_plot_options:badInput',...
            'Unknown plotting function: %s !',pf);
end
np=numel(varargin);

% check nargin
if(mod(np,2))
    error('seizmo:parse_seizmo_plot_options:badInput',...
        'Unpaired Option/Value!');
end

% options must be strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:parse_seizmo_plot_options:badInput',...
        'OPTION must be specified with a string!');
end

% loop over pairs
for i=1:2:np
    val=varargin{i+1};
    switch lower(varargin{i})
        case {'fgcolor' 'fgc' 'fg' 'foreground'}
            opt.FGCOLOR=val;
        case {'bgcolor' 'bgc' 'bg' 'background'}
            opt.BGCOLOR=val;
        case {'axis' 'ax' 'axes'}
            opt.AXIS=val;
        case {'cmap' 'colormap'}
            opt.CMAP=val;
        case {'xlabel'}
            opt.XLABEL=val;
        case {'ylabel'}
            opt.YLABEL=val;
        case {'title'}
            opt.TITLE=val;
        case {'xlim' 'xlimit' 'xlimits'}
            opt.XLIM=val;
        case {'ylim' 'ylimit' 'ylimits'}
            opt.YLIM=val;
        case {'linewidth' 'recwidth'}
            opt.LINEWIDTH=val;
        case {'linestyle' 'recstyle'}
            opt.LINESTYLE=val;
        case {'numcols' 'numcol' 'ncols' 'ncol' 'columns'}
            opt.NUMCOLS=val;
        case {'utc' 'abs' 'absolute'}
            opt.ABSOLUTE=val;
        case {'dateformat'}
            opt.DATEFORMAT=val;
        case {'normstyle'}
            opt.NORMSTYLE=val;
        case {'normmax'}
            opt.NORMMAX=val;
        case {'norm2yaxis'}
            opt.NORM2YAXIS=val;
        case {'namesonyaxis'}
            opt.NAMESONYAXIS=val;
        case {'ydir'}
            opt.YDIR=val;
        case {'xdir'}
            opt.XDIR=val;
        case {'fontsize'}
            opt.FONTSIZE=val;
        case {'align'}
            if(~islogical(val) || ~isscalar(val) || ~val)
                opt.ALIGN={};
            else
                opt.ALIGN={'align'};
            end
        case {'yfield'}
            opt.YFIELD=val;
        case {'xscale'}
            opt.XSCALE=val;
        case {'yscale'}
            opt.YSCALE=val;
        case {'ampscale'}
            opt.AMPSCALE=val;
        case {'cmp' 'spcmp' 'spectralcmp'}
            opt.SPECTRALCMP=val;
        case {'unwrap' 'unwrapphase'}
            opt.UNWRAP=val;
        case {'clucut' 'clustercutoff'}
            opt.CLUSTERCUTOFF=val;
        case {'clucutcolor' 'clustercutoffcolor'}
            opt.CLUSTERCUTOFFCOLOR=val;
        case {'othercolor'}
            opt.OTHERCOLOR=val;
        otherwise
            error('seizmo:parse_seizmo_plot_options:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% fix fg/bg
if(isempty(opt.FGCOLOR))
    if(isempty(opt.BGCOLOR))
        opt.FGCOLOR='w'; opt.BGCOLOR='k';
    else
        opt.FGCOLOR=invertcolor(opt.BGCOLOR,true);
    end
elseif(isempty(opt.BGCOLOR))
    if(isempty(opt.FGCOLOR))
        opt.FGCOLOR='w'; opt.BGCOLOR='k';
    else
        opt.BGCOLOR=invertcolor(opt.FGCOLOR,true);
    end
end

end
