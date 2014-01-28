function [opt]=parse_seizmo_plot_options(varargin)
%PARSE_SEIZMO_PLOT_OPTIONS    Unified SEIZMO plotting option parser

% Nov.  8, 2011 - added sft support
% Nov. 10, 2011 - fixed fg/bg fix code, sft defaults now setup
% Feb. 24, 2012 - pass char fg/bg inputs to name2rgb
% May  29, 2012 - pow2pad=0 by default
% Aug.  2, 2012 - added ftan support
% Jan. 30, 2013 - more distcut aliases
% Feb. 27, 2013 - parent allowed for specifying axis
% Jan. 27, 2014 - fix for octave handling of function handles

% todo:
% - value checking
% - documentation
% - support more plotting options and functions
% - per cmp linestyle, color, linewidth

% put function-dependent defaults at front of options list
pf=star69; % who is the calling function?
switch pf
    case {'drawmarkers' 'dm'}
        varargin=[{'m' true 'mac' 'g' 'mfc' 'r' 'moc' [1 .5 0] ...
            'mtc' 'y' 'mha' 'left' 'mva' 'top' 'mlw' 1 'mfs' 6 ...
            'mfw' 'normal' 'mfn' 'helvetica' 'mh' 50 'fm' 100 ...
            'fg' [] 'bg' []} varargin];
    case {'p0' 'plot0'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' 1 ...
            'ylabel' 1 'title' 1 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'abs' false 'dateformat' [] ...
            'namesonyaxis' false  'fontsize' 10  'markers' false ...
            'normstyle' 'group' 'normmax' 1/3 'norm2yaxis' true ...
            'xdir' 'normal' 'ydir' 'normal' 'xscale' 'linear' ...
            'yscale' 'linear' 'ampscale' 'linear' 'fontweight' 'bold'} ...
            varargin];
    case {'p1' 'plot1'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' 1 ...
            'ylabel' 1 'title' 1 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'ncols' [] 'abs' false 'dateformat' [] ...
            'xdir' 'normal' 'ydir' 'normal' 'fontsize' 10 'align' false ...
            'xscale' 'linear' 'yscale' 'linear' 'markers' false ...
            'fontweight' 'bold'} varargin];
    case {'p2' 'plot2'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' 1 ...
            'ylabel' 1 'title' 1 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'ncols' [] 'abs' false 'dateformat' [] ...
            'xdir' 'normal' 'ydir' 'normal' 'fontsize' 10 'align' false ...
            'xscale' 'linear' 'yscale' 'linear' 'fontweight' 'bold'} ...
            varargin];
    case {'plotclusters'}
        varargin=[{'fg' [] 'bg' [] 'poprng' [1 inf] 'clusters' [] ...
            'ax' [] 'ncols' [] 'align' false 'fontweight' 'bold' ...
            'fontsize' 10} varargin];
    case {'plotdendro'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' 1 ...
            'ylabel' 1 'title' 1 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'abs' false 'dateformat' [] 'markers' false ...
            'namesonyaxis' false  'fontsize' 10  'fontweight' 'bold' ...
            'normstyle' 'group' 'normmax' 1/3 'norm2yaxis' true ...
            'xdir' 'normal' 'ydir' 'normal' 'xscale' 'linear' ...
            'yscale' 'linear' 'ampscale' 'linear' 'cutoff' 0.2 ...
            'othercolor' [.5 .5 .5] 'cutoffcolor' 'r'} varargin];
    case {'pspec0' 'plotspectra0'}
        varargin=[{'fg' [] 'bg' [] 'unwrap' false 'cmp' 'a'} varargin];
    case {'pspec1' 'plotspectra1'}
        varargin=[{'fg' [] 'bg' [] 'unwrap' false 'cmp' 'a'} varargin];
    case {'pspec2' 'plotspectra2'}
        varargin=[{'fg' [] 'bg' [] 'unwrap' false 'cmp' 'a'} varargin];
    case {'recsec' 'recordsection'}
        varargin=[{'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' 'xlabel' 1 ...
            'ylabel' 1 'title' 1 'xlim' [] 'ylim' [] 'linewidth' 1 ...
            'linestyle' '-' 'abs' false 'dateformat' [] ...
            'namesonyaxis' false  'fontsize' 10 'yfield' 'gcarc' ...
            'normstyle' 'group' 'normmax' 1/3 'norm2yaxis' true ...
            'xdir' 'normal' 'ydir' 'normal' 'xscale' 'linear' ...
            'yscale' 'linear' 'ampscale' 'linear' 'markers' false ...
            'fontweight' 'bold'} varargin];
    case {'specsec' 'spectrasection'}
        varargin=[{'fg' [] 'bg' [] 'unwrap' false 'cmp' 'a'} varargin];
    case {'sft'}
        varargin=[{'u' '%' 'w' 2.5 'o' 75 'p2p' 0 'dbc' 'fire' 'hr' 0 ...
            'dbr' [] 'fr' [] 'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' ...
            'xlabel' 1 'ylabel' 1 'title' 1 'xlim' [] 'ylim' [] ...
            'linewidth' 1 'linestyle' '-' 'ncols' [] 'abs' false ...
            'dateformat' [] 'xdir' 'normal' 'ydir' 'normal' ...
            'fdir' 'normal' ...
            'fontsize' 10 'markers' false 'fontweight' 'bold'} varargin];
    case {'ftan'}
        varargin=[{'pf',@abs} {'nf' 100 'zc' 'fire' 'hr' 0 ...
            'zr' [] 'fr' [] 'fg' [] 'bg' [] 'ax' [] 'cmap' 'hsv' ...
            'xlabel' 1 'ylabel' 1 'title' 1 'xlim' [] 'ylim' [] ...
            'linewidth' 1 'linestyle' '-' 'ncols' [] 'abs' false ...
            'dateformat' [] 'xdir' 'normal' 'ydir' 'normal' ...
            'fdir' 'normal' ...
            'fontsize' 10 'markers' false 'fontweight' 'bold'} varargin];
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
        case {'axis' 'ax' 'axes' 'parent'}
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
        case {'namesonyaxis' 'nameonyaxis' 'yname' 'ynames'}
            opt.NAMESONYAXIS=val;
        case {'ydir'}
            opt.YDIR=val;
        case {'xdir'}
            opt.XDIR=val;
        case {'fontsize'}
            opt.FONTSIZE=val;
        case {'fontweight'}
            opt.FONTWEIGHT=val;
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
        case {'clucut' 'cutoff' 'distcut' 'distoff' 'dcut'}
            opt.CUTOFF=val;
        case {'clucutcolor' 'distcutcolor' 'cutoffcolor'}
            opt.CUTOFFCOLOR=val;
        case {'othercolor'}
            opt.OTHERCOLOR=val;
        case {'markers' 'showmarkers' 'drawmarkers' 'm'}
            opt.MARKERS=val;
        case {'markeracolor' 'mac'}
            opt.MARKERACOLOR=val;
        case {'markerfcolor' 'mfc'}
            opt.MARKERFCOLOR=val;
        case {'markerocolor' 'moc'}
            opt.MARKEROCOLOR=val;
        case {'markertcolor' 'mtc'}
            opt.MARKERTCOLOR=val;
        case {'markerhorzalign' 'mha'}
            opt.MARKERHORZALIGN=val;
        case {'markervertalign' 'mva'}
            opt.MARKERVERTALIGN=val;
        case {'markerlinewidth' 'mlw'}
            opt.MARKERLINEWIDTH=val;
        case {'markerfontsize' 'mfs'}
            opt.MARKERFONTSIZE=val;
        case {'markerfontweight' 'mfw'}
            opt.MARKERFONTWEIGHT=val;
        case {'markerfontname' 'mfn'}
            opt.MARKERFONTNAME=val;
        case {'markerheight' 'mh'}
            opt.MARKERHEIGHT=val;
        case {'flagmast' 'mastflags' 'fm' 'mf'}
            opt.FLAGMAST=val;
        case {'clusters'}
            opt.CLUSTERS=val;
        case {'poprng'}
            opt.POPRNG=val;
        case {'u' 'unit' 'units'}
            opt.UNITS=val;
        case {'w' 'window' 'width' 'length' 'len'}
            opt.WINDOW=val;
        case {'o' 'over' 'olap' 'overlap'}
            opt.OVERLAP=val;
        case {'p2p' 'padpower' 'pow2pad'}
            opt.POW2PAD=val;
        case {'dbc' 'dbcmap' 'dbcolor' 'dbcolormap'}
            opt.DBCMAP=val;
        case {'hr' 'hotrod'}
            opt.HOTROD=val;
        case {'dbr' 'dbrange'}
            opt.DBRANGE=val;
        case {'fr' 'frange' 'freqrange'}
            opt.FRANGE=val;
        case {'fdir' 'freqdir' 'frdir'}
            opt.FDIR=val;
        case {'numfreq' 'nfreq' 'nf'}
            opt.NUMFREQ=val;
        case {'postfunc' 'func' 'pf'}
            opt.POSTFUNC=val;
        case {'zcolormap' 'zcmap' 'zc'}
            opt.ZCMAP=val;
        case {'zrange' 'zrng' 'zr' 'zlim' 'zl'}
            opt.ZRANGE=val;
        otherwise
            error('seizmo:parse_seizmo_plot_options:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% fix fg/bg
if(ischar(opt.FGCOLOR)); opt.FGCOLOR=name2rgb(opt.FGCOLOR); end
if(ischar(opt.BGCOLOR)); opt.BGCOLOR=name2rgb(opt.BGCOLOR); end
if(isempty(opt.FGCOLOR))
    if(isempty(opt.BGCOLOR))
        opt.FGCOLOR='w'; opt.BGCOLOR='k';
    else
        opt.FGCOLOR=invertcolor(opt.BGCOLOR,true);
    end
elseif(isempty(opt.BGCOLOR))
    opt.BGCOLOR=invertcolor(opt.FGCOLOR,true);
end

end
