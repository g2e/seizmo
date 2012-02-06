function [varargout]=plotspectra2(varargin)
%PLOTSPECTRA2    Overlay plot of SEIZMO data record spectra
%
%    Usage:    plotspectra2(data)
%              plotspectra2(data1,data2) or plotspectra2(data1,...,dataN)
%              plotspectra2(...,'cmp',spectral_cmp,...)
%              plotspectra2(...,'option',value,...)
%              ax=plotspectra2(...)
%
%    Description:
%     PLOTSPECTRA2(DATA) draws the amplitude spectra of all non-xyz records
%     in SEIZMO struct DATA over one another in the same plot with the HSV
%     colormap determining the line colors.  The plot is drawn in a new
%     figure.
%
%     PLOTSPECTRA2(DATA1,DATA2) OR PLOTSPECTRA2(DATA1,...,DATAN) will draw
%     records with the same index in their dataset in the same subplot.  So
%     DATA1(3) and DATA2(3) will be plotted together in the 3rd subplot.
%     All datasets must have the same number of records or be scalar.  This
%     is extremely useful for plotting real vs synthetic data.
%
%     PLOTSPECTRA2(...,'CMP',SPECTRAL_CMP,...) changes the spectral
%     component being drawn to CMP.  CMP may be 'am', 'ph', 'rl', 'im', or
%     'pw'.
%
%     PLOTSPECTRA2(...,'OPTION',VALUE,...) sets certain plotting options to
%     do simple manipulation of the plots.  Available options are:
%      FGCOLOR    -- foreground color (axes, text, labels)
%      BGCOLOR    -- background color (does not set figure color)
%      AXIS       -- axes to plot in (must be multiple if multidataset)
%      COLORMAP   -- colormap for coloring data
%      XLABEL     -- x axis label
%      YLABEL     -- y axis label
%      TITLE      -- title
%      XLIM       -- x axis limits (tight by default)
%      YLIM       -- y axis limits (tight by default)
%      LINEWIDTH  -- line width of records (default is 1)
%      LINESTYLE  -- line style of records (can be char/cellstr array)
%      NUMCOLS    -- number of subplot columns
%      UTC        -- plot in absolute time if TRUE (UTC, no leap support)
%      DATEFORMAT -- date format used if ABSOLUTE (default is auto)
%      XDIR       -- 'normal' or 'reverse'
%      YDIR       -- 'normal' or 'reverse'
%      FONTSIZE   -- size of fonts in the axes
%      ALIGN      -- ignore label and tick overlaps when aligning subplots
%      XSCALE     -- 'linear' or 'log'
%      YSCALE     -- 'linear' or 'log'
%      CMP        -- spectral component string
%      UNWRAP     -- unwrap phase before plotting
%
%     AX=PLOTSPECTRA2(...) returns the handle for the axis drawn in.  This
%     is useful for more detailed plot manipulation.
%
%    Notes:
%
%    Examples:
%     % overlay the first 4 records' spectra:
%     plotspectra2(data(1:4))
%
%     % plot first 4 records' spectra against the next 4
%     plot2(data(1:4),data(5:8))
%
%    See also: PLOTSPECTRA0, PLOTSPECTRA1, SPECTRASECTION, PLOT0, PLOT1,
%              PLOT2, RECORDSECTION

%     Version History:
%        Aug. 15, 2010 - initial version
%        Dec. 21, 2011 - add power option
%        Feb.  6, 2012 - fixed error bug for cmp=pw, pw is now in dB so
%                        yscale is linear, better getheader usage
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  6, 2012 at 23:00 GMT

% todo:
% - autoticks in log scale is not useful (Matlab bug)

% check nargin
error(nargchk(1,inf,nargin));

% find seizmo structures
issz=false(nargin,1);
for i=1:nargin
    issz(i)=isseizmo(varargin{i},'dep');
end
nd=sum(issz);
if(~nd)
    error('seizmo:plot2:badInput',...
        'No datasets provided!');
end

% separate datasets from option/value pairs
data=varargin(issz);
varargin(issz)=[];

% require all datasets to have equal number of records
if(~isequalnumelorscalar(data{:}))
    error('seizmo:plot2:badInput',...
        'Datasets must have the same number of records!');
end

% extract options (this requires that we actually parse twice...oh well)
opt=parse_seizmo_plot_options(varargin{:});

% loop over datasets
for i=1:nd
    % get datatype
    iftype=getheader(data{i},'iftype id');
    time=strcmpi(iftype,'itime') | strcmpi(iftype,'ixy');
    spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
    
    % convert timeseries to spectra
    if(sum(time)); data{i}(time)=dft(data{i}(time)); end
    
    % extract component of choice
    yscale='linear';
    switch lower(opt.SPECTRALCMP)
        case {'a' 'am' 'amp' 'amplitude'}
            data{i}(time | spec)=keepam(data{i}(time | spec));
            yscale='log';
            spylabel='amp';
        case {'p' 'ph' 'phase'}
            data{i}(time | spec)=keepph(data{i}(time | spec));
            spylabel='phase';
            
            % unwrap phase if desired
            if(opt.UNWRAP)
                data{i}(time | spec)=unwrapphase(data{i}(time | spec));
            end
        case {'r' 'rl' 'real'}
            data{i}(time | spec)=keeprl(data{i}(time | spec));
            spylabel='real cmp';
        case {'i' 'im' 'imag' 'imaginary'}
            data{i}(time | spec)=keepim(data{i}(time | spec));
            spylabel='imag cmp';
        case {'pw' 'pow' 'power'}
            data{i}(time | spec)=keeppw(data{i}(time | spec));
            spylabel='power (dB)';
        otherwise
            error('seizmo:plotspectra2:badInput',...
                'Unknown spectral component: %s',opt.SPECTRALCMP);
    end
end

% pass to plot2 with a better options (can be overrode)
opt.AXIS=plot2(data{:},'xlabel','Freq (Hz)','ylabel',spylabel,...
    'xscale','log','yscale',yscale,varargin{:});

% output axes if wanted
if(nargout); varargout{1}=opt.AXIS; end

end
