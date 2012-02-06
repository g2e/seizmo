function [varargout]=spectrasection(data,varargin)
%SPECTRASECTION    Plots SEIZMO data record spectra in a record section
%
%    Usage:    spectrasection(data)
%              spectrasection(...,'cmp',spectral_cmp,...)
%              spectrasection(...,'option',value,...)
%              ax=spectrasection(...)
%
%    Description:
%     SPECTRASECTION(DATA) draws the amplitude spectra of all non-xyz
%     records in SEIZMO struct DATA spaced out by their 'gcarc' header
%     field values (degree distance from event location given by the EVLA &
%     EVLO header fields).  The records are normalized as a group with a
%     maximum amplitude range corresponding to a third of the y axis range.
%     The spectra are scaled as if in logarithmic space and are colored as
%     distinct colors in the HSV colormap.
%
%     SPECTRASECTION(...,'CMP',SPECTRAL_CMP,...) changes the spectral
%     component being drawn to CMP.  CMP may be 'am', 'ph', 'rl', 'im', or
%     'pw'.
%
%     SPECTRASECTION(...,'OPTION',VALUE,...) sets certain plotting options
%     to do simple manipulation of the plots.  Available options are:
%      FGCOLOR      -- foreground color (axes, text, labels)
%      BGCOLOR      -- background color (does not set figure color)
%      AXIS         -- axes to plot in
%      COLORMAP     -- colormap for coloring data
%      XLABEL       -- x axis label
%      YLABEL       -- y axis label
%      TITLE        -- title
%      XLIM         -- x axis limits (tight by default)
%      YLIM         -- y axis limits (tight by default)
%      LINEWIDTH    -- line width of records (default is 1)
%      LINESTYLE    -- line style of records (can be char/cellstr array)
%      NUMCOLS      -- number of subplot columns
%      UTC          -- plot in absolute time if TRUE (UTC, no leap support)
%      DATEFORMAT   -- date format used if ABSOLUTE (default is auto)
%      NORMSTYLE    -- normalize 'individually' or as a 'group'
%      NORMMAX      -- max value of normalized records
%      NORM2YAXIS   -- scale to yaxis range (NORMMAX is fraction of range)
%      XDIR         -- 'normal' or 'reverse'
%      YDIR         -- 'normal' or 'reverse'
%      FONTSIZE     -- size of fonts in the axes
%      YFIELD       -- header field for y-axis positioning of records
%      XSCALE       -- 'linear' or 'log'
%      YSCALE       -- 'linear' or 'log'
%      AMPSCALE     -- 'linear' or 'log'
%      CMP          -- spectral component string
%      UNWRAP       -- unwrap phase before plotting
%
%     AX=SPECTRASECTION(...) returns the handle for the axis drawn in.
%     This is useful for more detailed plot manipulation.
%
%    Notes:
%
%    Examples:
%     % see how the spectra morphs with distance
%     spectrasection(data)
%
%     % now with azimuth
%     spectrasection(data,'yfield','az')
%
%    See also: PLOTSPECTRA0, PLOTSPECTRA1, PLOTSPECTRA2, PLOT0, PLOT1,
%              PLOT2, RECORDSECTION

%     Version History:
%        Aug. 15, 2010 - initial version
%        Dec. 21, 2011 - add power option
%        Feb.  6, 2012 - better getheader usage, ampscale=linear for cmp=pw
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  6, 2012 at 23:00 GMT

% todo:
% - autoticks in log scale is not useful (Matlab bug)

% check nargin
error(nargchk(1,inf,nargin));

% check struct
error(seizmocheck(data,'dep'));

% extract options (this requires that we actually parse twice...oh well)
opt=parse_seizmo_plot_options(varargin{:});

% get datatype
iftype=getheader(data,'iftype id');
time=strcmpi(iftype,'itime') | strcmpi(iftype,'ixy');
spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');

% convert timeseries to spectra
if(sum(time)); data(time)=dft(data(time)); end

% extract component of choice
ampscale='linear';
switch lower(opt.SPECTRALCMP)
    case {'a' 'am' 'amp' 'amplitude'}
        data(time | spec)=keepam(data(time | spec));
        ampscale='log';
    case {'p' 'ph' 'phase'}
        data(time | spec)=keepph(data(time | spec));
        
        % unwrap phase if desired
        if(opt.UNWRAP)
            data(time | spec)=unwrapphase(data(time | spec));
        end
    case {'r' 'rl' 'real'}
        data(time | spec)=keeprl(data(time | spec));
    case {'i' 'im' 'imag' 'imaginary'}
        data(time | spec)=keepim(data(time | spec));
    case {'pw' 'pow' 'power'}
        data(time | spec)=keeppw(data(time | spec));
    otherwise
        error('seizmo:spectrasection:badInput',...
            'Unknown spectral component: %s',opt.SPECTRALCMP);
end

% pass to recordsection with a better options (can be overrode)
opt.AXIS=recordsection(data,'xlabel','Freq (Hz)',...
    'xscale','log','ampscale',ampscale,varargin{:});

% output axes if wanted
if(nargout); varargout{1}=opt.AXIS; end

end
