function [varargout]=plot_cmb_pdf(pf,field,y,p,fgcolor,bgcolor,ax)
%PLOT_CMB_PDF    Plot probability space of slowness or decay constant
%
%    Usage:    plot_cmb_pdf(pf,field)
%              plot_cmb_pdf(pf,field,y,p)
%              plot_cmb_pdf(pf,field,y,p,fgcolor,bgcolor)
%              plot_cmb_pdf(pf,field,y,p,fgcolor,bgcolor,ax)
%              pdf=plot_cmb_pdf(...)
%              [pdf,y,p]=plot_cmb_pdf(...)
%
%    Description:
%     PLOT_CMB_PDF(PF,FIELD) plots a PDF (probability density function) of
%     the measurements indicated by FIELD stored in the profile struct PF.
%     PF is expected to be in the format created by SLOWDECAYPROFILES or
%     SLOWDECAYPAIRS.  FIELD is one of the following:
%      'slow', 'cslow', 'decay', 'cdecay'
%     where the 'c' in front means that the measurement has corrections
%     applied (for 3D heterogeneous structure of the Earth or geometrical
%     spreading).  See SLOWDECAYPROFILES & SLOWDECAYPAIRS for more details.
%     The PDF is automatically scaled using the input data ranges (some
%     padding is added to the y-range to account for the errors in the
%     data).  The PDF is 100x100 pixels.
%
%     PLOT_CMB_PDF(PF,FIELD,Y,P) allows specifying the measurement values
%     Y and periods P at which the PDF is estimated.  PDF values are
%     calculated assuming the error for the data is gaussian (the given
%     error is 1 standard deviation and the measurement itself is the
%     mean).  Y may be a scalar (sets the number of points in an automatic-
%     ally determined range), a 2-element vector given as [MIN MAX], or a
%     regularly spaced vector of values (this is so pixel width can be
%     established for probability determination).  P inputs have similar
%     constraints.
%
%     PLOT_CMB_PDF(PF,FIELD,Y,P,FGCOLOR,BGCOLOR) adjusts the foreground
%     color FGCOLOR and the background color BGCOLOR.  The default is white
%     ('w') for the foreground and black ('k') for the background.  Note
%     that if one is specified and the other is not, an opposing color is
%     found using INVERTCOLOR.  The probability color scale is also changed
%     so the noise clip is at BGCOLOR.
%
%     PLOT_CMB_PDF(PF,FIELD,Y,P,FGCOLOR,BGCOLOR,AX) plots in the axes given
%     by handle AX.
%
%     PDF=PLOT_CMB_PDF(...) outputs the PDF matrix.  PDF is NYxNP in size
%     where NY is the number of points in the measurement value space and
%     NP is the number of points in period.  Use the next form if you want
%     the actual Y & P grid corresponding to PDF.
%
%     [PDF,Y,P]=PLOT_CMB_PDF(...) also returns the locations of the PDF
%     grid as measurement value (Y) and period (P).  See the Examples
%     section below for how to plot this information
%
%    Notes:
%     - Can be quite inaccurate when standard deviations are subpixel.  
%
%    Examples:
%     % A guide for how to plot the optional outputs:
%     [pdf,y,p]=plot_cmb_pdf(pf,'cslow');
%     fh=figure('color','k'); ax=axes('parent',fh);
%     imagesc(p,y,pdf,'parent',ax);
%     colormap(ax,fire);
%     set(ax,'ydir','normal','xcolor','w','ycolor','w');
%     xlabel(ax,'Period (sec)');
%     ylabel(ax,'Ray Parameter (s/^o)')
%     title(ax,'Ray Parameter Dispersion w/ 3D Corrections');
%     cb=colorbar('peer',ax);
%     ylabel(cb,'% Probability');
%     grid(ax,'on');
%
%    See also: PLOT_CMB_MEASUREMENTS, MAP_CMB_PROFILES, SLOWDECAYPAIRS,
%              SLOWDECAYPROFILES, IMAGESC

%     Version History:
%        Feb.  1, 2011 - initial version
%        Feb. 10, 2011 - reduced plot code redundancy, use period not freq,
%                        improved probability estimation, works with 1 pf
%        Feb. 17, 2011 - aesthetic touches
%        Mar. 30, 2011 - improve title and documentations
%        Apr. 22, 2011 - fixed docs, input arg order, pdf output transposed
%        Mar.  5, 2012 - axes handle input
%        Oct. 11, 2012 - drop corrections field requirement, period min/max
%                        input bugfix
%        Oct. 16, 2012 - color options
%        Feb. 26, 2013 - bugfix: precision issue with y interval detection
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 26, 2013 at 13:35 GMT

% todo:

% check nargin
error(nargchk(2,7,nargin));

% check input
reqfields={'gcdist','azwidth','slow','slowerr','decay','decayerr',...
    'cslow','cslowerr','cdecay','cdecayerr','cluster','kname','st','ev',...
    'delaz','synthetics','earthmodel','corrcoef','freq',...
    'phase','runname','dirname','time'};
if(~isstruct(pf) || any(~isfield(pf,reqfields)))
    error('seizmo:plot_cmb_pdf:badInput',...
        ['PF must be a struct with the fields:\n' ...
        sprintf('''%s'' ',reqfields{:}) '!']);
end

% check field
if(~ischar(field) || ...
        ~any(strcmpi(field,{'slow' 'cslow' 'decay' 'cdecay'})))
    error('seizmo:plot_cmb_pdf:badInput',...
        'FIELD must be ''SLOW'', ''CSLOW'', ''DECAY'', or ''CDECAY''!');
end
field=lower(field);

% get periods
p0=1./cell2mat({pf.freq}');

% get values & errors
y0=[pf.(field)].';
e0=[pf.([field 'err'])].';

% get min/max
pmin=min(p0(:)); pmax=max(p0(:));
ymin=min(y0-median(e0)); ymax=max(y0+median(e0));
ypad=(ymax-ymin)/4;

% default period/value ranges
% 100 pixels each way, scaled to data
if(nargin<3 || isempty(y)); y=linspace(ymin-ypad,ymax+ypad,100); end
if(nargin<4 || isempty(p)); p=linspace(pmin,pmax,100); end

% check colors
if(nargin<5);
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

% change char to something rgb
if(ischar(fgcolor)); fgcolor=name2rgb(fgcolor); end
if(ischar(bgcolor)); bgcolor=name2rgb(bgcolor); end

% default axes
if(nargin<7); ax=[]; end

% special 2-element range
if(numel(y)==2); y=linspace(y(1),y(2),100); end
if(numel(p)==2); p=linspace(p(1),p(2),100); end

% special positive scalar integer input
if(isscalar(y) && y==fix(y) && y>0); y=linspace(ymin-ypad,ymax+ypad,y); end
if(isscalar(p) && p==fix(p) && p>0); p=linspace(pmin,pmax,p); end

% check period/value
% - value needs to be regularly spaced
spacing=y(2)-y(1);
if(~isreal(p) || ~isvector(p) || any(p<=0))
    error('seizmo:plot_cmb_pdf:badInput',...
        'P must be a vector of positive real-valued periods in seconds!');
elseif(~isreal(y) || ~isvector(y) || any(abs(diff(y)-spacing)/spacing>.01))
    error('seizmo:plot_cmb_pdf:badInput',...
        'Y must be a vector of evenly spaced reals!');
end

% preallocate output
mtx=zeros(numel(p),numel(y));

% loop over period
for i=1:numel(p)
    % find period ranges encompassing current period
    idx=find(p(i)>=p0(:,2) & p(i)<=p0(:,1));
    nidx=numel(idx);
    
    % skip if none
    if(~nidx); continue; end
    
    % loop over values
    % - get values at edge & center of pixel so we can do
    %   trapezoidal integration for each pixel's probability
    for j=1:nidx
        mtx(i,:)=mtx(i,:) ...
            +2*gaussian(y,y0(idx(j)),e0(idx(j))) ...
            +gaussian(y-spacing/2,y0(idx(j)),e0(idx(j))) ...
            +gaussian(y+spacing/2,y0(idx(j)),e0(idx(j)));
    end
    
    % finish average & integration
    mtx(i,:)=mtx(i,:)/nidx/4*spacing*100;
end

% transpose pdf
mtx=mtx';

% output depends on number out
if(nargout>1)
    varargout={mtx y p};
elseif(nargout)
    varargout={mtx};
else
    % plot output
    if(isempty(ax))
        fh=figure('color',bgcolor);
        ax=axes('parent',fh);
        if(strcmp(bgcolor,'w') || isequal(bgcolor,[1 1 1]))
            colormap(ax,flipud(fire));
        elseif(strcmp(bgcolor,'k') || isequal(bgcolor,[0 0 0]))
            colormap(ax,fire);
        else
            if(ischar(bgcolor))
                bgcolor=name2rgb(bgcolor);
            end
            hsv=rgb2hsv(bgcolor);
            colormap(ax,hsvcustom(hsv));
        end
    end
    imagesc(p,y,mtx,'parent',ax);
    set(ax,'ydir','normal','color',bgcolor,...
        'xcolor',fgcolor,'ycolor',fgcolor);
    set(ax,'xminortick','on','yminortick','on');
    set(ax,'xgrid','on','ygrid','on','linewidth',1);
    xlabel(ax,'Period (s)','color',fgcolor);
    cb=colorbar('peer',ax,'xcolor',fgcolor,'ycolor',fgcolor);
    ylabel(cb,'% Probability','color',fgcolor);
    set(cb,'linewidth',1);
    switch field
        case 'slow'
            ylabel(ax,'Ray Parameter (s/^o)','color',fgcolor);
            title(ax,'Ray Parameter Dispersion','color',fgcolor);
        case 'cslow'
            ylabel(ax,'Ray Parameter (s/^o)','color',fgcolor);
            title(ax,['Ray Parameter Dispersion Corrected ' ...
                'for 3D Heterogeneity'],'color',fgcolor);
        case 'decay'
            ylabel(ax,'Decay Constant','color',fgcolor);
            title(ax,'Decay Constant Dispersion','color',fgcolor);
        case 'cdecay'
            ylabel(ax,'Decay Constant','color',fgcolor);
            title(ax,['Decay Constant Dispersion without ' ...
                     'Geometrical Spreading'],'color',fgcolor);
    end
end

end
