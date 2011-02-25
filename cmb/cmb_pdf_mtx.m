function [varargout]=cmb_pdf_mtx(pf,field,p,x)
%CMB_PDF_MTX    Plot probability space of slowness or decay constant
%
%    Usage:    cmb_pdf_mtx(pf,field)
%              cmb_pdf_mtx(pf,field,ti,xi)
%              mtx=cmb_pdf_mtx(...)
%              [ti,xi,mtx]=cmb_pdf_mtx(...)
%
%    Description:
%     CMB_PDF_MTX(PF,FIELD) plots a PDF (probability density function) of
%     the measurements indicated by FIELD stored in the profile struct PF.
%     PF is expected to be in the format created by SLOWDECAYPROFILES or
%     SLOWDECAYPAIRS.  FIELD is one of the following:
%      'slow', 'cslow', 'decay', 'cdecay'
%     where the 'c' in front means that the measurement has corrections
%     applied for 3D heterogeneous structure of the Earth.  See
%     SLOWDECAYPROFILES & SLOWDECAYPAIRS for more details.  The PDF is
%     automatically scaled using the input data ranges (some padding is
%     added to the y-range to account for the errors in the data).  The PDF
%     is 100x100 pixels.
%
%     CMB_PDF_MTX(PF,FIELD,TI,XI) allows specifying the period and
%     measurement points at which PDF estimation is performed.  Values are
%     estimated assuming the error is gaussian (the given error is 1
%     standard deviation and the measurement itself is the mean).  XI must
%     be regularly spaced (this is so pixel width can be established for
%     probability determination).  Both TI & XI must be vectors.
%
%     PDF=CMB_PDF_MTX(...) outputs the PDF matrix.  PDF is NFxNX in size
%     where NF is the number of frequency points and NX is the number of
%     points in the measurement value space.  Use the next form if you did
%     not supply FI & XI.
%
%     [TI,XI,PDF]=CMB_PDF_MTX(...) also returns the locations of the PDF in
%     period (TI) and measurement value (XI).  See the Examples section
%     below for how to plot this information
%
%    Notes:
%     - Can be quite inaccurate when errors are subpixel.  
%
%    Examples:
%     % How to plot the optional outputs:
%     [p,x,pdf]=cmb_pdf_mtx(pf,'cslow');
%     fh=figure('color','k'); ax=axes('parent',fh);
%     imagesc(p,x,pdf','parent',ax);
%     colormap(ax,fire);
%     set(ax,'ydir','normal','xcolor','w','ycolor','w');
%     xlabel(ax,'Period (sec)');
%     ylabel(ax,'Ray Parameter (s/^o)')
%     title(ax,'Ray Parameter Dispersion Corrected for 3D Heterogeneity');
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 17, 2011 at 13:35 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% check input
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

% get values & errors
x0=[pf.(field)].';
e0=[pf.([field 'err'])].';

% get min/max
pmin=min(p0(:)); pmax=max(p0(:));
xmin=min(x0-median(e0)); xmax=max(x0+median(e0));
xpad=(xmax-xmin)/4;

% default period/value ranges
% 100 pixels each way, scaled to data
if(nargin<3 || isempty(p)); p=linspace(pmin,pmax,100); end
if(nargin<4 || isempty(x)); x=linspace(xmin-xpad,xmax+xpad,100); end

% special positive scalar integer input
if(isscalar(p) && p==fix(p) && p>0); p=linspace(pmin,pmax,p); end
if(isscalar(x) && x==fix(x) && x>0); x=linspace(xmin-xpad,xmax+xpad,x); end

% check period/value
% - value needs to be regularly spaced
spacing=abs(unique(single(diff(x))));
if(~isreal(p) || ~isvector(p) || any(p<=0))
    error('seizmo:cmb_pdf_mtx:badInput',...
        'TI must be a vector of positive real-valued periods in seconds!');
elseif(~isreal(x) || ~isvector(x) || ~isscalar(spacing))
    error('seizmo:cmb_pdf_mtx:badInput',...
        'XI must be a vector of evenly spaced reals!');
end

% preallocate output
mtx=zeros(numel(p),numel(x));

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
            +2*gaussian(x,x0(idx(j)),e0(idx(j))) ...
            +gaussian(x-spacing/2,x0(idx(j)),e0(idx(j))) ...
            +gaussian(x+spacing/2,x0(idx(j)),e0(idx(j)));
    end
    
    % finish average & integration
    mtx(i,:)=mtx(i,:)/nidx/4*spacing*100;
end

% output depends on number out
if(nargout>1)
    varargout={p x mtx};
elseif(nargout)
    varargout={mtx};
else
    % plot output
    fh=figure('color','k'); ax=axes('parent',fh);
    imagesc(p,x,mtx','parent',ax);
    colormap(ax,fire);
    set(ax,'ydir','normal','color','k','xcolor','w','ycolor','w');
    set(ax,'xminortick','on','yminortick','on');
    set(ax,'xgrid','on','ygrid','on','linewidth',1);
    xlabel(ax,'Period (s)');
    cb=colorbar('peer',ax);
    ylabel(cb,'% Probability');
    set(cb,'linewidth',1);
    switch field
        case 'slow'
            ylabel(ax,'Ray Parameter (s/^o)');
            title(ax,'Ray Parameter Dispersion','color','w');
        case 'cslow'
            ylabel(ax,'Ray Parameter (s/^o)');
            title(ax,['Ray Parameter Dispersion Corrected ' ...
                'for 3D Heterogeniety'],'color','w');
        case 'decay'
            ylabel(ax,'Decay Constant');
            title(ax,'Decay Constant Dispersion','color','w');
        case 'cdecay'
            ylabel(ax,'Decay Constant');
            title(ax,['Decay Constant Dispersion Corrected ' ...
                     'for Geometrical Spreading'],'color','w');
    end
end

end
