function [varargout]=cmb_pdf_mtx(pf,field,f,x)
%CMB_PDF_MTX    Plot probability space of slowness or decay constant
%
%    Usage:    cmb_pdf_mtx(pf,field)
%              cmb_pdf_mtx(pf,field,fi,xi)
%              mtx=cmb_pdf_mtx(...)
%              [fi,xi,mtx]=cmb_pdf_mtx(...)
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
%     CMB_PDF_MTX(PF,FIELD,FI,XI) allows specifying the frequency and
%     measurement points at which PDF estimation is performed.  Values are
%     estimated assuming the error is gaussian (the given error is 1
%     standard deviation and the measurement itself is the mean).  XI must
%     be regularly spaced (this is so pixel width can be established for
%     probability determination).  Both FI & XI must be vectors.
%
%     PDF=CMB_PDF_MTX(...) outputs the PDF matrix.  PDF is NFxNX in size
%     where NF is the number of frequency points and NX is the number of
%     points in the measurement value space.  Use the next form if you did
%     not supply FI & XI.
%
%      See the Examples
%     section below for how to plot this information
%
%     [FI,XI,PDF]=CMB_PDF_MTX(...) also returns the locations of the PDF in
%     frequency (FI) and measurement value (XI).  See the Examples section
%     below for how to plot this information
%
%    Notes:
%
%    Examples:
%     % How to plot the output:
%     [f,x,pdf]=cmb_pdf_mtx(pf,'cslow');
%     fh=figure('color','k'); ax=axes('parent',fh);
%     imagesc(f,x,pdf','parent',ax);
%     colormap(ax,fire);
%     set(ax,'ydir','normal','xcolor','w','ycolor','w');
%     xlabel(ax,'Freq (Hz)');
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  1, 2011 at 13:35 GMT

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

% get frequencies
f0=cell2mat({pf.freq}');

% get values & errors
x0=[pf.(field)].';
e0=[pf.([field 'err'])].';

% get min/max
fmin=min(f0(:)); fmax=max(f0(:));
xmin=min(x0(:)); xmax=max(x0(:));
xpad=(xmax-xmin)/4;

% default freq/value ranges
% 100 pixels each way, scaled to data
if(nargin<3 || isempty(f)); f=linspace(fmin,fmax,100); end
if(nargin<4 || isempty(x)); x=linspace(xmin-xpad,xmax+xpad,100); end

% check freq/value
% - value needs to be regularly spaced
spacing=abs(unique(single(diff(x))));
if(~isreal(f) || ~isvector(f) || any(f<=0))
    error('seizmo:cmb_pdf_mtx:badInput',...
        'FI must be a vector of positive real-valued frequencies in Hz!');
elseif(~isreal(x) || ~isvector(x) || ~isscalar(spacing))
    error('seizmo:cmb_pdf_mtx:badInput',...
        'XI must be a vector of evenly spaced reals!');
end

% preallocate output
mtx=zeros(numel(f),numel(x));

% loop over freq
for i=1:numel(f)
    % find freq ranges encompassing current freq
    idx=find(f(i)>=f0(:,1) & f(i)<=f0(:,2));
    nidx=numel(idx);
    
    % skip if none
    if(~nidx); continue; end
    
    % loop over values
    for j=1:nidx
        mtx(i,:)=mtx(i,:)+100*gaussian(x,x0(idx(j)),e0(idx(j)));
    end
    
    % finish average
    mtx(i,:)=mtx(i,:)*spacing/nidx;
end

% output depends on number out
if(nargout>1)
    varargout={f x mtx};
elseif(nargout)
    varargout={mtx};
else
    % plot output
    switch field
        case 'slow'
            fh=figure('color','k'); ax=axes('parent',fh);
            imagesc(f,x,mtx','parent',ax);
            colormap(ax,fire);
            set(ax,'ydir','normal','xcolor','w','ycolor','w');
            xlabel(ax,'Freq (Hz)');
            ylabel(ax,'Ray Parameter (s/^o)');
            title(ax,'Ray Parameter Dispersion','color','w');
            cb=colorbar('peer',ax);
            ylabel(cb,'% Probability');
            grid(ax,'on');
        case 'cslow'
            fh=figure('color','k'); ax=axes('parent',fh);
            imagesc(f,x,mtx','parent',ax);
            colormap(ax,fire);
            set(ax,'ydir','normal','xcolor','w','ycolor','w');
            xlabel(ax,'Freq (Hz)');
            ylabel(ax,'Ray Parameter (s/^o)');
            title(ax,['Ray Parameter Dispersion Corrected ' ...
                'for 3D Heterogeniety'],'color','w');
            cb=colorbar('peer',ax);
            ylabel(cb,'% Probability');
            grid(ax,'on');
        case 'decay'
            fh=figure('color','k'); ax=axes('parent',fh);
            imagesc(f,x,mtx','parent',ax);
            colormap(ax,fire);
            set(ax,'ydir','normal','xcolor','w','ycolor','w');
            xlabel(ax,'Freq (Hz)');
            ylabel(ax,'Decay Constant');
            title(ax,'Decay Constant Dispersion','color','w');
            cb=colorbar('peer',ax);
            ylabel(cb,'% Probability');
            grid(ax,'on');
        case 'cdecay'
            fh=figure('color','k'); ax=axes('parent',fh);
            imagesc(f,x,mtx','parent',ax);
            colormap(ax,fire);
            set(ax,'ydir','normal','xcolor','w','ycolor','w');
            xlabel(ax,'Freq (Hz)');
            ylabel(ax,'Decay Constant');
            title(ax,['Decay Constant Dispersion Corrected ' ...
                     'for Geometrical Spreading'],'color','w');
            cb=colorbar('peer',ax);
            ylabel(cb,'% Probability');
            grid(ax,'on');
    end
end

end
