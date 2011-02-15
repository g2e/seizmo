function [fh]=plot_yrmo_volumes(band,cmp,dblim,zerodb)
%PLOT_YRMO_VOLUMES    Makes 4x3 grid of year-month fk spectra plots
%
%    Usage:    fh=plot_yrmo_volumes(band,cmp,dblim,zerodb)
%
%    Description:
%     FH=PLOT_YRMO_VOLUMES(BAND,CMP,DBLIM,ZERODB) reads fk volumes in
%     the current directory created by MAKE_YRMO_Z_VOLUMES or
%     MAKE_YRMO_HORZ_VOLUMES and creates a 4x3 set of plots in a figure for
%     each year available (one plot per month) and for each frequency band
%     in BAND.  Axes of missing months are not shown.  So if you created fk
%     volumes spanning 3 different years then this will create 3 figures
%     per frequency band.  BAND should be a Nx2 array of [LOW HIGH] in Hz.
%     CMP should be one of 'Z', 'R', or 'T'.  See PLOTFKMAP for details
%     about DBLIM & ZERODB.  There are defaults for each of inputs: BAND
%     is a series of frequency bands between .01 & .2 Hz, CMP is 'Z', and
%     ZERODB is 'median' and DBLIM is [0 6].
%
%    Notes:
%     - Figures are saved as:
%        fkyrmo_CMP_YEAR_BANDLOs-BANDHIs_ZERODB_DBLIM1db-DBLIM2db_orig.fig
%
%    Examples:
%     % Create fk volumes, create plots, and prepare for printing:
%     make_yrmo_z_volumes(stack_dir);
%     fh=plot_yrmo_volumes;
%     adjust_yrmo_plots(fh,'My Array Name');
%
%    See also: ADJUST_YRMO_PLOTS, MAKE_YRMO_HORZ_VOLUMES, 
%              MAKE_YRMO_Z_VOLUMES

%     Version History:
%        Feb. 14, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2011 at 16:05 GMT

% todo:

% check nargin
error(nargchk(0,4,nargin));

% defaults
if(nargin<1 || isempty(band))
    band=1./[100 50; 50 40; 40 30; 30 25; 25 20;
        20 15; 15 10; 10 7.5; 7.5 5];
end
if(nargin<2 || isempty(cmp)); cmp='z'; end
if(nargin<3 || isempty(dblim)); dblim=[0 6]; end
if(nargin<4 || isempty(zerodb)); zerodb='median'; end

% check inputs
if(~isreal(band) || size(band,2)~=2)
    error('seizmo:plot_yrmo_volumes:badInput',...
        'BAND must be a Nx2 real-valued array!');
end
if(~isstring(cmp) || ~ismember(lower(cmp),{'z' 'r' 't'}))
    error('seizmo:plot_yrmo_volumes:badInput',...
        'CMP must be one of ''Z'', ''R'', or ''T''!');
end

% find all volumes matching name scheme
% fkvol.(zrt).year.mo.mat
% ==> remove fkvol.(zrt).year.day.mat
files=xdir(['fkvol.' lower(cmp) '.*.*.mat']);
names={files.name}';
names(cellfun('size',names,2)~=19)=[];
mo0=char(names);
yr0=str2double(cellstr(mo0(:,9:12)));
mo0=str2double(cellstr(mo0(:,14:15)));
uyr=unique(yr0); uyrstr=cellstr(num2str(uyr));
nyr=numel(uyr);

% initialize plots and axes
nb=size(band,1);
fh=nan(nb,nyr); ax=cell(nb,nyr);
for i=1:nb
    for j=1:nyr
        userdata.year=uyr(j);
        userdata.band=band(i,:);
        userdata.cmp=cmp;
        userdata.dblim=dblim;
        userdata.zerodb=zerodb;
        fh(i,j)=figure('color','w','tag','fkyrmo','name',...
            ['YEAR ' uyrstr{j} '  ' num2str(1./band(i,2)) 's-' ...
            num2str(1./band(i,1)) 's'],'userdata',userdata);
        ax{i,j}=makesubplots(4,3,[],'align',...
            'parent',fh(i,j),'visible','off','tag','fkmap');
        drawnow;
    end
end

% month plot labels
month={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' ...
    'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

% read in yrmo volumes and plot
for i=1:numel(names)
    % read in volume
    vol=load(names{i});
    iy=yr0(i)==uyr;
    im=mo0(i);
    
    % loop over each band
    for j=1:nb
        % get fkmap
        map=fkvol2map(vol,band(j,:));
        
        % plot it
        set(ax{j,iy}(im),'visible','on');
        plotfkmap(map,dblim,zerodb,'k','w',ax{j,iy}(im));
        
        % fix title
        title(ax{j,iy}(im),[month{im} ' ' uyrstr{iy}]);
    end
    
    drawnow;
end

% save figures
for i=1:nb
    for j=1:nyr
        saveas(fh(i,j),['fkyrmo_' lower(cmp) '_' uyrstr{j} '_' ...
            num2str(1./band(i,2)) 's-' num2str(1./band(i,1)) 's_' ...
            zerodb '_' num2str(dblim(1)) 'db-' num2str(dblim(2)) ...
            'db_orig.fig'],'fig');
    end
end

end
