function [fh]=plot_monthly_volumes(band,cmp,dblim,zerodb)
%PLOT_MONTHLY_VOLUMES    Makes 4x3 grid of monthly fk spectra plots
%
%    Usage:    fh=plot_monthly_volumes(band,cmp,dblim,zerodb)
%
%    Description:
%     FH=PLOT_MONTHLY_VOLUMES(BAND,CMP,DBLIM,ZERODB) reads fk volumes in
%     the current directory created by MAKE_MONTHLY_Z_VOLUMES or
%     MAKE_MONTHLY_HORZ_VOLUMES and creates a 4x3 set of plots (one plot
%     per month) in a figure (one figure per band in BAND).  BAND should be
%     a Nx2 array of [LOW HIGH] in Hz.  CMP should be one of 'Z', 'R', or
%     'T'.  See PLOTFKMAP for details about DBLIM & ZERODB.  There are
%     defaults for each of inputs: BAND is a series of frequency bands
%     between .01 & .2 Hz, CMP is 'Z', and ZERODB is 'median' and DBLIM is
%     [0 6].
%
%    Notes:
%     - Figures are saved as:
%        fkmonthly_CMP_BANDLOs-BANDHIs_ZERODB_DBLIM1db-DBLIM2db_orig.fig
%
%    Examples:
%     % Create fk volumes, create plots, and prepare for printing:
%     make_monthly_z_volumes(stack_dir);
%     fh=plot_monthly_volumes;
%     adjust_monthly_plots(fh,'My Array Name');
%
%    See also: MAKE_MONTHLY_HORZ_VOLUMES, MAKE_MONTHLY_Z_VOLUMES

%     Version History:
%        Oct. 10, 2010 - initial version
%        Oct. 11, 2010 - several bug fixes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 11, 2010 at 16:05 GMT

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
    error('seizmo:plot_monthly_volumes:badInput',...
        'BAND must be a Nx2 real-valued array!');
end
if(~isstring(cmp) || ~ismember(lower(cmp),{'z' 'r' 't'}))
    error('seizmo:plot_monthly_volumes:badInput',...
        'CMP must be one of ''Z'', ''R'', or ''T''!');
end

% initialize plots and axes
nb=size(band,1);
fh=nan(nb,1); ax=cell(nb,1); tax=ax;
for i=1:nb
    userdata.band=band(i,:);
    userdata.cmp=cmp;
    userdata.dblim=dblim;
    userdata.zerodb=zerodb;
    fh(i)=figure('color','w','tag','fkmonthly','name',...
        [num2str(1./band(i,2)) 's-' num2str(1./band(i,1)) 's'],...
        'userdata',userdata);
    ax{i}=makesubplots(4,3,[],'align','parent',fh(i));
    tax{i}=ax{i}';
    drawnow;
end

% month plot labels
month={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' ...
    'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

% read in monthly volumes and plot
for i=1:12
    % read in volume
    vol=load(['fkvol.' lower(cmp) '.' num2str(i,'%02d') '.mat']);
    
    % loop over each band
    for j=1:nb
        % get fkmap
        map=fkvol2map(vol,band(j,:));
        
        % plot it
        plotfkmap(map,dblim,zerodb,'k','w',ax{j}(i));
        
        % fix title
        title(ax{j}(i),month{i});
    end
    
    drawnow;
end

% save figures
for i=1:nb
    saveas(fh(i),['fkmonthly_' lower(cmp) '_' num2str(1./band(i,2)) ...
        's-' num2str(1./band(i,1)) 's_' zerodb '_' num2str(dblim(1)) ...
        'db-' num2str(dblim(2)) 'db_orig.fig'],'fig');
end

end
