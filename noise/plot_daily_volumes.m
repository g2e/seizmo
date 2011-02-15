function [fh]=plot_daily_volumes(band,cmp,dblim,zerodb)
%PLOT_DAILY_VOLUMES    Makes weekly grids of daily fk spectra plots
%
%    Usage:    fh=plot_daily_volumes(band,cmp,dblim,zerodb)
%
%    Description:
%     FH=PLOT_DAILY_VOLUMES(BAND,CMP,DBLIM,ZERODB) reads fk volumes in
%     the current directory created by MAKE_DAILY_Z_VOLUMES or
%     MAKE_DAILY_HORZ_VOLUMES and creates a set of plots in a figure for
%     each week available and for each frequency band in BAND.  Axes of
%     missing days in a week are not shown.  So if you created fk volumes
%     spanning 3 weeks then this will create 3 figures per frequency band.
%     BAND should be a Nx2 array of [LOW HIGH] in Hz.  CMP should be one of
%     'Z', 'R', or 'T'.  See PLOTFKMAP for details about DBLIM & ZERODB.
%     There are defaults for each of inputs: BAND is a series of frequency
%     bands between .01 & .2 Hz, CMP is 'Z', and ZERODB is 'median' and
%     DBLIM is [0 6].
%
%    Notes:
%     - Figures are saved as:
%        fkdaily_CMP_YEAR.DAY-YEAR.DAY_BANDLOs-BANDHIs_ (rest on next line)
%         ZERODB_DBLIM1db-DBLIM2db_orig.fig
%
%    Examples:
%     % Create fk volumes, create plots, and prepare for printing:
%     make_daily_z_volumes(stack_dir);
%     fh=plot_daily_volumes;
%     adjust_daily_plots(fh,'My Array Name');
%
%    See also: ADJUST_DAILY_PLOTS, MAKE_DAILY_HORZ_VOLUMES, 
%              MAKE_DAILY_Z_VOLUMES

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
    error('seizmo:plot_daily_volumes:badInput',...
        'BAND must be a Nx2 real-valued array!');
end
if(~isstring(cmp) || ~ismember(lower(cmp),{'z' 'r' 't'}))
    error('seizmo:plot_daily_volumes:badInput',...
        'CMP must be one of ''Z'', ''R'', or ''T''!');
end

% month plot labels
month=['Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; ...
    'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'];

% find all volumes matching name scheme
% fkvol.(zrt).year.day.mat
% ==> remove fkvol.(zrt).year.mo.mat
files=xdir(['fkvol.' lower(cmp) '.*.*.mat']);
names={files.name}';
names(cellfun('size',names,2)~=20)=[];
dy0=char(names);
yr0=str2double(cellstr(dy0(:,9:12)));
dy0=str2double(cellstr(dy0(:,14:16)));
cal0=doy2cal([yr0 dy0]);
cal0str=strcat(num2str(cal0(:,1)),'-',month(cal0(:,2),:),...
    '-',num2str(cal0(:,3),'%02d'));
sr0=datenum(cal0);
wk0=floor((sr0-2)/7);
[dow0,dow0str]=weekday(sr0);
uwk=unique(wk0);
nwk=numel(uwk);

% get start/end dates for each week
stdate=serial2gregorian(uwk*7+2,'doydate');
endate=serial2gregorian(uwk*7+8,'doydate');
stdatestr=strcat(num2str(stdate(:,1)),'.',num2str(stdate(:,2),'%03d'));
endatestr=strcat(num2str(endate(:,1)),'.',num2str(endate(:,2),'%03d'));

% initialize plots and axes
nb=size(band,1);
fh=nan(nb,nwk); ax=cell(nb,nwk);
for i=1:nb
    for j=1:nwk
        userdata.week=[stdatestr(j,:) '-' endatestr(j,:)];
        userdata.band=band(i,:);
        userdata.cmp=cmp;
        userdata.dblim=dblim;
        userdata.zerodb=zerodb;
        fh(i,j)=figure('color','w','tag','fkdaily','name',...
            [userdata.week '  ' num2str(1./band(i,2)) ...
            's-' num2str(1./band(i,1)) 's'],'userdata',userdata);
        ax{i,j}=makesubplots(3,3,[],'align',...
            'parent',fh(i,j),'visible','off','tag','fkmap');
        drawnow;
    end
end

% read in daily volumes and plot
for i=1:numel(names)
    % read in volume
    vol=load(names{i});
    iw=wk0(i)==uwk;
    id=dow0(i);
    
    % make title string
    % Mon 2006-Jan-20 (020)
    titstr=[dow0str(i,:) ' ' cal0str(i,:) ' (' num2str(dy0(i),'%03d') ')'];
    
    % loop over each band
    for j=1:nb
        % get fkmap
        map=fkvol2map(vol,band(j,:));
        
        % plot it
        set(ax{j,iw}(id),'visible','on');
        plotfkmap(map,dblim,zerodb,'k','w',ax{j,iw}(id));
        
        % fix title
        title(ax{j,iw}(id),titstr);
    end
    
    drawnow;
end

% save figures
for i=1:nb
    for j=1:nwk
        saveas(fh(i,j),['fkdaily_' lower(cmp) '_' stdatestr(j,:) '-' ...
            endatestr(j,:) '_' num2str(1./band(i,2)) 's-' ...
            num2str(1./band(i,1)) 's_' zerodb '_' num2str(dblim(1)) ...
            'db-' num2str(dblim(2)) 'db_orig.fig'],'fig');
    end
end

end
