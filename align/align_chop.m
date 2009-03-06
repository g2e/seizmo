function [data,CONF,good]=align_chop(data,CONF)
%CHOP    Window data

% check input
error(nargchk(2,2,nargin))

% 
good=1;

% window around signal
if(CONF.WINDOW)
    disp('WINDOWING DATA')
    
    % will the user refine the window
    if(CONF.USERWIN)
        % default window (aligns and focuses on arrivaltimes) 
        data=cutim(data,CONF.TTFIELD,CONF.DEFWIN(1),CONF.TTFIELD,CONF.DEFWIN(2),CONF.FILL,CONF.FILLER);
        
        % remove mean
        data=rmean(data);
        
        % user window
        [data]=userwindow(data);
    else
        % auto-window (for the brave)
        data=cutim(data,CONF.TTFIELD,CONF.SIGWIN(1),CONF.TTFIELD,CONF.SIGWIN(2),CONF.FILL,CONF.FILLER);
    end
end
    % check window lengths are equal
    if (~isempty(find(gh(data,'npts')~=gh(data(1),'npts'),1)))
        error('npts not equal for all data')
    end
    CONF.NPTS=gh(data(1),'npts');
    
    % write window info
    if(exist([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window'],'file'))
        delete([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window']);
    end
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window'],...
    ['interactive = ' num2str(CONF.USERWIN)],'');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window'],...
    ['delta = ' num2str(CONF.DELTA)],'-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window'],...
    ['npts = ' num2str(CONF.NPTS)],'-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window'],...
    ['length = ' num2str(CONF.NOWWIN(2)-CONF.NOWWIN(1))],...
    '-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window'],...
    ['tapertype = ' num2str(CONF.TAPERTYPE)],'-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window'],...
    ['taperhw = ' num2str(CONF.TAPERHW)],'-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window'],...
    ['taperopt = ' num2str(CONF.TAPEROPT)],'-append','delimiter','');
pad=ones(length(data),5)*32;
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.window'],...
    [char({data.name}') pad num2str(gh(data,'b')) ...
    pad num2str(gh(data,'e'))],'-append','delimiter','');

end
