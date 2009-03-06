function [T,pop,big,small]=groupwaves(corr_grid,data)
%GROUPWAVES    Perform cluster analysis on correlated seismograms

% globals
global CONF GCOLOR

% check nargin
error(nargchk(2,2,nargin))

% check data structure
if(~isstruct(data))
    error('sac data is not a structure')
elseif(~isfield(data,'name') || ~isfield(data,'head') || ~isfield(data,'amp'))
    error('sac data structure does not have proper fields')
elseif(isfield(data,'time'))
    error('function does not work with unevenly sampled data')
elseif(~isvector(data))
    error('sac data structure not a vector')
end

% number of seismograms
nseis=length(data);

% check if a fraction
if(CONF.MINSIG<1)
    CONF.MINSIG=CONF.MINSIG*nseis;
end

% need 3 for statistical measures
if(CONF.MINSIG<3)
    CONF.MINSIG=3;
end

% check if nseis is low
if(nseis<CONF.MINSIG)
    disp('WARNING: number of signals is less than min group size')
end

% make vector from grid
lt=tril(corr_grid,-1);  % only need lower triangle
clear corr_grid
corr_vec=lt(lt~=0).';   % if a correlation = 0 things will crash
clear lt

% calculate heirarchial linkage
Z=linkage(1-corr_vec,CONF.CMETHOD);
clear corr_vec

% cluster analysis loop
done=2;
while (done==2)
    % user-defined dissimilarity cutoff
    if (CONF.USERCLUSTER)
        % dendrogram/waveforms
        [perm,h,sh]=plotdendrogram(Z,[data.amp]);
        drawnow;
        
        mymenu={'Please choose a cutoff in the dendrogram to indicate';
                'the maximum dissimilarity link for clusters.';
                '  '
                'USAGE';
                '================================';
                'Left click  -=> indicates cutoff position';
                'Right click -=> finalizes the cutoff';
                '================================'};
        menu(mymenu,'OK');
        
        % loop until right click
        button=1;
        while (button~=3)
            
            % user picks
            figure(h);
            subplot(sh(1));
            [x,y,button]=ginput(1);
            
            % action on left click
            if (button==1)
                CONF.MAXDIST=x;
                
                % redraw plot (to show new grouping)
                delete(sh);
                [perm,h,sh]=plotdendrogram(Z,[data.amp],h);
            end
        end
    end
    
    % ensure distance measure range: 0-1
    if(CONF.MAXDIST<0)
        CONF.MAXDIST=0;
    elseif(CONF.MAXDIST>1)
        CONF.MAXDIST=1;
    end

    % cluster based on user input
    T=cluster(Z,'cutoff',CONF.MAXDIST,'criterion','distance');
    
    % sort out groups
    ng=max(T);
    [pop,big,small]=grouper(T);
    
    % user-defined min pop cut
    if (CONF.USERCLUSTER)
        % individual coloring to group coloring
        iperm(perm)=1:nseis;
        color=GCOLOR(iperm,:);
        GCOLOR=zeros(ng,3);
        for i=1:ng
            GCOLOR(i,:)=color(find(i==T,1),:);
        end
        clear color perm iperm
        
        % pop/group plot
        [h2,sh2]=plotclusters(T,pop,big,small,[data.amp]);
        drawnow;
        
        mymenu={'Please choose a cutoff for the minimum population';
                'size for clusters to be considered.';
                '  '
                'USAGE';
                '================================';
                'Left click  -=> indicates cutoff position';
                'Right click -=> finalizes the cutoff';
                '================================'};
        menu(mymenu,'ok');
        
        % loop until right click
        button=1;
        while (button~=3)
            
            % user picks
            figure(h2(1))
            subplot(sh2)
            [x,y,button]=ginput(1);
            
            % action on left click
            if (button==1)
                CONF.MINSIG=y;
                
                % replot groups
                close(h2);
                [pop,big,small]=grouper(T);
                [h2,sh2]=plotclusters(T,pop,big,small,[data.amp]);
            end
        end
    else
        % make some group colors (for autocluster)
        GCOLOR=hsv(ng);
    end
    
    % asking if content with results - if so, exit = if not, redo cluster analysis
    if (CONF.USERCLUSTER)
        done=menu('Keep groups?','YES','NO - REDO');
        if (done==2)
            close(h);
            close(h2);
        end
    else
        done=1;
    end
end

% write group info
if(exist([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group'],'file'))
    delete([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group']);
end
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group'],...
    ['interactive = ' num2str(CONF.USERCLUSTER)],'');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group'],...
    ['method = ' CONF.CMETHOD],'-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group'],...
    ['max dist = ' num2str(CONF.MAXDIST)],'-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group'],...
    ['min grp size = ' num2str(CONF.MINSIG)],'-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group'],...
    ['num grps = ' num2str(ng)],'-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group'],...
    ['big grps = ' num2str(length(big))],'-append','delimiter','');
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group'],...
    ['small grps = ' num2str(length(small))],'-append','delimiter','');
pad=ones(length(data),5)*32;
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.group'],...
    [char({data.name}') pad num2str(T)],...
    '-append','delimiter','');
if(exist([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.linkage'],'file'))
    delete([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.linkage']);
end
dlmwrite([strtrim(CONF.PHASE) '.' num2str(CONF.ID) '.linkage'],Z,' ');


end
