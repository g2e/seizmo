function [grp,ax]=popcut(data,grp)
%POPCUT    Cut low POP values interactively
%
%    Usage:    [grp,ax]=popcut(data,grp)
%
%    Description:
%     [GRP,AX]=POPCUT(DATA,GRP) creates an interactive plot of cluster
%     population in the struct GRP where the user is expected to
%     indicate the cutoff for removing low population groups.  The default
%     cutoff is set by GRP.popcut and is altered when the user left clicks
%     the plot.  Below the interactive population plot are subplots showing
%     the waveforms of the clusters with an acceptable population.  To
%     complete the cut and exit, the user must middle-click the plot.  The
%     logical indices of the groups above the cluster population cutoff
%     are returned in GRP.good.  AX is an array of the axes handles in the
%     interactive plot.  The first handle is the population plot and the
%     remaining handles are for the individual cluster waveform plots.
%
%    Notes:
%
%    Examples:
%     % Cluster some data interactively then look at the populations:
%     grp=usercluster(data,cg);
%     grp=popcut(data,grp);
%
%     % Change the initial population cutoff to 10 members:
%     grp.popcut=10;
%     grp=popcut(data,grp);
%
%    See also: ARRCUT, AMPCUT, ERRCUT, SNRCUT, PLOTPOP, PLOTCLUSTERS,
%              USERCLUSTER, SELECTCLUSTERS

%     Version History:
%        Sep. 18, 2010 - initial version
%        Sep. 21, 2010 - altered inputs/outputs
%        Jan.  6, 2011 - proper ginput handling, use key2zoompan
%        Mar. 31, 2011 - more columns than rows is more appealing
%        Apr.  3, 2012 - use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 23:55 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% number of records
nrecs=numel(data);

% check cluster struct
if(~isstruct(grp) ...
        || any(~ismember({'T' 'color' 'popcut'},fieldnames(grp))) ...
        || numel(grp.T)~=nrecs ...
        || ~isequal([max(grp.T) 3],size(grp.color)))
    error('seizmo:plotpop:badInput',...
        'GRP must be a struct with fields T & color!');
end

% cutoff
cutoff=grp.popcut;

% check cutoff
if(~isreal(cutoff) || ~isscalar(cutoff) || cutoff<0)
    error('seizmo:popcut:badInput',...
        'CUTOFF must be a positive value!');
end

% get populations
pop=histc(grp.T,1:max(grp.T));

% groups above cutoff
good=pop>=cutoff;
fgood=find(good);
ngrp=sum(good);

% initialize figure
fh=figure('color','k');

% only draw waveform subplots if some meet cutoff
if(ngrp)
    % figure out subplot indexing
    nrows=round(sqrt(ngrp));
    ncols=ceil(ngrp/nrows);
    anrows=ceil(nrows*3/2);
    sidx=ncols*(anrows-nrows)+1;
    ax(2:1+ngrp)=makesubplots(anrows,ncols,sidx:sidx+ngrp-1,...
        'parent',fh,'color','k');
    
    % make waveform plots
    ax(2:1+ngrp)=plotclusters(data,grp,...
        'clusters',fgood,'ax',ax(2:1+ngrp));
end

% make the top plot second
% - we do this b/c creating a single subplot with makesubplots
%   somehow turns on the 'align' flag ... still figuring this one out
ax(1)=makesubplots(3,1,1,'parent',fh,'color','k');

% plot pop
plotpop(grp,ax(1));

% add in cutoff bar
x=get(ax(1),'xlim');
hold(ax(1),'on');
hcut=plot(x,[cutoff cutoff],'r--','linewidth',2);
hold(ax(1),'off');

% adjust tick labels
set(ax(1),'xtick',fgood,'xticklabel',fgood);

% title
title(ax(1),['Left Click = Change Cutoff, ' ...
    'Middle Click = Implement Cut, ' ...
    'Cutoff = ' num2str(cutoff)],...
    'fontsize',10,'fontweight','bold','color','w');

% let user adjust the limits
unhappy=true;
while(unhappy)
    % get population cutoff
    axis(ax(1));
    try
        [x,y,button]=ginput(1);
    catch
        % plot closed - break out
        break;
    end
    
    % act based on button
    switch button
        case 1
            % skip if not correct axis
            if(ax(1)~=gca); continue; end
            
            % get cutoff
            cutoff=abs(y);
            
            % adjust limits
            set(hcut,'ydata',[cutoff cutoff]);
            
            % reset title
            set(get(ax(1),'title'),'string',...
                ['Left Click = Change Cutoff, ' ...
                'Middle Click = Implement Cut, ' ...
                'Cutoff = ' num2str(cutoff)]);
            
            % new groups above cutoff
            good=pop>=cutoff;
            ngrp=sum(good);
            fgood=find(good);
            
            % reset xticks
            set(ax(1),'xtick',fgood,'xticklabel',fgood);
            
            % delete cluster plots
            delete(ax(2:end));
            ax(2:end)=[];
            
            % skip if no groups
            if(~ngrp); continue; end
            
            % figure out subplot indexing
            nrows=round(sqrt(ngrp));
            ncols=ceil(ngrp/nrows);
            anrows=ceil(nrows*3/2);
            sidx=ncols*(anrows-nrows)+1;
            ax(2:1+ngrp)=makesubplots(anrows,ncols,sidx:sidx+ngrp-1,...
                'parent',fh,'color','k');
            
            % make waveform plots
            ax(2:1+ngrp)=plotclusters(data,grp,...
                'clusters',fgood,'ax',ax(2:1+ngrp));
        case 2
            % skip if not correct axis
            if(ax(1)~=gca); continue; end
            
            % break out
            unhappy=false;
        otherwise
            key2zoompan(button);
    end
end

% add to struct
grp.popcut=cutoff;
grp.good=good;

end
