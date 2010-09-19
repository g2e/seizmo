function [good,varargout]=popcut(data,grp,cutoff)
%POPCUT    Cut low POP values interactively
%
%    Usage:    good=popcut(data,grp)
%              good=popcut(data,grp,cutoff)
%              [good,cutoff]=popcut(...)
%              [good,cutoff,ax]=popcut(...)
%
%    Description:
%     GOOD=POPCUT(DATA,GRP) creates an interactive plot of cluster
%     population in the struct GRP where the user is expected to
%     indicate the cutoff for removing low population groups.  The default
%     cutoff is a population of 3 and is altered when the user left clicks
%     the plot.  Below the interactive population plot are subplots showing
%     the waveforms of the clusters with an acceptable population.  To
%     complete the cut and exit, the user must middle-click the plot.  The
%     indices of the good groups are returned in GOOD.
%
%     GOOD=POPCUT(DATA,GRP,CUTOFF) adjusts the default population cutoff to
%     CUTOFF.  The default POP cutoff is 3.
%
%     [GOOD,CUTOFF]=POPCUT(...) returns the final cutoff used for removing
%     the outliers.
%
%     [GOOD,CUTOFF,AX]=POPCUT(...) returns the axes of the plots.  The
%     first handle is the population plot and the remaining handles are for
%     the individual cluster waveform plots.
%
%    Notes:
%
%    Examples:
%     % Cluster some data interactively then look at the populations:
%     grp=usercluster(data,cg);
%     popcut(data,grp);
%
%    See also: ARRCUT, AMPCUT, ERRCUT, SNRCUT, PLOTPOP, PLOTCLUSTERS,
%              USERCLUSTER

%     Version History:
%        Sep. 18, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 18, 2010 at 20:00 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% defaults
if(nargin<3 || isempty(cutoff)); cutoff=3; end

% check data structure
versioninfo(data,'dep');

% number of records
nrecs=numel(data);

% check cluster struct
if(~isstruct(grp) || any(~ismember({'T' 'color'},fieldnames(grp))) ...
        || numel(grp.T)~=nrecs || ~isequal([max(grp.T) 3],size(grp.color)))
    error('seizmo:plotpop:badInput',...
        'GRP must be a struct with fields T & color!');
end

% check cutoff
if(~isreal(cutoff) || ~isscalar(cutoff) || cutoff<0)
    error('seizmo:popcut:badInput',...
        'CUTOFF must be a positive value!');
end

% get populations
pop=histc(grp.T,1:max(grp.T));

% groups above cutoff
good=find(pop>=cutoff);
ngrp=numel(good);

% initialize figure
fh=figure('color','k');

% only draw waveform subplots if some meet cutoff
if(ngrp)
    % figure out subplot indexing
    ncols=fix(sqrt(ngrp));
    nrows=ceil(ngrp/ncols);
    anrows=ceil(nrows*3/2);
    sidx=ncols*(anrows-nrows)+1;
    ax(2:1+ngrp)=makesubplots(anrows,ncols,sidx:sidx+ngrp-1,...
        'parent',fh,'color','k');
    
    % make waveform plots
    ax(2:1+ngrp)=plotclusters(data,grp,'clusters',good,'ax',ax(2:1+ngrp));
end

% make the top plot second
% - we do this b/c creating a single subplot with makesubplots
%   somehow turns on the 'align' flag ... still figuring this one out
ax(1)=makesubplots(3,1,1,'parent',fh,'color','k');

% plot pop
plotpop(grp,ax(1))

% add in cutoff bar
x=get(ax(1),'xlim');
hold(ax(1),'on');
hcut=plot(x,[cutoff cutoff],'r--','linewidth',2);
hold(ax(1),'off');

% adjust tick labels
set(ax(1),'xtick',good,'xticklabel',good);

% title
title(ax(1),['Left Click = Change Cutoff, ' ...
    'Middle Click = Implement Cut, ' ...
    'Cutoff = ' num2str(cutoff)],...
    'fontsize',10,'fontweight','bold','color','w');

% let user adjust the limits
unhappy=true;
while(unhappy)
    axis(ax(1));
    [x,y,button]=ginput(1);
    switch button
        case 1
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
            good=find(pop>=cutoff);
            ngrp=numel(good);
            
            % reset xticks
            set(ax(1),'xtick',good,'xticklabel',good);
            
            % delete cluster plots
            delete(ax(2:end));
            ax(2:end)=[];
            
            % skip if no groups
            if(~ngrp); continue; end
            
            % figure out subplot indexing
            ncols=fix(sqrt(ngrp));
            nrows=ceil(ngrp/ncols);
            anrows=ceil(nrows*3/2);
            sidx=ncols*(anrows-nrows)+1;
            ax(2:1+ngrp)=makesubplots(anrows,ncols,sidx:sidx+ngrp-1,...
                'parent',fh,'color','k');
            
            % make waveform plots
            ax(2:1+ngrp)=plotclusters(data,grp,...
                'clusters',good,'ax',ax(2:1+ngrp));
        case 2
            good=find(pop>cutoff);
            unhappy=false;
    end
end

% output if desired
if(nargout>1); varargout={cutoff ax}; end

end
