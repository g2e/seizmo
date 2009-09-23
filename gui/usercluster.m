function [cutoff,color,perm]=usercluster(data,Z,cutoff,varargin)
%USERCLUSTER    Interactively cluster SEIZMO records
%
%    Usage:    data=usercluster(data)
%              data=usercluster(data,Z)
%              data=usercluster(data,Z,'field',value,...)
%              [data,cutoff]=usercluster(...)
%              [data,cutoff,fh]=usercluster(...)
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also: userwindow, usertaper, plotdendro, cluster, linkage

%     Version History:
%        Sep. 22, 2009 - rewrite and added documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 22, 2009 at 07:30 GMT

% todo:

% cluster analysis loop
while (1)
    % default limit
    %cutoff=0.03;
    
    % dendrogram/waveforms
    [perm,color,fh,sfh]=plotdendro(Z,data,'treelimit',cutoff,varargin{:});
    drawnow;
    
    mymenu={'Please choose a cutoff in the dendrogram to indicate';
            'the maximum dissimilarity link for clusters.';
            '  '
            'USAGE';
            '================================';
            'Left click   -=> indicates cutoff position';
            'Middle click -=> finalizes the cutoff';
            '================================'};
    menu(mymenu,'OK');
    
    % loop until right click
    button=1;
    while (button~=2)
        % user picks
        figure(fh);
        subplot(sfh(1));
        [x,y,button]=ginput(1);
        
        % action on left click
        if (button==1)
            cutoff=x;
            
            % redraw plot (to show new grouping)
            delete(sfh);
            [perm,color,fh,sfh]=plotdendro(Z,data,'treelimit',cutoff,...
                'fighandle',fh,varargin{:});
        end
    end
    
    % check satisfaction
    done=menu('Keep groups?','YES','NO - TRY AGAIN','NO - DIE!');
    if(done==1)
        return;
    elseif(done==2)
        close(fh);
    elseif(done==3)
        % DEATH!
        error('user seppuku')
    end
end

end
