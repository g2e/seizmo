function [grp,fh]=usercluster(data,cg,cutoff,method,criterion,varargin)
%USERCLUSTER    Interactively cluster SEIZMO records
%
%    Usage:    grp=usercluster(data)
%              grp=usercluster(data,cg)
%              grp=usercluster(data,cg,cutoff)
%              grp=usercluster(data,cg,cutoff,method)
%              grp=usercluster(data,cg,cutoff,method,criterion)
%              grp=usercluster(data,cg,cutoff,method,criterion,...
%                              'field1',value1,...,'fieldN',valueN)
%              [grp,fh]=usercluster(...)
%
%    Description: GRP=USERCLUSTER(DATA) cross correlates the records in
%     DATA and passes the maximum correlation values to clustering routines
%     that facilitate rapid segregation of the dataset with just a few
%     mouse clicks.  GRP is a struct containing several fields providing
%     the parameters used in the analysis and some related info.  The
%     layout is as follows:
%       GRP.cutoff    = clustering limit
%       GRP.method    = linkage method
%       GRP.criterion = cluster formation method
%       GRP.perm      = permutation indices to match dendrogram ordering
%       GRP.color     = coloring for each cluster in the plot
%       GRP.T         = cluster index for each record
%
%     GRP=USERCLUSTER(DATA,CG) passes in the correlation matrix CG.  CG
%     should be N(N-1)/2 x 1 where N is the number of records in DATA.  CG
%     may be empty ([]) to trigger cross-correlation of records in DATA.
%
%     GRP=USERCLUSTER(DATA,CG,CUTOFF) sets the initial dissimilarity cutoff
%     for forming clusters to CUTOFF.  This should be a value between 0 and
%     1.  The default value is 0.2.
%
%     GRP=USERCLUSTER(DATA,CG,CUTOFF,METHOD) alters the method for forming
%     the dendrogram.  See the function LINKAGE for details.  The default
%     is 'average'.
%
%     GRP=USERCLUSTER(DATA,CG,CUTOFF,METHOD,CRITERION) sets the method for
%     forming clusters to CRITERION.  The default value is 'distance'.
%     Another valid value is 'inconsistent' although it is NOT RECOMMENDED
%     as it will be inconsistent with the dendrogram.  See CLUSTER for more
%     info.
%
%     GRP=USERCLUSTER(DATA,CG,CUTOFF,METHOD,CRITERION,'FIELD',VALUE,...)
%     passes field and value pairs on to PLOTDENDRO for further
%     customization.
%
%     [GRP,FH]=USERCLUSTER(...) returns the figure handle in FH.
%
%    Notes:
%
%    Examples:
%     Cluster starting with a dissimilarity cutoff of 0.05:
%      grp=usercluster(data,[],0.05);
%
%    See also: USERWINDOW, USERTAPER, USERALIGN, PLOTDENDRO,
%              CLUSTER, LINKAGE, PDIST, INCONSISTENT, DENDROGRAM

%     Version History:
%        Sep. 25, 2009 - rewrite and added documentation
%        Mar.  1, 2010 - minor doc update, better comments, updated for
%                        newer checking methods
%        Mar. 12, 2010 - pretty text menu for Octave
%        Mar. 17, 2010 - redesigned interface, several more options,
%                        grp.color lists group colors
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 17, 2010 at 16:15 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end
if(nargin>5 && ~mod(nargin,2))
    error('seizmo:usercluster:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
versioninfo(data,'dep');

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt clustering
try
    % correlate if no second input
    if(nargin==1 || isempty(cg))
        disp('CORRELATING DATASET (MAY TAKE A FEW MINUTES)')
        peaks=correlate(data,'npeaks',1);
        cg=peaks.cg;
    end

    % check correlation grid
    cg=squeeze(cg); sz=size(cg);
    if(sz(2)~=1 || numel(sz)>2)
        error('seizmo:usercluster:badCorrGrid',...
            'CG does not have proper dimensions!');
    end

    % default dendrogram cutoff
    if(nargin<3 || isempty(cutoff)); cutoff=0.2; end
    if(~isreal(cutoff) || ~isscalar(cutoff))
        error('seizmo:usercluster:badCutoff',...
            'CUTOFF must be a real-valued scalar in the range 0 to 1!');
    end
    if(cutoff<0); cutoff=0; end
    if(cutoff>1); cutoff=1; end

    % default/check method
    methods={'average' 'centroid' 'complete' ...
        'median' 'single' 'ward' 'weighted'};
    if(nargin<4 || isempty(method)); method='average'; end
    if(~any(strcmpi(method,methods)))
        error('seizmo:usercluster:badCriterion',...
            ['METHOD must be one of the following:\n' ...
            sprintf('''%s'' ',methods{:})]);
    end

    % default/check criterion
    criterions={'inconsistent' 'distance'};
    if(nargin<5 || isempty(criterion)); criterion='distance'; end
    if(~any(strcmpi(criterion,criterions)))
        error('seizmo:usercluster:badCriterion',...
            ['CRITERION must be one of the following:\n' ...
            sprintf('''%s'' ',criterions{:})]);
    end
    
    % save parameters
    grp.cutoff=cutoff;
    grp.method=method;
    grp.criterion=criterion;
    
    % outer loop - only breaks free on user command
    happy_user=false; fh=-1; cfh=-1; oldfh=[]; save=false;
    while(~happy_user)
        % list of extra plots
        oldfh=[oldfh(ishandle(oldfh)) cfh(ishandle(cfh)) fh(ishandle(fh))];
        
        % get choice from user
        choice=menu('CHANGE CLUSTERING SETTINGS?',...
            ['DISSIMILARITY CUTOFF (' num2str(grp.cutoff) ')'],...
            ['LINKAGE METHOD (' upper(grp.method) ')'],...
            ['CLUSTERING CRITERION (' upper(grp.criterion) ')'],...
            'CHECK LINKAGE FAITHFULNESS TO OBSERVATIONS',...
            'NO, CLUSTER AND RETURN INFO','CRASH!');
        
        % act on user choice
        switch choice
            case 1 % cutoff
                % get linkage
                Z=linkage(1-cg.',grp.method);
                
                % plot dendrogram/waveforms
                if(~save || ~ishandle(fh))
                    [grp.perm,grp.color,fh,sfh]=plotdendro(data,Z,...
                        varargin{:},'treelimit',grp.cutoff);
                end
                
                % menu telling user how to interactively adjust cutoff
                prompt={'+-------------------------------------------------------+'
                    '|               Welcome to SEIZMO''s interactive clustering function              |'
                    '+-------------------------------------------------------+'
                    '|                                                                                                               |'
                    '|                                            MOUSE USAGE                                             |'
                    '|                                                                                                               |'
                    '|    LEFT CLICK                      MIDDLE CLICK                      RIGHT CLICK   |'
                    '+-------------+          +-------------+          +--------------+'
                    '|     Set Cluster                      Finalize Limit                                              |'
                    '|         Limit                                                                                              |'
                    '+-------------------------------------------------------+'};
                % way cooler menu -- if only matlab gui's used fixed width
                if(strcmpi(getapplication,'OCTAVE'))
                    prompt={'+-------------------------------------------------------+'
                            '|  Welcome to SEIZMO''s interactive clustering function  |'
                            '+-------------------------------------------------------+'
                            '|                                                       |'
                            '|                     MOUSE USAGE:                      |'
                            '|                                                       |'
                            '|   LEFT CLICK        MIDDLE CLICK         RIGHT CLICK  |'
                            '+---------------+   +--------------+    +---------------+'
                            '|  Set Cluster       Finalize Limit                     |'
                            '|     Limit                                             |'
                            '+-------------------------------------------------------+'};
                end
                menu(prompt,'I''M READY!');
                
                % loop until right click
                button=1;
                while (button~=2)
                    % bring plot to focus (redraw if closed)
                    if(~ishandle(fh))
                        % redraw figure
                        [grp.perm,grp.color,fh,sfh]=plotdendro(data,Z,...
                            varargin{:},'treelimit',grp.cutoff);
                    else
                        % bring figure to focus
                        figure(fh);
                    end
                    subplot(sfh(1));
                    
                    % get user click/key
                    try
                        [x,y,button]=ginput(1);
                    catch
                        % user closed window - break from loop
                        button=2;
                    end

                    % action on left click
                    if (button==1)
                        grp.cutoff=x;

                        % redraw subplots (to show new grouping)
                        delete(sfh);
                        [grp.perm,grp.color,fh,sfh]=plotdendro(data,Z,...
                            varargin{:},'treelimit',grp.cutoff,...
                            'fighandle',fh);
                    end
                end
                
                % indicate plot is current
                save=true;
            case 2 % method
                % change method
                choice=menu('SELECT A LINKAGE METHOD',...
                    ['CURRENT (' upper(grp.method) ')'],...
                    methods{:});
                if(choice>1);
                    grp.method=methods{choice-1};
                    save=false;
                end
            case 3 % criterion
                % change criterion
                choice=menu('SELECT A CLUSTERING CRITERION',...
                    ['CURRENT (' upper(grp.method) ')'],...
                    criterions{:});
                if(choice>1)
                    grp.criterion=criterions{choice-1};
                    save=false;
                end
            case 4 % cophenetic inspection
                % get linkage
                Z=linkage(1-cg.',grp.method);
                
                % now get cophenetic stats
                [C,D]=cophenet(Z,1-cg.');
                
                % calculate Spearman's rank correlation
                rho=corr(1-cg,D','type','spearman');
                
                % plot against one another
                cfh=figure;
                plot([0 1],[0 1],'r','linewidth',3);
                hold on
                plot(1-cg.',D,'.');
                xlabel('TRUE DISTANCE (1 - Corr. Coeff.)');
                ylabel('COPHENETIC DISTANCE');
                title({['LINKAGE METHOD: ''' upper(grp.method) '''']
                       ['SPEARMAN''S RANK COEFFICIENT: ' num2str(rho)]});
            case 5 % cluster
                % get linkage
                Z=linkage(1-cg.',grp.method);
                
                % delete all old cluster plots & plot new if needed
                if(save)
                    close(oldfh(ishandle(oldfh) & oldfh~=fh));
                else
                    close(oldfh(ishandle(oldfh)));
                    [grp.perm,grp.color,fh]=plotdendro(data,Z,...
                        varargin{:},'treelimit',grp.cutoff);
                end
                
                
                % get clusters
                grp.T=cluster(Z,'cutoff',grp.cutoff,...
                    'criterion',grp.criterion);
                [idx,idx]=ismember(1:max(grp.T),grp.T(grp.perm));
                grp.color=grp.color(idx,:);
                
                % break loop
                happy_user=true;
            case 6 % crash
                error('seizmo:usertaper:killYourSelf',...
                    'User demanded Seppuku!');
        end
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
