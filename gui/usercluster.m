function [grp,oldax]=usercluster(data,cg,distcut,method,crit,pcut,varargin)
%USERCLUSTER    Interactively cluster SEIZMO records
%
%    Usage:    grp=usercluster(data)
%              grp=usercluster(data,cg)
%              grp=usercluster(data,cg,distcut)
%              grp=usercluster(data,cg,distcut,method)
%              grp=usercluster(data,cg,distcut,method,criterion)
%              grp=usercluster(data,cg,distcut,method,criterion,popcut)
%              grp=usercluster(data,cg,distcut,method,criterion,popcut,...
%                              'field1',value1,...,'fieldN',valueN)
%              [grp,ax]=usercluster(...)
%
%    Description:
%     GRP=USERCLUSTER(DATA) cross correlates the records in DATA and passes
%     the maximum correlation values to clustering routines that facilitate
%     rapid segregation of the dataset with just a few mouse clicks.  GRP
%     is a struct containing several fields providing the parameters used
%     in the analysis and some related info.  The layout is as follows:
%       GRP.distcut   = clustering distance limit
%       GRP.method    = linkage method
%       GRP.criterion = cluster formation method
%       GRP.pop       = cluster population
%       GRP.popcut    = "good" cluster population threshhold
%       GRP.good      = "good" clusters (logical array)
%       GRP.perm      = permutation indices to match dendrogram ordering
%       GRP.color     = coloring for each cluster in the plot
%       GRP.T         = cluster index for each record
%
%     GRP=USERCLUSTER(DATA,CG) passes in the correlation matrix CG.  CG
%     should be N(N-1)/2 x 1 where N is the number of records in DATA.  CG
%     may be empty ([]) to trigger cross-correlation of records in DATA.
%
%     GRP=USERCLUSTER(DATA,CG,DISTCUT) sets the initial dissimilarity cut
%     for forming clusters to DISTCUT.  This should be a value between 0
%     and 1.  The default value is 0.2.
%
%     GRP=USERCLUSTER(DATA,CG,DISTCUT,METHOD) alters the method for forming
%     the dendrogram.  See the function LINKAGE for details.  The default
%     is 'average'.
%
%     GRP=USERCLUSTER(DATA,CG,DISTCUT,METHOD,CRITERION) sets the method for
%     forming clusters to CRITERION.  The default value is 'distance'.
%     Another valid value is 'inconsistent' although it is NOT RECOMMENDED
%     as it will be inconsistent with the dendrogram.  See CLUSTER for more
%     info.
%
%     GRP=USERCLUSTER(DATA,CG,DISTCUT,METHOD,CRITERION,POPCUT) sets the
%     default population cutoff.  The default is set at 3 otherwise. 
%
%     GRP=USERCLUSTER(DATA,CG,DISTCUT,METHOD,CRITERION,POPCUT,...
%     'FIELD1',VALUE1,...,'FIELDN',VALUEN) passes field and value pairs on
%     to PLOTDENDRO for further customization.
%
%     [GRP,AX]=USERCLUSTER(...) returns the axes handles in AX.
%
%    Notes:
%
%    Examples:
%     % Cluster starting with a dissimilarity cutoff of 0.05:
%     grp=usercluster(data,[],0.05);
%
%    See also: USERWINDOW, USERTAPER, USERALIGN, PLOTDENDRO, CORRELATE
%              CLUSTER, LINKAGE, PDIST, INCONSISTENT, DENDROGRAM
%              POPCUT, PLOTPOP, PLOTCLUSTERS, SELECTCLUSTERS

%     Version History:
%        Sep. 25, 2009 - rewrite and added documentation
%        Mar.  1, 2010 - minor doc update, better comments, updated for
%                        newer checking methods
%        Mar. 12, 2010 - pretty text menu for Octave
%        Mar. 17, 2010 - redesigned interface, several more options,
%                        grp.color lists group colors
%        Mar. 22, 2010 - make sure input CG sizes up
%        Mar. 24, 2010 - minor whitespace fix
%        Apr. 21, 2010 - replace crash button with exit (but still crash)
%        May   7, 2010 - button to draw cluster map
%        June 26, 2010 - mapclusters deprecated (uses mapstations)
%        Aug. 26, 2010 - update for axes plotting output, checkheader fix
%        Sep. 21, 2010 - added POPCUT & SELECTRECORDS into the fold,
%                        changed CUTOFF to DISTCUT
%        Sep. 22, 2010 - a host of bug fixes
%        Sep. 30, 2010 - added pop field, fixed nasty bug
%        Oct.  1, 2010 - fixed no left click distcut bug
%        Oct.  2, 2010 - return all axes handles
%        Oct.  6, 2010 - combine clusters option
%        Oct. 10, 2010 - combine clusters now in ADJUSTCLUSTERS, drop exit
%                        and crash option, event grid for map
%        Jan.  6, 2011 - use key2zoompan
%        Jan. 13, 2011 - finally fixed cluster map warning issue
%        Feb. 25, 2011 - fixed bug in criterion selection menu
%        Mar. 31, 2011 - prompt once quick fix
%        Apr.  3, 2012 - minor doc update, use seizmocheck
%        Jan. 30, 2013 - update for new correlate
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 30, 2013 at 10:00 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>6 && mod(nargin,2))
    error('seizmo:usercluster:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data,'dep'));

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);

% prompt state variable is persistent
persistent showprompt
if(isempty(showprompt)); showprompt=true; end

% attempt clustering
try
    % correlate if no second input
    if(nargin==1 || isempty(cg))
        disp('CORRELATING DATASET (MAY TAKE A FEW MINUTES)')
        peaks=correlate(data,'mcxc','noauto','normxc','peaks');
        cg=peaks.cg;
    end
    
    % number of records
    nrecs=numel(data);

    % check correlation grid
    cg=squeeze(cg); sz=size(cg);
    if(sz(2)~=1 || numel(sz)>2 || (nrecs^2-nrecs)/2~=sz(1))
        error('seizmo:usercluster:badCorrGrid',...
            'CG does not have proper dimensions!');
    end

    % default dendrogram dissimilarity cutoff
    if(nargin<3 || isempty(distcut)); distcut=0.2; end
    if(~isreal(distcut) || ~isscalar(distcut))
        error('seizmo:usercluster:badCutoff',...
            'DISTCUT must be a real-valued scalar in the range 0 to 1!');
    end
    if(distcut<0); distcut=0; end
    if(distcut>1); distcut=1; end

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
    if(nargin<5 || isempty(crit)); crit='distance'; end
    if(~any(strcmpi(crit,criterions)))
        error('seizmo:usercluster:badCriterion',...
            ['CRITERION must be one of the following:\n' ...
            sprintf('''%s'' ',criterions{:})]);
    end
    
    % default/check population cutoff
    if(nargin<6 || isempty(pcut)); pcut=3; end
    if(~isreal(pcut) || ~isscalar(pcut) || pcut<0)
        error('seizmo:usercluster:badPopulationCutoff',...
            'POPCUT must be a positive real-valued scalar!');
    end
    
    % event location
    [evla,evlo]=getheader(data,'evla','evlo');
    ev=unique([evla evlo],'rows');
    
    % save parameters
    grp.distcut=distcut;
    grp.method=method;
    grp.criterion=crit;
    grp.pop=[];
    grp.popcut=pcut;
    grp.good=[];
    
    % outer loop - only breaks free on user command
    happy_user=false;
    ax=[-1 -1]; % dendrogram & waveform axes
    cx=-1;      % cophenetic distance vs actual distance
    mx=-1;      % cluster map
    px=[];      % popcut axes
    sx=[];      % selectclusters axes
    oldax=[];   % axes from previous loop
    save=false;
    while(~happy_user)
        % get linkage
        Z=linkage(1-cg.',grp.method);
        
        % delete all old cluster plots & plot dendrogram if needed
        if(~save)
            oldfh=get(oldax(ishandle(oldax)),'parent');
            if(iscell(oldfh)); oldfh=cell2mat(oldfh); end
            close(oldfh);
            [grp.perm,grp.color,ax]=plotdendro(data,Z,...
                varargin{:},'distcut',grp.distcut);
            
            % get clusters
            grp.T=cluster(Z,'cutoff',grp.distcut,...
                'criterion',grp.criterion);
            
            % fix color matrix
            [idx,idx]=ismember(1:max(grp.T),grp.T(grp.perm));
            grp.color=grp.color(idx,:);
            
            % get good clusters
            grp.pop=histc(grp.T,1:max(grp.T));
            grp.good=grp.pop>=grp.popcut;
        end
        
        % get choice from user
        choice=menu('CLUSTERING OPTIONS',...
            ['DISSIMILARITY CUTOFF (' num2str(grp.distcut) ')'],...
            ['POPULATION CUTOFF (' num2str(grp.popcut) ')'],...
            'MANUALLY PICK GOOD CLUSTERS',...
            ['LINKAGE METHOD (' upper(grp.method) ')'],...
            ['CLUSTERING CRITERION (' upper(grp.criterion) ')'],...
            'CHECK LINKAGE FAITHFULNESS',...
            'VIEW CLUSTER MAP',...
            'RETURN');
        
        % act on user choice
        switch choice
            case 1 % distcut
                % plot dendrogram/waveforms
                if(any(~ishandle(ax)))
                    [grp.perm,grp.color,ax]=plotdendro(data,Z,...
                        varargin{:},'distcut',grp.distcut);
                else
                    % redraw subplots (to get perm & color)
                    [grp.perm,grp.color,ax]=plotdendro(data,Z,...
                        varargin{:},'distcut',grp.distcut,'ax',ax);
                end
                
                % menu telling user how to interactively adjust distcut
                if(showprompt)
                    prompt={'+-------------------------------------------------------+'
                        '|               Welcome to SEIZMO''s interactive clustering function              |'
                        '+-------------------------------------------------------+'
                        '|                                                                                                               |'
                        '|                                            MOUSE USAGE                                             |'
                        '|                                                                                                               |'
                        '|    LEFT CLICK                      MIDDLE CLICK                      RIGHT CLICK   |'
                        '+-------------+          +-------------+          +--------------+'
                        '|     Set Distance                    Finalize Limit                                              |'
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
                            '|  Set Distance       Finalize Limit                    |'
                            '|     Limit                                             |'
                            '+-------------------------------------------------------+'};
                    end
                    menu(prompt,'I''M READY!');
                    showprompt=false;
                end
                
                % loop until right click
                button=1;
                while (button~=2)
                    % bring plot to focus (redraw if closed)
                    if(any(~ishandle(ax)))
                        % redraw figure
                        [grp.perm,grp.color,ax]=plotdendro(data,Z,...
                            varargin{:},'distcut',grp.distcut);
                    end
                    axes(ax(1));
                    
                    % get user click/key
                    try
                        [x,y,button]=ginput(1);
                    catch
                        % user closed window - break from loop
                        break;
                    end
                    
                    % skip if not dendrogram axis
                    if(ax(1)~=gca); continue; end

                    % action on left click
                    if (button==1)
                        grp.distcut=x;

                        % redraw subplots (to show new grouping)
                        [grp.perm,grp.color,ax]=plotdendro(data,Z,...
                            varargin{:},'distcut',grp.distcut,'ax',ax);
                    else
                        key2zoompan(button,ax(1));
                    end
                end
            
                % get clusters
                grp.T=cluster(Z,'cutoff',grp.distcut,...
                    'criterion',grp.criterion);
                
                % fix color matrix
                [idx,idx]=ismember(1:max(grp.T),grp.T(grp.perm));
                grp.color=grp.color(idx,:);
            
                % get good clusters
                grp.pop=histc(grp.T,1:max(grp.T));
                grp.good=grp.pop>=grp.popcut;
                
                % indicate plot is current
                save=true;
            case 2 % popcut
                % cut clusters by population
                if(any(ishandle(px)))
                    px=px(ishandle(px));
                    close(get(px(1),'parent'));
                end
                [grp,px]=popcut(data,grp);
                save=true;
            case 3 % select clusters manually
                % cut clusters manually
                if(any(ishandle(sx)))
                    sx=sx(ishandle(sx));
                    close(get(sx(1),'parent'));
                end
                [grp,sx]=selectclusters(data,grp);
                sx=sx';
                save=true;
            case 4 % method
                % change method
                choice=menu('SELECT A LINKAGE METHOD',...
                    ['CURRENT (' upper(grp.method) ')'],...
                    methods{:});
                if(choice>1);
                    grp.method=methods{choice-1};
                    save=false;
                end
            case 5 % criterion
                % change criterion
                choice=menu('SELECT A CLUSTERING CRITERION',...
                    ['CURRENT (' upper(grp.criterion) ')'],...
                    criterions{:});
                if(choice>1)
                    grp.criterion=criterions{choice-1};
                    save=false;
                end
            case 6 % cophenetic inspection
                % get cophenetic stats
                [C,D]=cophenet(Z,1-cg.');
                
                % calculate Spearman's rank correlation
                rho=corr(1-cg,D','type','spearman');
                
                % plot against one another
                if(~ishandle(cx))
                    cfh=figure('color','k');
                    cx=axes('parent',cfh);
                end
                plot(cx,[0 1],[0 1],'r','linewidth',3);
                hold(cx,'on');
                plot(cx,1-cg.',D,'.');
                hold(cx,'off');
                set(cx,'color','k','xcolor','w','ycolor','w');
                xlabel(cx,'TRUE DISTANCE (1 - Corr. Coeff.)');
                ylabel(cx,'COPHENETIC DISTANCE');
                title(cx,{['LINKAGE METHOD: ''' upper(grp.method) '''']
                    ['SPEARMAN''S RANK COEFFICIENT: ' num2str(rho)]},...
                    'color','w');
                save=true;
            case 7 % map
                % map clusters
                if(~ishandle(mx))
                    mfh=figure('color','k');
                    mx=axes('parent',mfh);
                end
                mx=mapstations(data,'ax',mx);
                if(size(ev,1)==1 && ~any(isnan(ev)))
                    mapeventgrid(mx,ev(1),ev(2));
                end
                set(findobj(mx,'tag','stations'),...
                    'cdata',grp.color(grp.T,:));
                drawnow;
                movekids(findobj(mx,'tag','stations'),'front');
                movekids(findobj(mx,'tag','events'),'front');
                save=true;
            case 8 % break loop
                happy_user=true;
        end
        
        % list of valid plots from previous run
        oldax=[oldax(ishandle(oldax)) cx(ishandle(cx)) ...
            px(ishandle(px)) sx(ishandle(sx)) ...
            mx(ishandle(mx)) ax(ishandle(ax))];
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
