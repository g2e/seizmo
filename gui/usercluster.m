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
%       GRP.method    = cluster distance method
%       GRP.criterion = cluster formation method
%       GRP.perm      = permutation indices for matching plot ordering
%       GRP.color     = record coloring from plot
%       GRP.T         = group index for each record
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
%     Another valid value is 'inconsistent'.  See CLUSTER for more info.
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2010 at 01:05 GMT

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
    if(nargin<4 || isempty(method)); method='average'; end

    % default/check criterion
    if(nargin<5 || isempty(criterion)); criterion='distance'; end
    if(~any(strcmpi(criterion,{'inconsistent' 'distance'})))
        error('seizmo:usercluster:badCriterion',...
            'CRITERION must be ''inconsistent'' or ''distance''!');
    end

    % get linkage
    Z=linkage(1-cg.',method);

    % save parameters
    grp.cutoff=cutoff;
    grp.method=method;
    grp.criterion=criterion;

    % cluster analysis loop
    disp('WAITING FOR USER INPUT OF CLUSTERING CUTOFF')
    while (1)
        % dendrogram/waveforms
        [grp.perm,grp.color,fh,sfh]=plotdendro(data,Z,varargin{:},...
            'treelimit',grp.cutoff);
        drawnow;
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
        %{
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
        %}
        menu(prompt,'OK');

        % loop until right click
        button=1;
        while (button~=2)
            % user picks
            figure(fh);
            subplot(sfh(1));
            [x,y,button]=ginput(1);

            % action on left click
            if (button==1)
                grp.cutoff=x;

                % redraw plot (to show new grouping)
                delete(sfh);
                [grp.perm,grp.color,fh,sfh]=plotdendro(data,Z,varargin{:},...
                    'treelimit',grp.cutoff,'fighandle',fh);
            end
        end

        % check satisfaction
        done=menu('KEEP GROUPS?','YES','NO - TRY AGAIN','NO - CRASH!');
        switch done
            case 1 % rainbow's end
                disp(['DISSIMILARITY LIMIT SET TO: ' num2str(grp.cutoff)])
                grp.T=cluster(Z,'cutoff',grp.cutoff,'criterion',criterion);
                seizmocheck_state(oldseizmocheckstate);
                return;
            case 2 % never never quit
                close(fh);
            case 3 % i bear too great a shame to go on
                error('seizmo:usercluster:killYourSelf',...
                    'User Demanded Seppuku!')
        end
    end
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
