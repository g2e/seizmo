function [data,grp,arr,pol,units]=adjustclusters(data,grp,arr,pol)
%ADJUSTCLUSTERS    Combine & Adjust timing/polarity of clusters graphically
%
%    Usage:    [data,grp,arr,pol,units]=adjustclusters(data,grp,arr,pol)
%
%    Description:
%     [DATA,GRP,ARR,POL,UNITS]=ADJUSTCLUSTERS(DATA,GRP,ARR,POL) presents an
%     interface for adjusting, combining, and splitting clusters given by
%     SEIZMO struct DATA and cluster struct GRP (see USERCLUSTER for more
%     details).  ARR & POL are vectors equal in size to DATA that contain
%     the timing and polarity of the records.  UNITS is a NRECSx1 matrix of
%     integer values ranging from -2 to 2 that indicate how many times the
%     data was/needs to be integrated -- negative numbers indicate
%     differentiation.
%
%    Notes:
%     - Combining/Splitting clusters will reset manual "good" cluster
%       selection to those above the population cutoff given by GRP.POPCUT.
%
%    Examples:
%     % Use ADJUSTCLUSTERS after USERCLUSTER to allow for advanced
%     % manipulation of clustered waveforms:
%     grp=usercluster(data,xc);
%     [data,grp,arr,pol]=adjustclusters(data,grp,0,1);
%
%    See also: USERCLUSTER, USERALIGN, SELECTCLUSTERS, TIMESHIFT, MULTIPLY

%     Version History:
%        Oct. 10, 2010 - initial version
%        Nov. 20, 2010 - fix manual picking bug
%        Dec. 12, 2010 - added split clusters by polarity, add Note
%        Dec. 17, 2010 - added groundunits code calls
%        Jan. 13, 2011 - added in map clusters option, fixed many bugs
%        Jan. 14, 2011 - fixed unassigned units bug
%        Jan. 16, 2011 - fix combine clusters not removing 0 pop clusters
%                        that are last in the set (ie make sure max(grp.T)
%                        points to a non-zero group), close cluster map
%        Mar. 31, 2011 - fix plot closing breakage when only 1 cluster
%        Apr.  3, 2012 - use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 10:00 GMT

% todo:

% check nargin
error(nargchk(4,4,nargin));

% check data structure
error(seizmocheck(data,'dep'));
data=checkheader(data);

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% present interface
try
    % number of records
    nrecs=numel(data);
    
    % default units
    units=zeros(nrecs,1);
    
    % check grp struct
    if(~isstruct(grp) ...
            || any(~ismember({'T' 'color' 'good' 'popcut'},...
            fieldnames(grp))) || numel(grp.T)~=nrecs ...
            || ~isequal([max(grp.T) 3],size(grp.color)))
        error('seizmo:selectclusters:badInput',...
            'GRP must be a struct with fields T, color, good & popcut!');
    end
    
    % check arr & pol
    if(~isreal(arr) || ~any(numel(arr)==[1 nrecs]))
        error('seizmo:adjustclusters:badInput',...
            'ARR must be a real-valued vector in seconds!');
    elseif(~isreal(pol) || ~any(numel(pol)==[1 nrecs]) || any(abs(pol)~=1))
        error('seizmo:adjustclusters:badInput',...
            'POL must be a vector of 1s & -1s!');
    end
    
    % event location
    [evla,evlo]=getheader(data,'evla','evlo');
    ev=unique([evla evlo],'rows');
    
    % expand scalars
    if(isscalar(arr)); arr(1:nrecs,1)=arr; end
    if(isscalar(pol)); pol(1:nrecs,1)=pol; end
    
    % outer loop
    happy_user=false; mx=-1;
    while(~happy_user)
        % get choice from user
        choice=menu('ADJUST/COMBINE CLUSTERS OPTIONS',...
            'COMBINE CLUSTERS',...
            'SPLIT CLUSTERS BY RECORD POLARITIES',...
            'ADJUST GROUND UNITS',...
            'FLIP CLUSTER POLARITY',...
            'TIME SHIFT CLUSTER',...
            'MANUALLY PICK GOOD CLUSTERS',...
            'MAP CLUSTERS',...
            'RETURN');
        
        % act on choice
        switch choice
            case 1 % combine
                % loop until user says
                happy=false;
                while(~happy)
                    % let user select clusters
                    tmp=grp;
                    tmp.good=false;
                    [tmp,tmpax]=selectclusters(data,tmp);
                    if(any(ishandle(tmpax)))
                        fh=get(tmpax(ishandle(tmpax)),'parent');
                        if(iscell(fh)); fh=cell2mat(fh); end
                        close(unique(fh));
                    end
                    cgrps=find(tmp.good);
                    
                    % skip if <2 selected
                    if(numel(cgrps)<2); happy=true; continue; end
                    
                    % who is in these groups
                    crecs=ismember(grp.T,cgrps);
                    
                    % combined colormap
                    ccmap=grp.color(grp.T(crecs),:);
                    
                    % show a combined plot
                    tmpax=plot2(data(crecs),'cmap',ccmap);
                    
                    % ask if that is ok or not (cannot undo!)
                    choice=menu({'Combine Clusters?' ...
                        'Note 1: Cannot be Undone!' ...
                        'Note 2: Resets good cluster selection!'},...
                        'YES','NO','RETRY');
                    
                    % close tmp plot
                    close(get(tmpax,'parent'));
                    
                    % action
                    switch choice
                        case 1 % combine
                            % change grp.T to lowest 2+ group
                            % - if no 2+ then lowest group
                            % - do we shift numbers? no
                            % - truncate off 0 pop groups on end
                            % - do we lose manual grp selection? yes
                            mingrp=min(cgrps(grp.pop(cgrps)>1));
                            if(isempty(mingrp))
                                mingrp=min(cgrps);
                                % make group randomly whiter
                                grp.color(mingrp,:)=...
                                    (1+rand).*grp.color(mingrp,:);
                            end
                            grp.T(crecs)=mingrp;
                            grp.pop=histc(grp.T,1:max(grp.T));
                            grp.good=grp.pop>=grp.popcut;
                            grp.color((max(grp.T)+1):end,:)=[];
                            happy=true;
                        case 2 % don't
                            happy=true;
                    end
                end
            case 2 % polarity split
                grp=split_clusters_by_polarity(grp,pol);
            case 3 % ground units
                [data,grp,units,tmpax]=...
                    selectclusters_and_groundunits(data,grp);
                if(any(ishandle(tmpax)))
                    close(unique(cell2mat(get(...
                        tmpax(ishandle(tmpax)),'parent'))));
                end
            case 4 % flip
                % let user select clusters
                tmp=grp;
                tmp.good=false;
                [tmp,tmpax]=selectclusters(data,tmp);
                if(any(ishandle(tmpax)))
                    fh=get(tmpax(ishandle(tmpax)),'parent');
                    if(iscell(fh)); fh=cell2mat(fh); end
                    close(unique(fh));
                end
                cgrps=find(tmp.good);
                
                % skip if none selected
                if(isempty(cgrps)); continue; end
                
                % who is in these groups
                crecs=ismember(grp.T,cgrps);
                
                % flip polarity
                pol(crecs)=pol(crecs).*(-1);
                data(crecs)=multiply(data(crecs),-1);
            case 5 % timeshift
                % let user select clusters
                tmp=grp;
                tmp.good=false;
                [tmp,tmpax]=selectclusters(data,tmp);
                cgrps=find(tmp.good);
                
                % skip if none selected
                if(isempty(cgrps)); continue; end
                
                % who is in these groups
                crecs=ismember(grp.T,cgrps);
                
                % get time shift
                shift=0;
                tmp=inputdlg(...
                    'Shift time of clusters by how much (in sec)? [0]:',...
                    'Time Shift Clusters',1,{'0'});
                if(~isempty(tmp))
                    try
                        tmp=str2double(tmp{:});
                        if(isscalar(tmp) && isreal(tmp))
                            shift=tmp;
                        end
                    catch
                        % do not change
                    end
                end
                if(any(ishandle(tmpax)))
                    fh=get(tmpax(ishandle(tmpax)),'parent');
                    if(iscell(fh)); fh=cell2mat(fh); end
                    close(unique(fh));
                end
                
                % apply time shift
                arr(crecs)=arr(crecs)+shift;
                data(crecs)=timeshift(data(crecs),shift);
            case 6 % select clusters manually
                % cut clusters manually
                [grp,tmpax]=selectclusters(data,grp);
                if(any(ishandle(tmpax)))
                    fh=get(tmpax(ishandle(tmpax)),'parent');
                    if(iscell(fh)); fh=cell2mat(fh); end
                    close(unique(fh));
                end
            case 7 % map clusters
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
            case 8 % break loop
                happy_user=true;
        end
    end
    
    % close map figure
    if(ishandle(mx)); close(get(mx,'parent')); end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end
