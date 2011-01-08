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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 17, 2010 at 10:00 GMT

% todo:

% check nargin
error(nargchk(4,4,nargin));

% check data structure
versioninfo(data,'dep');
data=checkheader(data);

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% present interface
try
    % number of records
    nrecs=numel(data);
    
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
    
    % expand scalars
    if(isscalar(arr)); arr(1:nrecs,1)=arr; end
    if(isscalar(pol)); pol(1:nrecs,1)=pol; end
    
    % outer loop
    happy_user=false;
    while(~happy_user)
        % get choice from user
        choice=menu('ADJUST/COMBINE CLUSTERS OPTIONS',...
            'COMBINE CLUSTERS',...
            'SPLIT CLUSTERS BY RECORD POLARITIES',...
            'ADJUST GROUND UNITS',...
            'FLIP CLUSTER POLARITY',...
            'TIME SHIFT CLUSTER',...
            'MANUALLY PICK GOOD CLUSTERS',...
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
                    close(get(tmpax(1),'parent'));
                    cgrps=find(tmp.good);
                    
                    % who is in these groups
                    crecs=ismember(grp.T,cgrps);
                    
                    % combined colormap
                    ccmap=grp.color(grp.T(crecs),:);
                    
                    % show a combined plot
                    tmpax=plot2(data(crecs),'cmap',ccmap);
                    
                    % ask if that is ok or not (cannot undo!)
                    choice=menu({'Combine Clusters?' ...
                        'Note 1: Cannot be Undone!' ...
                        'Note 2: Undoes manually selected good clusters!'},...
                        'YES','NO','RETRY');
                    
                    % close tmp plot
                    close(get(tmpax,'parent'));
                    
                    % action
                    switch choice
                        case 1 % combine
                            % change grp.T, grp.pop to all lowest
                            % - do we shift numbers? no
                            % - do we lose manual grp selection? yes
                            grp.T(crecs)=min(cgrps);
                            grp.pop=histc(grp.T,1:max(grp.T));
                            grp.good=grp.pop>=grp.popcut;
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
                close(get(tmpax(1),'parent'));
            case 4 % flip
                % let user select clusters
                tmp=grp;
                tmp.good=false;
                [tmp,tmpax]=selectclusters(data,tmp);
                close(get(tmpax(1),'parent'));
                cgrps=find(tmp.good);
                
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
                close(get(tmpax(1),'parent'));
                
                % apply time shift
                arr(crecs)=arr(crecs)+shift;
                data(crecs)=timeshift(data(crecs),shift);
            case 6 % select clusters manually
                % cut clusters manually
                [grp,tmpax]=selectclusters(data,grp);
                close(get(tmpax(1),'parent'));
            case 7 % break loop
                happy_user=true;
        end
    end
    
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
