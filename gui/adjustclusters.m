function [data,grp,arr,pol]=adjustclusters(data,grp,arr,pol)
%ADJUSTCLUSTERS    Combine & Adjust timing/polarity of clusters graphically
%
%    Usage:    [data,grp,arr,pol]=adjustclusters(data,grp,arr,pol)
%
%    Description:
%     [DATA,GRP,ARR,POL]=ADJUSTCLUSTERS(DATA,GRP,ARR,POL) presents an
%     interface for adjusting and combining clusters given by SEIZMO struct
%     DATA and cluster struct GRP (see USERCLUSTER for more details).  ARR
%     & POL are vectors equal in size to DATA that contain the timing and
%     polarity of the records.
%
%    Notes:
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 10, 2010 at 10:00 GMT

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
            case 2 % flip
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
            case 3 % timeshift
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
            case 4 % select clusters manually
                % cut clusters manually
                if(any(ishandle(sx)))
                    sx=sx(ishandle(sx));
                    close(get(sx(1),'parent'));
                end
                [grp,tmpax]=selectclusters(data,grp);
                close(get(tmpax(1),'parent'));
            case 5 % break loop
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
