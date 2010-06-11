function [data]=rotate_correlations(data)
%ROTATE_CORRELATIONS    Rotates EE/EN/NE/NN correlations to RR/RT/TR/TT
%
%    Usage:    data=rotate_correlations(data)
%
%    Description: DATA=ROTATE_CORRELATIONS(DATA) will rotate horizontal
%     correlogram sets in SEIZMO struct DATA to the radial and transverse
%     directions (defined by the azimuths between each pair of stations).
%     Note that correlograms are expected to be between Eest and North
%     components (make sure to rotate all records to this system prior to
%     correlation).  This is only compatible with correlograms generated
%     by CORRELATE (due to header field setup).  Also note that this will
%     only return correlogram sets for each unique station pair.
%
%    Notes:
%     - Correlogram .name fields are altered.
%
%    Header changes:
%     Master & Slave Fields may be switched.
%     KCMPNM & KT3 are changed to end with R or T.
%     USER0 & USER1 reflect the stream pair indices.
%
%    Examples:
%     
%
%    See also: CORRELATE, REVERSE_CORRELATIONS

%     Version History:
%        June 10, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 10, 2010 at 22:30 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check data
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_B','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'MULTIPLE_NPTS','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt rotation
try
    % verbosity
    verbose=seizmoverbose(false);
    
    % get component names, orientation code, drop verticals
    [mcmp,scmp]=getheader(data,'kt3','kcmpnm');
    moc=char(mcmp); moc=cellstr(moc(:,3));
    soc=char(scmp); soc=cellstr(soc(:,3));
    vert=strcmpi('Z',moc) | strcmpi('Z',soc);
    data(vert)=[]; moc(vert)=[]; soc(vert)=[];
    
    % skip if nothing left
    if(~numel(data))
        % toggle checking back
        seizmocheck_state(oldseizmocheckstate);
        checkheader_state(oldcheckheaderstate);
        return;
    end
    
    % make sure all remaining orientation codes are E/N
    if(~all(strcmpi([moc; soc],'E') | strcmpi([moc; soc],'N')))
        error('seizmo:rotate_correlations:NEonly',...
            'Horizontal correlograms must be EE/EN/NE/NN only!');
    end
    
    % get station indices
    % - well really we just get stream names and do a unique sort
    [ksname,kt0,kt1,kt2,kt3]=getheader(data,...
        'kname','kt0','kt1','kt2','kt3');
    kmname=[kt0 kt1 kt2 kt3];
    kmname=strnlen(kmname,8); % force 8 character fields
    ksname=strnlen(ksname,8); % force 8 character fields
    kstream=[kmname; ksname];
    kstream(:,4)=strnlen(kstream(:,4),2);
    [idx,idx,idx]=unique(strcat(kstream(:,1),kstream(:,2),...
        kstream(:,3),kstream(:,4)),'first');
    mi=idx(1:end/2);
    si=idx(end/2+1:end);
    
    % drop autocorrelations
    ac=mi==si;
    data(ac)=[];
    
    % re-get station indices
    kmname(ac,:)=[];
    ksname(ac,:)=[];
    kstream=[kmname; ksname];
    kstream(:,4)=strnlen(kstream(:,4),2);
    [idx,idx,idx]=unique(strcat(kstream(:,1),kstream(:,2),...
        kstream(:,3),kstream(:,4)),'first');
    mi=idx(1:end/2);
    si=idx(end/2+1:end);
    
    % reverse correlations where master<slave
    rc=mi>si;
    data(rc)=reverse_correlations(data(rc));
    [kmname(rc,:),ksname(rc,:)]=deal(ksname(rc,:),kmname(rc,:));
    
    % remove repeats
    [idx,idx]=unique(strcat(kmname(:,1),kmname(:,2),kmname(:,3),...
        kmname(:,4),ksname(:,1),ksname(:,2),ksname(:,3),ksname(:,4)),...
        'first');
    kmname=kmname(idx,:);
    ksname=ksname(idx,:);
    data=data(idx);
    
    % remove incomplete sets
    k2=[kmname ksname];
    k2(:,[4 8])=strnlen(k2(:,[4 8]),2);
    [idx,idx,idx]=unique(strcat(k2(:,1),k2(:,2),k2(:,3),k2(:,4),k2(:,5),...
        k2(:,6),k2(:,7),k2(:,8)),'first');
    bad=false(numel(idx));
    for i=1:max(idx)
        if(sum(idx==i)~=4)
            bad(idx==i)=true;
        end
    end
    data(bad)=[];
    kmname(bad,:)=[];
    ksname(bad,:)=[];
    
    % reorder one last time (to quads of EE, EN, NE, NN)
    k2=[kmname ksname];
    k2(:,[4 8])=strnlen(k2(:,[4 8]),2);
    moc=char(kmname(:,4)); moc=cellstr(moc(:,3));
    soc=char(ksname(:,4)); soc=cellstr(soc(:,3));
    k2=[k2 moc soc];
    [idx,idx]=sort(strcat(k2(:,1),k2(:,2),k2(:,3),k2(:,4),k2(:,5),...
        k2(:,6),k2(:,7),k2(:,8),k2(:,9),k2(:,10)));
    k2=k2(idx,:);
    data=data(idx);
    kmname=kmname(idx,:);
    ksname=ksname(idx,:);
    
    % get master/slave indices one last time
    kstream=[kmname; ksname];
    kstream(:,4)=strnlen(kstream(:,4),2);
    [idx,idx,idx]=unique(strcat(kstream(:,1),kstream(:,2),...
        kstream(:,3),kstream(:,4)),'first');
    mi=idx(1:end/2);
    si=idx(end/2+1:end);
    
    % get azis
    [az,baz]=getheader(data,'az','baz');
    az=az*pi/180;
    baz=baz*pi/180;
    
    % combine quads (memory expensive)
    horz=joinrecords(data(1:4:end),data(2:4:end),...
                     data(3:4:end),data(4:4:end));
    
    % loop over quads, rotating to RR/RT/TR/TT
    % - see Lin et al 2008
    % - note NE/NN order is switched
    nquads=numel(horz);
    if(verbose)
        disp('Rotating Horizontal Correlogram Set(s)');
        print_time_left(0,nquads);
    end
    for i=1:nquads
        data((i-1)*4+1)=seizmofun(horz(i),...
            @(x)(-sin(az(i))*sin(baz(i)).*x(:,1)...
            -sin(az(i))*cos(baz(i)).*x(:,2)...
            -cos(az(i))*sin(baz(i)).*x(:,3)...
            -cos(az(i))*cos(baz(i)).*x(:,4)));
        data((i-1)*4+2)=seizmofun(horz(i),...
            @(x)(-sin(az(i))*cos(baz(i)).*x(:,1)...
            +sin(az(i))*sin(baz(i)).*x(:,2)...
            -cos(az(i))*cos(baz(i)).*x(:,3)...
            +cos(az(i))*sin(baz(i)).*x(:,4)));
        data((i-1)*4+3)=seizmofun(horz(i),...
            @(x)(-cos(az(i))*sin(baz(i)).*x(:,1)...
            -cos(az(i))*cos(baz(i)).*x(:,2)...
            +sin(az(i))*sin(baz(i)).*x(:,3)...
            +sin(az(i))*cos(baz(i)).*x(:,4)));
        data((i-1)*4+4)=seizmofun(horz(i),...
            @(x)(-cos(az(i))*sin(baz(i)).*x(:,1)...
            +cos(az(i))*sin(baz(i)).*x(:,2)...
            +sin(az(i))*cos(baz(i)).*x(:,3)...
            -sin(az(i))*sin(baz(i)).*x(:,4)));
        if(verbose); print_time_left(i,nquads); end
    end
    
    % fix headers
    kmcmp=k2(:,4);
    kmcmp([1:4:end 2:4:end])=strcat(kmcmp([1:4:end 2:4:end]),'R');
    kmcmp([3:4:end 4:4:end])=strcat(kmcmp([3:4:end 4:4:end]),'T');
    kscmp=k2(:,8);
    kscmp([1:4:end 3:4:end])=strcat(kscmp([1:4:end 3:4:end]),'R');
    kscmp([2:4:end 4:4:end])=strcat(kscmp([2:4:end 4:4:end]),'T');
    data=changeheader(data,'kt3',kmcmp,'kcmpnm',kscmp,...
        'user0',mi,'user1',si);
    
    % fix names
    k2=strtrim(k2);
    name=strcat('CORR_-_ROTATED_-_MASTER_-_',k2(:,1),'.',k2(:,2),'.',...
        k2(:,3),'.',kmcmp,'_-_SLAVE_-_',k2(:,5),'.',k2(:,6),'.',k2(:,7),...
        '.',kscmp);
    [data.name]=deal(name{:});
    
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
