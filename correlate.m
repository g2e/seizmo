function [data]=correlate(data1,varargin)
%CORRELATE    Compute cross correlations for SAClab data records
%
%    Description: CORRELATE(DATA) cross correlates every possible pairing
%     of records in DATA.  Output is a SAClab dataset of the correlograms.
%     
%     CORRELATE(DATA1,DATA2) cross correlates every record in DATA1 against
%     every record in DATA2 (meaning that records in DATA1 are the 'master'
%     records).
%
%     CORRELATE(...,'normxc',false) leaves the correlograms unnormalized
%     (they are normalized by their autocorrelations by default).
%
%    Header Changes:  B, E, DEPMEN, DEPMAX, DEPMIN, KT0-3
%     The master record for each correlation is used for remaining header
%     values.  Slave record's KNETWK, KSTNM, KHOLE, KCMPNM are pushed into
%     fields KT0-3.  B, E are the lag range in seconds adjusted for the
%     differential beginning times of the records. All records are checked
%     that they share the same reference time and delta.
%
%    Usage: correlograms=correlate(data)
%           correlograms=correlate(data,'normxc',false)
%           correlograms=correlate(data1,data2)
%           correlograms=correlate(data1,data2,'normxc',false)
%
%    Examples:
%     Equivalent to 'correlate' in SAC:
%     correlograms=correlate(data(1),data,'normxc',false)
%
%    See also: convolve, dft, idft

% split based on npeaks (0 or 1+)           check
% check datsets                             check
% peaks get added to fields .cg .lg .pg     check
% correlograms replace dataset
% figure out lags
%   lag_slave + b_slave - b_master

% check data structure
error(seischk(data1,'dep'))

% 1 or 2 datasets
nrecs1=numel(data1);
if(nargin>1 && isseis(varargin{1},'dep')); 
    ndata=2; nrecs2=numel(varargin{1});
else
    ndata=1;
end

% how many peaks
npeaks=0;
for i=1:length(varargin)
    if(ischar(varargin{i}) && strcmpi('npeaks',varargin{i}))
        npeaks=varargin{i+1};
    end
end

% check npeaks
if(~isscalar(npeaks) && ~isnumeric(npeaks) && fix(npeaks)~=npeaks)
    error('SAClab:correlate:badInput',...
        'Option NPEAKS must be a scalar integer')
end

% split based on npeaks
% npeaks>0 ==> assign grids to new fields
if(npeaks)
    % split based on ndatasets
    if(ndata==2)
        % check datasets
        delta=check_data(data1,varargin{1});
        % get relative start times (disregarding absolute timing)
        b1=gh(data1,'b');
        b2=gh(varargin{1},'b');
        bdiff=b2(:,ones(nrecs1,1)).'-b1(:,ones(nrecs2,1));
        % extract records
        data1=combo(data1);
        varargin{1}=combo(varargin{1});
        % get correlation peaks
        [data.cg,data.lg,data.pg]=mcxc(data1,varargin{:});
        % adjust lags
        s=size(data.lg); s(ndims(s):4)=1;
        data.lg=data.lg*delta+bdiff(:,:,ones(s(3),1),ones(s(4),1));
    else
        % check datasets
        delta=check_data(data1);
        % get relative start times (disregarding absolute timing)
        b1=gh(data1,'b');
        bdiff=(b1(:,ones(nrecs1,1)).'-b1(:,ones(nrecs1,1))).';
        bdiff=bdiff(tril(true(nrecs1),-1)); % extract upper triangle
        % extract records
        data1=combo(data1);
        % get correlation peaks
        [data.cg,data.lg,data.pg]=mcxc(data1,varargin{:});
        % adjust lags
        s=size(data.lg); s(ndims(s):4)=1;
        data.lg=data.lg*delta+bdiff(:,:,ones(s(3),1),ones(s(4),1));
    end
% npeaks==0 ==> replace dataset with correlograms
else
    data=old_corr(data1,varargin{:});
    
    % the new stuff (won't run unless the next line is commented out)
    if(1); return; end
    % split based on ndatasets
    if(ndata==2)
        % what do we need to do
        %  - setup data headers for correlograms
        %  - permute the correlograms crazy like
        %  - pump in correlogram data
        %  - cut it up right
        %  - final header updates
    else
        % what do we need to do
        %  - setup data headers for correlograms
        %  - permute the correlograms crazy like
        %  - pump in correlogram data
        %  - cut it up right
        %  - final header updates
        
        % check datasets
        delta=check_data(data1);
        % get relative start times (disregarding absolute timing)
        b1=gh(data1,'b');
        bdiff=(b1(:,ones(nrecs1,1)).'-b1(:,ones(nrecs1,1))).';
        bdiff=bdiff(tril(true(nrecs1),-1)); % extract upper triangle
        % extract records
        data1=combo(data1);
        % get correlation peaks
        [data.cg,data.lg,data.pg]=mcxc(data1,varargin{:});
        % adjust lags
        s=size(data.lg); s(ndims(s):4)=1;
        data.lg=data.lg*delta+bdiff(:,:,ones(s(3),1),ones(s(4),1));
    end
end

end

function [delta]=check_data(data1,data2)
% assure datasets are capable of being sensibly cross correlated

% if 2 datasets, datasets must be consistent with one another
if(nargin==2)
    % check delta,iftype,leven
    [delta]=gh([data1; data2],'delta');
    iftype=genumdesc([data1; data2],'iftype');
    if(~isscalar(unique(delta)))
        error('SAClab:correlate:mixedSampleRates',...
            'Mixed sample rates not allowed');
    elseif(any(~strcmp(iftype,'Time Series File')...
            & ~strcmp(iftype,'General X vs Y file')))
        error('SAClab:correlate:illegalOperation',...
            'Illegal operation on non-Time Series File')
    end
    if(any(~strcmp(glgc([data1; data2],'leven'),'true')))
        error('SAClab:correlate:evenlySpacedOnly',...
            'Illegal operation on unevenly spaced data');
    end
    
    % check ncmp
    for i=1:length(data1)
        if(size(data1(i).x,2)~=1)
            error('SAClab:correlate:tooManyComponents',...
                'Multi-component data is not supported!')
        end
    end
    for i=1:length(data2)
        if(size(data1(2).x,2)~=1)
            error('SAClab:correlate:tooManyComponents',...
                'Multi-component data is not supported!')
        end
    end
% if 1 dataset, needs to be consistent with itself
else
    % check delta,iftype,leven
    [delta]=gh(data1,'delta');
    iftype=genumdesc(data1,'iftype');
    if(~isscalar(unique(delta)))
        error('SAClab:correlate:mixedSampleRates',...
            'Mixed sample rates not allowed');
    elseif(any(~strcmp(iftype,'Time Series File')...
            & ~strcmp(iftype,'General X vs Y file')))
        error('SAClab:correlate:illegalOperation',...
            'Illegal operation on non-Time Series File')
    end
    if(any(~strcmp(glgc(data1,'leven'),'true')))
        error('SAClab:correlate:evenlySpacedOnly',...
            'Illegal operation on unevenly spaced data');
    end
    
    % check ncmp
    for i=1:length(data1)
        if(size(data1(i).x,2)~=1)
            error('SAClab:correlate:tooManyComponents',...
                'Multi-component data is not supported!')
        end
    end
end

% return common sample rate
delta=delta(1);

end

function [data]=old_corr(data1,varargin)

%%%%%%%%%%%%%%%%%%%
%%% OLD VERSION %%%
%%%%%%%%%%%%%%%%%%%

% check data structure
error(seischk(data1,'x'))

% cross correlate 2 data sets or 1
if(nargin>1 && isseis(varargin{1},'x'))
    % two data sets
    data2=varargin{1};
    varargin(1)=[];
    
    % checks (delta,iftype,leven)
    [delta,nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec]=...
        gh([data1; data2],'delta','nzyear','nzjday',...
        'nzhour','nzmin','nzsec','nzmsec');
    iftype=genumdesc([data1; data2],'iftype');
    if(~isscalar(unique(delta)))
        error('SAClab:correlate:mixedSampleRates',...
            'Mixed sample rates not allowed');
    elseif(~isscalar(unique(nzyear)) || ~isscalar(unique(nzjday)) ...
            || ~isscalar(unique(nzhour)) || ~isscalar(unique(nzmin)) ...
            || ~isscalar(unique(nzsec)) || ~isscalar(unique(nzmsec)))
        error('SAClab:correlate:mixedReferenceTimes',...
            'Mixed reference times not allowed');
    elseif(any(~strcmp(iftype,'Time Series File')...
            & ~strcmp(iftype,'General X vs Y file')))
        error('SAClab:correlate:illegalOperation',...
            'Illegal operation on non-Time Series File')
    end
    if(any(~strcmp(glgc([data1; data2],'leven'),'true')))
        error('SAClab:correlate:evenlySpacedOnly',...
            'Illegal operation on unevenly spaced data');
    end
    
    % size'em up
    s1=length(data1);
    s2=length(data2);
    k=s2*s1;
    
    % preallocate data
    data(k,1)=data2(s2);
    
    % some header values
    b1=gh(data1,'b');
    [knetwk,kstnm,khole,kcmpnm,b2]=...
        gh(data2,'knetwk','kstnm','khole','kcmpnm','b');
    
    % loop through every pair (starts with last pair)
    for i=s2:-1:1
        for j=s1:-1:1
            % copy header
            data(k)=data2(i);
            
            % cross correlate
            [data(k).x,lg]=mcxc(data1(j).x,data2(i).x,varargin{:});
            data(k).x=permute(data(k).x,[3 2 1]);
            lg=lg*delta(1)+b2(i)-b1(j);
            
            % update header
            data(k)=ch(data(k),'b',lg(1),...
                'e',lg(end),...
                'npts',length(lg),...
                'depmen',mean(data(k).x),...
                'depmin',min(data(k).x),...
                'depmax',max(data(k).x),...
                'kt0',knetwk{j},'kt1',kstnm{j},...
                'kt2',khole{j},'kt3',kcmpnm{j});
            
            % increment
            k=k-1;
        end
    end
else
    % one data set
    
    % checks (delta,iftype,leven,ref time)
    [delta,nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec]=...
        gh(data1,'delta','nzyear','nzjday','nzhour',...
        'nzmin','nzsec','nzmsec');
    iftype=genumdesc(data1,'iftype');
    if(~isscalar(unique(delta)))
        error('SAClab:correlate:mixedSampleRates',...
            'Mixed sample rates not allowed');
    elseif(~isscalar(unique(nzyear)) || ~isscalar(unique(nzjday)) ...
            || ~isscalar(unique(nzhour)) || ~isscalar(unique(nzmin)) ...
            || ~isscalar(unique(nzsec)) || ~isscalar(unique(nzmsec)))
        error('SAClab:correlate:mixedReferenceTimes',...
            'Mixed reference times not allowed');
    elseif(any(~strcmp(iftype,'Time Series File')...
            & ~strcmp(iftype,'General X vs Y file')))
        error('SAClab:correlate:illegalOperation',...
            'Illegal operation on non-Time Series File')
    end
    if(any(~strcmp(glgc(data1,'leven'),'true')))
        error('SAClab:correlate:evenlySpacedOnly',...
            'Illegal operation on unevenly spaced data');
    end
    
    % size'em up
    s=length(data1);
    k=(s*(s-1))/2;
    
    % preallocate data
    data(k,1)=data1(s-1);
    
    % some header values
    [knetwk,kstnm,khole,kcmpnm,b]=...
        gh(data1,'knetwk','kstnm','khole','kcmpnm','b');
    
    % loop through every pair (starts with last pair)
    for i=s-1:-1:1
        for j=s:-1:i+1
            % copy header
            data(k)=data1(i);
            
            % cross correlate
            [data(k).x,lg]=mcxc(data1(j).x,data1(i).x,varargin{:});
            data(k).x=permute(data(k).x,[3 2 1]);
            
            % get relative timing
            lg=lg*delta(1)+b(i)-b(j);
            
            % update header
            data(k)=ch(data(k),'b',lg(1),...
                'e',lg(end),...
                'npts',length(lg),...
                'depmen',mean(data(k).x),...
                'depmin',min(data(k).x),...
                'depmax',max(data(k).x),...
                'kt0',knetwk{j},'kt1',kstnm{j},...
                'kt2',khole{j},'kt3',kcmpnm{j});
            
            % increment
            k=k-1;
        end
    end
end

end
