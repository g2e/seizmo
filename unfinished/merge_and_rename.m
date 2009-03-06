function [data]=merge_and_rename(data,varargin)
%MERGE    Merge SEIZMO records
%
%
%

% todo:
%   - uneven support - just toss together and sort after?
%   - make description
%   - add rest of gap/lap methods

% check nargin
if(mod(nargin-1,2))
    error('seizmo:merge:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data,'dep'))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% defaults
option.TOLERANCE=0.2;
option.OVERLAP='sequential';
option.GAP='sequential';
option.FILLMETHOD='spline';
option.SEQUENTIALMETHOD='longer';
option.TRUNCATEMETHOD='longer';
option.FILLER=0;
option.TIMING='tai';
option.USEABSOLUTETIMING=true;
option.REQUIREDCHARFIELDS={'knetwk' 'kstnm' 'khole' 'kcmpnm'};
option.REQUIREDREALFIELDS={'delta' 'cmpinc' 'cmpaz'};

% get options from SEIZMO global
global SEIZMO
fields=fieldnames(option);
for i=fields
    try
        option.(i)=SEIZMO.MERGE.(i);
    catch
        % do nothing
    end
end

% parse options
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:merge:badInput',...
            'Options must be specified as a strings!');
    end
    if((isnumeric(varargin{i+1}) && ~isscalar(varargin{i+1}))...
        || (ischar(varargin{i+1}) && size(varargin{i+1},1)~=1))
        error('seizmo:merge:badInput',...
            'Bad value for option %s !',varargin{i});
    end
    switch lower(varargin{i})
        case 'tolerance'
            if(~isnumeric(varargin{i+1}))
                error('seizmo:merge:badInput',...
                    'TOLERANCE must be a scalar number!');
            end
            option.TOLERANCE=varargin{i+1};
        case 'overlap'
            if(~ischar(varargin{i+1}) || ~any(strcmpi(varargin{i+1},...
                    {'sequential' 'truncate'})))
                error('seizmo:merge:badInput',...
                    ['OVERLAP option must be '...
                    '''sequential'' or ''truncate''!']);
            end
            option.OVERLAP=lower(varargin{i+1});
        case 'gap'
            if(~ischar(varargin{i+1}) || ~any(strcmpi(varargin{i+1},...
                    {'sequential' 'fill'})))
                error('seizmo:merge:badInput',...
                    ['GAP option must be '...
                    '''sequential'' or ''fill''!']);
            end
            option.GAP=lower(varargin{i+1});
        case 'fillmethod'
            if(~ischar(varargin{i+1}) || ~any(strcmpi(varargin{i+1},...
                    {'value' 'nearest' 'linear' 'pchip' 'spline'})))
                error('seizmo:merge:badInput',...
                    ['FILLMETHOD option must be ''value'' '...
                    '''nearest'' ''linear'' ''pchip'' or ''spline''!']);
            end
            option.FILLMETHOD=lower(varargin{i+1});
        case 'sequentialmethod'
            if(~ischar(varargin{i+1}) || ~any(strcmpi(varargin{i+1},...
                    {'first' 'last' 'shorter' 'longer'})))
                error('seizmo:merge:badInput',...
                    ['SEQUENTIALMETHOD option must be '...
                    '''first'' ''last'' ''shorter'' or ''longer''!']);
            end
            option.SEQUENTIALMETHOD=lower(varargin{i+1});
        case 'truncatemethod'
            if(~ischar(varargin{i+1}) || ~any(strcmpi(varargin{i+1},...
                    {'first' 'last' 'shorter' 'longer'})))
                error('seizmo:merge:badInput',...
                    ['TRUNCATEMETHOD option must be '...
                    '''first'' ''last'' ''shorter'' or ''longer''!']);
            end
            option.TRUNCATEMETHOD=lower(varargin{i+1});
        case 'filler'
            if(~isnumeric(varargin{i+1}))
                error('seizmo:merge:badInput',...
                    'FILLER must be a scalar number!');
            end
            option.FILLER=varargin{i+1};
        case 'timing'
            if(~ischar(varargin{i+1}) || ~any(strcmpi(varargin{i+1},...
                    {'tai' 'utc'})))
                error('seizmo:merge:badInput',...
                    ['TIMING option must be '...
                    '''tai'' or ''utc''!']);
            end
            option.TIMING=lower(varargin{i+1});
        case 'useabsolutetiming'
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:merge:badInput',...
                    'USEABSOLUTETIMING option must be a logical!');
            end
            option.USEABSOLUTETIMING=varargin{i+1};
        case 'requiredcharfields'
            if(~iscellstr(varargin{i+1}))
                error('seizmo:merge:badInput',...
                    ['REQUIREDCHARFIELDS option must be '...
                    'a cellstr array of header fields!']);
            end
            option.REQUIREDCHARFIELDS=varargin{i+1};
        case 'requiredrealfields'
            if(~iscellstr(varargin{i+1}))
                error('seizmo:merge:badInput',...
                    ['REQUIREDREALFIELDS option must be '...
                    'a cellstr array of header fields!']);
            end
            option.REQUIREDREALFIELDS=varargin{i+1};
        otherwise
            error('seizmo:merge:badInput',...
                'Unknown option: %s !',varargin{i});
    end
end

% get header fields
ncmp=getncmp(data);
if(option.USEABSOLUTETIMING)
    [b,e,delta,npts,depmin,depmax,depmen,...
        nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec]=getheader(data,...
        'b','e','delta','npts','depmin','depmax','depmen',...
        'nzyear','nzjday','nzhour','nzmin','nzsec','nzmsec');
else
    [b,e,delta,npts,depmin,depmax,depmen]=...
        getheader(data,'b','e','delta','npts','depmin','depmax','depmen');
end
szreal=size(option.REQUIREDREALFIELDS); reqreal=cell(szreal);
szchar=size(option.REQUIREDCHARFIELDS); reqchar=cell(szchar);
if(prod(szreal)~=0)
    [reqreal{:}]=getheader(data,option.REQUIREDREALFIELDS{:});
end
if(prod(szchar)~=0)
    [reqchar{:}]=getheader(data,option.REQUIREDCHARFIELDS{:});
end
iftype=getenumid(data,'iftype');
leven=strcmp(getlgc(data,'leven'),'true');

% require timeseries and general x vs y
if(any(~strcmp(iftype,'itime') & ~strcmp(iftype,'ixy')))
    error('seizmo:merge:badRecordType',...
        'Records must be Time Series or General X vs Y !');
end

% get start and end of records in absolute time
nrecs=numel(data);
if(option.USEABSOLUTETIMING)
    if(strcmp(option.TIMING,'utc'))
        ab=gregorian2modserial(utc2tai(...
            [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+b]));
        ae=gregorian2modserial(utc2tai(...
            [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+e]));
    else
        ab=gregorian2modserial(...
            [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+b]);
        ae=gregorian2modserial(...
            [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000+e]);
    end
else
    ab=[zeros(nrecs,1) b];
    ae=[zeros(nrecs,1) e];
end

% change real to char
for i=1:prod(szreal)
    reqreal{i}=num2str(reqreal{i},'%16.16e');
end

% make groups (require at least leven and ncmp to be the same)
[f,h,h]=unique(char(strcat(strcat('',reqchar{:}),'_',...
    strcat('',reqreal{:}),'_',num2str(leven),'_',num2str(ncmp))),...
    'rows');

% get some header info for renaming
[knetwk,kstnm,khole,kcmpnm]=...
    gh(data,'knetwk','kstnm','khole','kcmpnm');

% fix empty khole
empty=strcmp(khole,'');
khole(empty)={'__'};
data(empty)=ch(data(empty),'khole','__');

% loop through each group
destroy=false(nrecs,1);
for i=1:size(f,1)
    % get group member indices
    gidx=find(h==i);
    ng=numel(gidx);
    
    % rename
    names=strcat(knetwk(gidx(1)),'.',kstnm(gidx(1)),'.',khole(gidx(1)),...
        '.',kcmpnm(gidx(1)),'.',num2str((1:ng)','%02d'));
    [data(gidx).name]=deal(names{:});
    
    % no records to merge with
    if(ng==1); continue; end
    
    % find duplicates
    dups=(ab(gidx,ones(ng,1))>ab(gidx,ones(ng,1)).' ...
        | (ab(gidx,ones(ng,1))==ab(gidx,ones(ng,1)).' ...
        & ab(gidx,2*ones(ng,1))>=ab(gidx,2*ones(ng,1)).')) ...
        & (ae(gidx,ones(ng,1))<ae(gidx,ones(ng,1)).' ...
        | (ae(gidx,ones(ng,1))==ae(gidx,ones(ng,1)).' ...
        & ae(gidx,2*ones(ng,1))<=ae(gidx,2*ones(ng,1)).'));
    dupsu=(dups-dups.')>0;    % only delete if other not deleted
    dupsu(tril(true(ng)))=0;  % upper triangle
    dups(triu(true(ng)))=0;   % lower triangle
    dups=sum(dups+dupsu,2)>0; % logical indices in group
    
    % remove duplicates
    destroy(gidx(dups))=true;
    gidx(dups)=[];
    ng=numel(gidx);
    
    % merge
    % looped so multiple merges are sane
    % if a vectorized solution exists please let me know :D
    for j=1:(ng-1)
        if(destroy(gidx(j))); continue; end
        for k=2:ng
            if(destroy(gidx(j))); continue; end
            if(destroy(gidx(k))); continue; end
            
            % handle uneven
            if(~leven(j))
                % what to do?
                % - pass to mseq as if perfect
                % - merge .ind
                % - sort by .ind
                % - if diff .ind==0 whine
                continue;
            end
            
            % only find mergible (sequential within tolerance)
            diff12=delta(gidx(j))+(ae(gidx(j),1)-ab(gidx(k),1))*86400 ...
                +(ae(gidx(j),2)-ab(gidx(k),2));
            diff21=delta(gidx(j))+(ae(gidx(k),1)-ab(gidx(j),1))*86400 ...
                +(ae(gidx(k),2)-ab(gidx(j),2));
            if(diff12==0)
                % perfectly mergible, 1st file first
                [data(gidx([j k])),destroy(gidx([j k])),...
                    ab(gidx([j k]),:),ae(gidx([j k]),:),...
                    npts(gidx([j k])),depmin(gidx([j k])),...
                    depmax(gidx([j k])),depmen(gidx([j k]))]=...
                    mseq(data(gidx([j k])),ab(gidx([j k]),:),...
                    ae(gidx([j k]),:),npts(gidx([j k])),...
                    depmin(gidx([j k])),depmax(gidx([j k])),...
                    depmen(gidx([j k])),diff12,option,1);
            elseif(diff21==0)
                % perfectly mergible, 2nd file first
                [data(gidx([j k])),destroy(gidx([j k])),...
                    ab(gidx([j k]),:),ae(gidx([j k]),:),...
                    npts(gidx([j k])),depmin(gidx([j k])),...
                    depmax(gidx([j k])),depmen(gidx([j k]))]=...
                    mseq(data(gidx([j k])),ab(gidx([j k]),:),...
                    ae(gidx([j k]),:),npts(gidx([j k])),...
                    depmin(gidx([j k])),depmax(gidx([j k])),...
                    depmen(gidx([j k])),diff21,option,2);
            elseif(diff12>0 && diff12<option.TOLERANCE)
                % mergible gap, 1st file first
                [data(gidx([j k])),destroy(gidx([j k])),...
                    ab(gidx([j k]),:),ae(gidx([j k]),:),...
                    npts(gidx([j k])),depmin(gidx([j k])),...
                    depmax(gidx([j k])),depmen(gidx([j k]))]=...
                    mgap(data(gidx([j k])),ab(gidx([j k]),:),...
                    ae(gidx([j k]),:),npts(gidx([j k])),...
                    depmin(gidx([j k])),depmax(gidx([j k])),...
                    depmen(gidx([j k])),diff12,option,1);
            elseif(diff21>0 && diff21<option.TOLERANCE)
                % mergible gap, 2nd file first
                [data(gidx([j k])),destroy(gidx([j k])),...
                    ab(gidx([j k]),:),ae(gidx([j k]),:),...
                    npts(gidx([j k])),depmin(gidx([j k])),...
                    depmax(gidx([j k])),depmen(gidx([j k]))]=...
                    mgap(data(gidx([j k])),ab(gidx([j k]),:),...
                    ae(gidx([j k]),:),npts(gidx([j k])),...
                    depmin(gidx([j k])),depmax(gidx([j k])),...
                    depmen(gidx([j k])),diff21,option,2);
            elseif(diff12<0 && abs(diff12)<option.TOLERANCE)
                % mergible overlap, 1st file first
                [data(gidx([j k])),destroy(gidx([j k])),...
                    ab(gidx([j k]),:),ae(gidx([j k]),:),...
                    npts(gidx([j k])),depmin(gidx([j k])),...
                    depmax(gidx([j k])),depmen(gidx([j k]))]=...
                    mlap(data(gidx([j k])),ab(gidx([j k]),:),...
                    ae(gidx([j k]),:),npts(gidx([j k])),...
                    depmin(gidx([j k])),depmax(gidx([j k])),...
                    depmen(gidx([j k])),diff12,option,1);
            elseif(diff21<0 && abs(diff21)<option.TOLERANCE)
                % mergible overlap, 2nd file first
                [data(gidx([j k])),destroy(gidx([j k])),...
                    ab(gidx([j k]),:),ae(gidx([j k]),:),...
                    npts(gidx([j k])),depmin(gidx([j k])),...
                    depmax(gidx([j k])),depmen(gidx([j k]))]=...
                    mlap(data(gidx([j k])),ab(gidx([j k]),:),...
                    ae(gidx([j k]),:),npts(gidx([j k])),...
                    depmin(gidx([j k])),depmax(gidx([j k])),...
                    depmen(gidx([j k])),diff21,option,2);
            elseif(diff12<0 && abs(diff12)<option.TOLERANCE)
            end
        end
    end
end

% get relative times from absolute
if(option.USEABSOLUTETIMING)
    if(strcmp(option.TIMING,'utc'))
        az=gregorian2modserial(utc2tai(...
            [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000]));
        b=(ab(:,1)-az(:,1))*86400+(ab(:,2)-az(:,2));
        e=(ae(:,1)-az(:,1))*86400+(ae(:,2)-az(:,2));
    else
        az=gregorian2modserial(...
            [nzyear nzjday nzhour nzmin nzsec+nzmsec/1000]);
        b=(ab(:,1)-az(:,1))*86400+(ab(:,2)-az(:,2));
        e=(ae(:,1)-az(:,1))*86400+(ae(:,2)-az(:,2));
    end
else
    b=ab(:,2);
    e=ae(:,2);
end

% update header
data=changeheader(data,'b',b,'e',e,'delta',delta,'npts',npts,...
    'depmin',depmin,'depmax',depmax,'depmen',depmen);

% remove unwanted data
data(destroy)=[];

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end

function [data,destroy,ab,ae,npts,depmin,depmax,depmen]=...
    mseq(data,ab,ae,npts,depmin,depmax,depmen,diff,option,first)
%MSEQ    Merge sequential records

% who's header do we keep
destroy=false(2,1);
last=3-first;
long=find(npts==max(npts),1,'first');
short=3-long;
switch option.SEQUENTIALMETHOD
    case 'first'
        keep=first;
        kill=last;
        ab(keep,:)=ab(first,:);
        ae(keep,:)=ae(last,:)-[0 diff];
    case 'last'
        keep=last;
        kill=first;
        ab(keep,:)=ab(first,:)+[0 diff];
        ae(keep,:)=ae(last,:);
    case 'longer'
        keep=long;
        kill=short;
        if(keep==first)
            ab(keep,:)=ab(first,:);
            ae(keep,:)=ae(last,:)-[0 diff];
        else
            ab(keep,:)=ab(first,:)+[0 diff];
            ae(keep,:)=ae(last,:);
        end
    otherwise % shorter
        keep=short;
        kill=long;
        if(keep==first)
            ab(keep,:)=ab(first,:);
            ae(keep,:)=ae(last,:)-[0 diff];
        else
            ab(keep,:)=ab(first,:)+[0 diff];
            ae(keep,:)=ae(last,:);
        end
end

% merge
data(keep).dep=[data(first).dep; data(last).dep];

% update header fields
depmin(keep)=min(data(keep).dep(:));
depmax(keep)=max(data(keep).dep(:));
depmen(keep)=mean(data(keep).dep(:));

% update npts, destroy
npts(keep)=sum(npts);
destroy(kill)=true;
    
end


function [data,destroy,ab,ae,npts,depmin,depmax,depmen]=...
    mgap(data,ab,ae,npts,depmin,depmax,depmen,diff,option,first)
%MGAP    Merge gaps

switch option.GAP
    case 'sequential'
        [data,destroy,ab,ae,npts,depmin,depmax,depmen]=...
            mseq(data,ab,ae,npts,depmin,depmax,depmen,...
            diff,option,first);
    otherwise % fill
        
end

end


function [data,destroy,ab,ae,npts,depmin,depmax,depmen]=...
    mlap(data,ab,ae,npts,depmin,depmax,depmen,diff,option,first)
%MLAP    Merge overlaps

switch option.OVERLAP
    case 'sequential'
        [data,destroy,ab,ae,npts,depmin,depmax,depmen]=...
            mseq(data,ab,ae,npts,depmin,depmax,depmen,...
            diff,option,first);
    otherwise % truncate
        
end

end
