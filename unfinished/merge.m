function [data]=merge(data,varargin)
%MERGE    Merge SEIZMO records
%
%
%

%     Version History:
%        Dec.  6, 2008 - initial version
%        Dec.  8, 2008 - more options
%        Mar.  5, 2009 - 
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  8, 2008 at 21:20 GMT

% todo:
%   - uneven support - just toss together and sort after?
%   - make description
%   - add rest of gap/lap methods
%   - handle delta

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
option.TOLERANCE=0.02; % seconds, any positive number
option.OVERLAP='sequential'; % sequential/truncate
option.GAP='sequential'; % sequential/interpolate/fill
option.INTERPOLATE='spline'; % spline/pchip/linear/nearest
option.ADJUST='shorter'; % longer/shorter/first/last
option.ADJUSTUNITS='intervals'; % seconds/intervals
option.ADJUSTMAX=0.5; % interval: 0-0.5 , seconds: 0+
option.FILLER=0; % any number
option.TIMING='utc'; % utc/tai
option.USEABSOLUTETIMING=true; % true/false
option.REQUIREDCHARFIELDS={'knetwk' 'kstnm' 'khole' 'kcmpnm'};
option.REQUIREDREALFIELDS={'delta' 'cmpinc' 'cmpaz'};

% get options from SEIZMO global
global SEIZMO
fields=fieldnames(SEIZMO.MERGE);
for i=fields
    option.(i)=SEIZMO.MERGE.(i);
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:merge:badInput',...
            'Options must be specified as a strings!');
    end
    option.(upper(varargin{i}))=varargin{i+1};
end

% check options
fields=fieldnames(option);
for i=fields
    % get value of field and do a basic check
    value=option.i;
    if((isnumeric(value) && ~isscalar(value))...
        || (ischar(value) && size(value,1)~=1))
        error('seizmo:merge:badInput',...
            'Bad value for option %s !',value);
    end
    
    % specific checks
    switch lower(i)
        case 'tolerance'
            if(~isnumeric(value))
                error('seizmo:merge:badInput',...
                    'TOLERANCE must be a scalar number!');
            end
        case 'overlap'
            if(~ischar(value) || ~any(strcmpi(value,...
                    {'sequential' 'truncate'})))
                error('seizmo:merge:badInput',...
                    ['OVERLAP option must be '...
                    '''sequential'' or ''truncate''!']);
            end
        case 'gap'
            if(~ischar(value) || ~any(strcmpi(value,...
                    {'sequential' 'fill'})))
                error('seizmo:merge:badInput',...
                    ['GAP option must be '...
                    '''sequential'' or ''fill''!']);
            end
        case 'interpolate'
            if(~ischar(value) || ~any(strcmpi(value,...
                    {'nearest' 'linear' 'pchip' 'spline'})))
                error('seizmo:merge:badInput',...
                    ['INTERPOLATE option must be '...
                    '''nearest'' ''linear'' ''pchip'' or ''spline''!']);
            end
        case 'adjust'
            if(~ischar(value) || ~any(strcmpi(value,...
                    {'first' 'last' 'shorter' 'longer'})))
                error('seizmo:merge:badInput',...
                    ['ADJUST option must be '...
                    '''first'' ''last'' ''shorter'' or ''longer''!']);
            end
        case 'adjustunits'
            if(~ischar(value) || ~any(strcmpi(value,...
                    {'seconds' 'intervals'})))
                error('seizmo:merge:badInput',...
                    ['ADJUSTUNITS option must be '...
                    '''seconds'' or ''intervals''!']);
            end
        case 'adjustmax'
            if(~isnumeric(value))
                error('seizmo:merge:badInput',...
                    'ADJUSTMAX must be a scalar number!');
            end
        case 'filler'
            if(~isnumeric(value))
                error('seizmo:merge:badInput',...
                    'FILLER must be a scalar number!');
            end
        case 'timing'
            if(~ischar(value) || ~any(strcmpi(value,...
                    {'tai' 'utc'})))
                error('seizmo:merge:badInput',...
                    'TIMING option must be ''tai'' or ''utc''!');
            end
        case 'useabsolutetiming'
            if(~islogical(value) || ~isscalar(value))
                error('seizmo:merge:badInput',...
                    'USEABSOLUTETIMING option must be a logical!');
            end
        case 'requiredcharfields'
            if(~iscellstr(value))
                error('seizmo:merge:badInput',...
                    ['REQUIREDCHARFIELDS option must be '...
                    'a cellstr array of header fields!']);
            end
        case 'requiredrealfields'
            if(~iscellstr(value))
                error('seizmo:merge:badInput',...
                    ['REQUIREDREALFIELDS option must be '...
                    'a cellstr array of header fields!']);
            end
        otherwise
            error('seizmo:merge:badInput',...
                'Unknown option: %s !',i);
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

% loop through each group
destroy=false(nrecs,1);
for i=1:size(f,1)
    % get group member indices
    gidx=find(h==i);
    ng=numel(gidx);
    
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
    
    % handle uneven
    if(~leven(j))
        % what to do?
        % - pass to mseq as if perfect
        % - merge .ind
        % - sort by .ind
        % - if diff .ind==0 whine
        continue;
    end
    
    % check for merge between every pair
    %   looped so multiple merges are sane
    %   if a vectorized solution exists please let me know :D
    for j=1:(ng-1)
        % make sure record is not deleted from previous merge
        if(destroy(gidx(j))); continue; end
        for k=2:ng
            % make sure these records are not deleted
            if(destroy(gidx(j))); break; end
            if(destroy(gidx(k))); continue; end
            
            % get time gap
            diff12=delta(gidx(j))+(ae(gidx(j),1)-ab(gidx(k),1))*86400 ...
                +(ae(gidx(j),2)-ab(gidx(k),2));
            diff21=delta(gidx(j))+(ae(gidx(k),1)-ab(gidx(j),1))*86400 ...
                +(ae(gidx(k),2)-ab(gidx(j),2));
            
            % only find mergible (sequential within tolerance)
            if(diff12==0)
                % perfectly mergible, 1st file first
                mergefunc=@mseq; tdiff=diff12; lead=1;
            elseif(diff21==0)
                % perfectly mergible, 2nd file first
                mergefunc=@mseq; tdiff=diff21; lead=2;
            elseif(diff12>0 && diff12<option.TOLERANCE)
                % mergible gap, 1st file first
                mergefunc=@mgap; tdiff=diff12; lead=1;
            elseif(diff21>0 && diff21<option.TOLERANCE)
                % mergible gap, 2nd file first
                mergefunc=@mgap; tdiff=diff21; lead=2;
            elseif(diff12<0 && abs(diff12)<option.TOLERANCE)
                % mergible overlap, 1st file first
                mergefunc=@mlap; tdiff=diff12; lead=1;
            elseif(diff21<0 && abs(diff21)<option.TOLERANCE)
                % mergible overlap, 2nd file first
                mergefunc=@mlap; tdiff=diff21; lead=2;
            else
                % data not mergible
                continue;
            end
            
            % merge the records
            %
            % why so many variables?
            % limits to one CHANGEHEADER call which speeds things up a lot
            [data(gidx([j k])),destroy(gidx([j k])),...
                ab(gidx([j k]),:),ae(gidx([j k]),:),delta(gidx([j k])),...
                npts(gidx([j k])),depmin(gidx([j k])),...
                depmax(gidx([j k])),depmen(gidx([j k]))]=...
                mergefunc(data(gidx([j k])),ab(gidx([j k]),:),...
                ae(gidx([j k]),:),delta(gidx([j k])),npts(gidx([j k])),...
                depmin(gidx([j k])),depmax(gidx([j k])),...
                depmen(gidx([j k])),tdiff,option,lead);
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


function [data,destroy,ab,ae,delta,npts,depmin,depmax,depmen]=...
    mseq(data,ab,ae,delta,npts,depmin,depmax,depmen,diff,option,first)
%MSEQ    Merge sequential records

% who's header do we keep
destroy=false(2,1);
last=3-first;
long=find(npts==max(npts),1,'first');
short=3-long;
switch option.ADJUST
    case 'first'
        keep=last;
        kill=first;
        ab(keep,:)=ab(last,:)-[0 npts(first)*delta(last)];
    case 'last'
        keep=first;
        kill=last;
        ae(keep,:)=ae(first,:)+[0 npts(last)*delta(first)];
    case 'longer'
        keep=short;
        kill=long;
        if(keep==first)
            ae(keep,:)=ae(first,:)+[0 npts(last)*delta(first)];
        else
            ab(keep,:)=ab(last,:)-[0 npts(first)*delta(last)];
        end
    otherwise % shorter
        keep=long;
        kill=short;
        if(keep==first)
            ae(keep,:)=ae(first,:)+[0 npts(last)*delta(first)];
        else
            ab(keep,:)=ab(last,:)-[0 npts(first)*delta(last)];
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


function [data,destroy,ab,ae,delta,npts,depmin,depmax,depmen]=...
    mgap(data,ab,ae,delta,npts,depmin,depmax,depmen,diff,option,first)
%MGAP    Merge gaps

switch option.GAP
    case 'sequential'
        [data,destroy,ab,ae,npts,depmin,depmax,depmen]=...
            mseq(data,ab,ae,npts,depmin,depmax,depmen,...
            diff,option,first);
    case 'interpolate'
        % what to do?
        % - adjust the record or interpolate new samples
        % - interpolate the gap
    case 'fill'
        % what to do?
        % - adjust the record or interpolate new samples
        % - fill the gap with the filler
end

end


function [data,destroy,ab,ae,delta,npts,depmin,depmax,depmen]=...
    mlap(data,ab,ae,delta,npts,depmin,depmax,depmen,diff,option,first)
%MLAP    Merge overlaps

switch option.OVERLAP
    case 'sequential'
        % just shift the option.ADJUST record to make sequential
        [data,destroy,ab,ae,npts,depmin,depmax,depmen]=...
            mseq(data,ab,ae,npts,depmin,depmax,depmen,...
            diff,option,first);
    case 'truncate'
        % what to do?
        % - decide how many samples to delete
        % - adjust the remaining or interpolate
        
        % how many samples to delete?
        last=3-first;
        long=find(npts==max(npts),1,'first');
        short=3-long;
        switch option.ADJUST
            case 'first'
                keep=last;
                kill=first;
                
                % can we adjust?
                
                
                ab(keep,:)=ab(last,:)-[0 npts(first)*delta(last)];
            case 'last'
                keep=first;
                kill=last;
                ae(keep,:)=ae(first,:)+[0 npts(last)*delta(first)];
            case 'longer'
                keep=short;
                kill=long;
                if(keep==first)
                    ae(keep,:)=ae(first,:)+[0 npts(last)*delta(first)];
                else
                    ab(keep,:)=ab(last,:)-[0 npts(first)*delta(last)];
                end
            otherwise % shorter
                keep=long;
                kill=short;
                if(keep==first)
                    ae(keep,:)=ae(first,:)+[0 npts(last)*delta(first)];
                else
                    ab(keep,:)=ab(last,:)-[0 npts(first)*delta(last)];
                end
        end
        
        
        round(abs(diff())/delta())
end

end


function []=truncate()
%

end
