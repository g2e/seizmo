function [data,removed]=removeduplicates_and_rename(data,varargin)
%REMOVEDUPLICATES    Remove duplicate SEIZMO records
%
%    Description: REMOVEDUPLICATES(DATA) returns the SEIZMO record dataset
%     DATA without any duplicate records or partial pieces based on the
%     start (B) and end (E) fields.  This is only valid for Time Series and
%     X vs Y records.  Uses the reference time fields (NZ*) to allow
%     finding duplicates that do not share the same reference position.
%
%     REMOVEDUPLICATES(...,'TIMING','UTC') takes UTC leapseconds into
%     account.  The default option is 'TAI', which assumes records use a
%     timing standard with no leapseconds.
%
%     REMOVEDUPLICATES(...,'USEABSOLUTETIMING',FALSE) will not account for
%     differences in the reference position of records (it just blindly
%     assumes they all share the same reference position).
%
%     REMOVEDUPLICATES(...,'REQUIREDCHARFIELDS',FIELDS) sets the required
%     character header fields (such as 'kcmpnm') that must be equal between
%     two records before they will even be checked as duplicates. FIELDS is
%     a cellstr array of field names.  Default list is 'knetwk', 'kstnm',
%     'khole', and 'kcmpnm'.
%
%     REMOVEDUPLICATES(...,'REQUIREDREALFIELDS',FIELDS) sets the required
%     numerical header fields (such as 'delta') that must be equal between
%     two records before they will even be checked as duplicates. FIELDS is
%     a cellstr array of field names.  Default list is 'delta', 'cmpinc',
%     and 'cmpaz'.
%
%     [DATA,REMOVED]=REMOVEDUPLICATES(...) also returns a listing of the
%     indices of the records removed in REMOVED.  These indices are
%     relative to the input dataset, not the output dataset (obviously
%     because the records are no longer in the output dataset!).
%
%    Notes:
%     - MERGE removes duplicates and merges records
%     - duplicates are required to have identical LEVEN and NCMP fields
%
%    Tested on: Matlab r2007b
%
%    Header changes: NONE (does use CHECKHEADER though)
%
%    Usage:    data=removeduplicates(data)
%              data=removeduplicates(data,'timing','utc')
%              data=removeduplicates(data,'useabsolutetiming',true|false)
%              data=removeduplicates(data,'requiredcharfields',fields)
%              data=removeduplicates(data,'requiredrealfields',fields)
%              [data,removed]=removeduplicates(...)
%
%    Examples:
%     Datasets sometimes include duplicates and partial pieces.  Usually
%     this is associated with data pulled from a DHI server.
%      data=removeduplicates(data)
%
%    See also: merge, removedeadrecords

%     Version History:
%        Dec.  6, 2008 - initial version
%        Dec.  8, 2008 - more options
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  8, 2008 at 21:20 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:removeduplicates:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% defaults
option.TIMING='tai';
option.USEABSOLUTETIMING=true;
option.REQUIREDCHARFIELDS={'knetwk' 'kstnm' 'khole' 'kcmpnm'};
option.REQUIREDREALFIELDS={'delta' 'cmpinc' 'cmpaz'};

% get options from SEIZMO global
global SEIZMO
fields=fieldnames(option);
for i=fields
    try
        option.(i)=SEIZMO.REMOVEDUPLICATES.(i);
    catch
        % do nothing
    end
end

% parse options
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:removeduplicates:badInput',...
            'Options must be specified as a strings!');
    end
    if((isnumeric(varargin{i+1}) && ~isscalar(varargin{i+1}))...
        || (ischar(varargin{i+1}) && size(varargin{i+1},1)~=1))
        error('seizmo:removeduplicates:badInput',...
            'Bad value for option %s !',varargin{i});
    end
    switch lower(varargin{i})
        case 'timing'
            if(~ischar(varargin{i+1}) || ~any(strcmpi(varargin{i+1},...
                    {'tai' 'utc'})))
                error('seizmo:removeduplicates:badInput',...
                    ['TIMING option must be '...
                    '''tai'' or ''utc''!']);
            end
            option.TIMING=lower(varargin{i+1});
        case 'useabsolutetiming'
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:removeduplicates:badInput',...
                    'USEABSOLUTETIMING option must be a logical!');
            end
            option.USEABSOLUTETIMING=varargin{i+1};
        case 'requiredcharfields'
            if(~iscellstr(varargin{i+1}))
                error('seizmo:removeduplicates:badInput',...
                    ['REQUIREDCHARFIELDS option must be '...
                    'a cellstr array of header fields!']);
            end
            option.REQUIREDCHARFIELDS=varargin{i+1};
        case 'requiredrealfields'
            if(~iscellstr(varargin{i+1}))
                error('seizmo:removeduplicates:badInput',...
                    ['REQUIREDREALFIELDS option must be '...
                    'a cellstr array of header fields!']);
            end
            option.REQUIREDREALFIELDS=varargin{i+1};
        otherwise
            error('seizmo:removeduplicates:badInput',...
                'Unknown option: %s !',varargin{i});
    end
end

% get header fields
ncmp=getncmp(data);
if(option.USEABSOLUTETIMING)
    [b,e,nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec]=getheader(data,...
        'b','e','nzyear','nzjday','nzhour','nzmin','nzsec','nzmsec');
else
    [b,e]=getheader(data,'b','e');
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
    error('seizmo:removeduplicates:badRecordType',...
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

% loop through each group
destroy=false(numel(data),1);
for i=1:size(f,1)
    % get group member indices
    gidx=find(h==i);
    ng=numel(gidx);
    
    % rename
    names=strcat(knetwk(gidx(1)),'.',kstnm(gidx(1)),'.',khole(gidx(1)),...
        '.',kcmpnm(gidx(1)),'.',num2str((1:ng)','%02d'));
    [data(gidx).name]=deal(names{:});
    
    % no records to removeduplicates with
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
end

% remove unwanted data
data(destroy)=[];
removed=find(destroy);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end