function [data]=addarrivals(data,varargin)
%ADDARRIVALS    Adds the indicated phases to the SEIZMO data header
%
%    Usage:    data=addarrivals(data)
%              data=addarrivals(...,'model',modelname,...)
%              data=addarrivals(...,'phases',phaselist,...)
%              data=addarrivals(...,'fields',indices,...)
%
%    Description: ADDARRIVALS(DATA) calls TAUPTIME to insert phase arrival
%     times into the headers of records in SEIZMO struct DATA.  The model
%     used to find phase arrival times is IASP91.  The phase list is
%     'ttall', which grabs nearly all phases of interest (more than can be
%     put into the SAC data header).  The first 10 arrivals from this list
%     are inserted into the t, kt, and user fields (t gets arrival times,
%     kt gets the phase names, and user gets the ray parameters).  Event
%     depth and distance are determined from the header fields evdp and
%     gcarc or dist or stla+stlo+evla+evlo.
%
%     ADDARRIVALS(...,'MODEL',MODELNAME,...) sets the model used in the
%     arrival time calculation to MODELNAME.  MODELNAME must be a cellstr
%     or char array of one model per record or a single model for all
%     records.  Accepts lots of different 1D models (see TauP for how to
%     add more models).  By default the MODELNAME is 'iasp91'.
%
%     ADDARRIVALS(...,'PHASES',PHASELIST,...) sets the phases to be added
%     to the header.  PHASELIST must be a cellstr or char array with one
%     comma-separated phase list per record or a single comma-separated
%     list for all records.  Accepts lots of different phases (see TauP for
%     a list and conventions).  By default PHASELIST is 'ttall', which
%     returns a large list similar to Brian Kennett's ttimes program would
%     when told to list all.
%
%     ADDARRIVALS(...,'FIELDS',INDICES,...) indicates the header field
%     indices of the t,kt,user fields to put arrivals in.  The default is
%     0:9.  The indices are for all records and can not be individually
%     set.
%
%    Notes:
%     - TauP was written by:
%       H. Philip Crotwell, Thomas J. Owens, Jeroen Ritsema
%     - matTaup was written by:
%       Qin Li
%
%    Header changes: T, KT, USER
%
%    Examples:
%     Fill t, kt, user fields 6:9 with several S phases:
%      data=addarrivals(data,'phases','tts+','fields',6:9);
%
%    See also: getarrival, tauptime, changeheader

%     Version History:
%        June 29, 2009 - initial version
%        Aug. 25, 2009 - description update (forgot fields option)
%        Sep.  2, 2009 - now uses tauptime
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  2, 2009 at 11:30 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:addarrivals:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% toggle off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% grab header setup
[h,vi]=versioninfo(data);

% default options
option.MODEL='iasp91';
option.PHASES='ttall';
option.FIELDS=0:9;

% get options from SEIZMO global
me=mfilename;
global SEIZMO
try
    fields=fieldnames(SEIZMO.(me));
    for i=1:numel(fields)
        if(~isempty(SEIZMO.(me).(fields{i})))
            option.(fields{i})=SEIZMO.(me).(fields{i});
        end
    end
catch
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:addarrivals:badInput',...
            'Options must be specified as a strings!');
    end
    if(~isempty(varargin{i+1}))
        option.(upper(varargin{i}))=varargin{i+1};
    end
end

% check options
nrecs=numel(data);
fields=fieldnames(option);
for i=1:numel(fields)
    % specific checks
    switch lower(fields{i})
        case 'model'
            if(iscellstr(option.(fields{i})))
                option.(fields{i})=char(option.(fields{i}));
            end
            if(~ischar(option.(fields{i})) ...
                    || ~any(size(option.(fields{i}),1)==[1 nrecs]))
                error('seizmo:addarrivals:badInput',...
                    ['MODEL must be a cellstr/char array with one\n'...
                    'model per record or a single model for all!']);
            end
            if(size(option.(fields{i}),1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
            option.(fields{i})=cellstr(option.(fields{i}));
        case 'phases'
            if(iscellstr(option.(fields{i})))
                option.(fields{i})=char(option.(fields{i}));
            end
            if(isempty(option.(fields{i})) || ...
                    ~ischar(option.(fields{i})) ...
                    || ~any(size(option.(fields{i}),1)==[1 nrecs]))
                error('seizmo:addarrivals:badInput',...
                    ['PHASES must be a cellstr/char array with one\n'...
                    'comma-separated phase list per record or a\n'...
                    'single comma-separated phase list for all!']);
            end
            if(size(option.(fields{i}),1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
            option.(fields{i})=cellstr(option.(fields{i}));
        case 'fields'
            if(~isempty(option.(fields{i})) && (...
                    numel(option.(fields{i}))>10 || ...
                    any(fix(option.(fields{i}))~=option.(fields{i})) || ...
                    any(option.(fields{i})<0 | option.(fields{i})>9)))
                error('seizmo:addarrivals:badInput',...
                    ['FIELDS must be a index array with values\n'...
                    'from 0 to 9 indicating the kt,t,user fields\n'...
                    'to be overwrote.  List is common to all records!']);
            end
    end
end

% get relevant header info
[evdp,gcarc,dist,stla,stlo,evla,evlo,o,t,kt,user]=getheader(data,...
    'evdp','gcarc','dist','stla','stlo','evla','evlo','o','t','kt','user');

% loop over records adding info to header
idx=option.FIELDS;
for i=1:nrecs
    % check header info
    if(evdp(i)==h(vi(i)).undef.ntype)
        warning('seizmo:addarrivals:badEVDP',...
            'Record: %d\nEVDP undefined! Treating as zero.',i);
        evdp(i)=0;
    end
    if(o(i)==h(vi(i)).undef.ntype)
        warning('seizmo:addarrivals:badO',...
            'Record: %d\nO field undefined! Treating as zero.',i);
        o(i)=0;
    end
    if(gcarc(i)~=h(vi(i)).undef.ntype)
        location={'deg' gcarc(i)};
    elseif(dist(i)~=h(vi(i)).undef.ntype)
        location={'km' dist(i)};
    elseif(all([stla(i) stlo(i) evla(i) evlo(i)]~=h(vi(i)).undef.ntype))
        location={'sta' [stla(i) stlo(i)] 'evt' [evla(i) evlo(i)]};
    else
        error('seizmo:addarrivals:badLocation',...
            ['Record: %d\nGCARC, DIST, or '...
            'STLA+STLO+EVLA+EVLO must be set to get arrivals!'],i)
    end
    
    % get arrivals
    arrivals=tauptime('mod',option.MODEL{i},'h',evdp(i)/1000,...
        'ph',option.PHASES{i},location{:});
    
    % add arrivals
    for j=1:min(numel(idx),numel(arrivals))
        t(i,idx(j)+1)=arrivals(j).time+o(i);
        kt{i,idx(j)+1}=arrivals(j).phaseName;
        user(i,idx(j)+1)=arrivals(j).rayParam;
    end
end

% update header
data=changeheader(data,'t',t,'kt',kt,'user',user);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
