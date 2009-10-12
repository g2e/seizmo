function [data]=changename(data,varargin)
%CHANGENAME    Change the filename of SEIZMO records
%
%    Usage:    data=changename(data,...,'name',name,...)
%              data=changename(data,...,'prepend',string,...)
%              data=changename(data,...,'append',string,...)
%              data=changename(data,...,'delete',string,...)
%              data=changename(data,...,'delete',{string1 ... stringN},...)
%              data=changename(data,...,'change',{orig replacement},...)
%
%    Description: CHANGENAME(DATA,...,'NAME',NAME,...) sets the name field
%     of DATA to NAME.  The name field is the filename associated with the
%     record.  If DATA contains more than 1 record and NAME is just a
%     single string, the name field for all records are changed to NAME.
%     Otherwise NAME must be a char/cellstr array with 1 row/element for
%     each record in DATA.
%
%     CHANGENAME(DATA,...,'PREPEND',STRING,...) prepends STRING to the
%     name field of DATA.  If DATA contains more than 1 record and STRING
%     is just a single string, STRING is prepended to the name field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.
%
%     CHANGENAME(DATA,...,'APPEND',STRING,...) appends STRING to the name
%     field of DATA.  If DATA contains more than 1 record and STRING is
%     just a single string, STRING is appended to the name field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.
%
%     CHANGENAME(DATA,...,'DELETE',STRING,...) deletes STRING from the name
%     field of DATA.  If DATA contains more than 1 record and STRING is
%     just a single string, STRING is deleted from the name field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.  To delete multiple parts from
%     the name field, see the next DELETE option usage.
%
%     CHANGENAME(DATA,...,'DELETE',{STRING1 ... STRINGN},...) allows
%     deleting multiple strings from the name field of DATA.  Multi-line
%     entries are not allowed (use one string per cell if working on
%     multiple records where each row corresponds to a separate record).  
%
%     CHANGENAME(DATA,...,'CHANGE',{ORIGINAL REPLACEMENT},...) replaces all
%     occurances of the string ORIGINAL in the name field of DATA with the
%     string REPLACEMENT.  If DATA contains more than 1 record and ORIGINAL
%     and REPLACEMENT are just single strings, then the same string set is
%     applied to every record.  Otherwise a string pair must be given for
%     each record (requires a cellstr array).  To change multiple portions
%     of the name field for records in DATA, just expand the number of
%     columns in the cellstr array.
%
%    Notes:
%     - CHANGENAME does NOT check the validity of input strings! 
%     - The order in which the options are implemented:
%        name, prepend, append, delete, change
%       This may come in handy for some more complicated cases.
%     - Using an option more than once replaces the previous option value
%       with the later one.
%
%    Examples:
%     Prepend and append to all filenames:
%      data=changename(data,'prepend','zzz','append','.new')
%
%     Change the name of a record:
%      data(3)=changename(data(3),'name','mynewfilename')
%
%     Delete multiple parts from all filenames:
%      data=changename(data,'delete',{'__' '.merged' 'SAC'})
%
%     Replace certain portions of all filenames:
%      data=changename(data,'change',{'..' '.__.' 'part1' 'part2'})
%
%    See also: CHANGEBYTEORDER, CHANGEPATH, WRITEHEADER, WRITESEIZMO

%     Version History:
%        May  28, 2009 - initial version
%        May  29, 2009 - fix empty option handling, allow scalar expansion
%                        for name
%        Sep.  5, 2009 - FILENAME and FILE also point to NAME option, error
%                        on unknown option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  5, 2009 at 05:55 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:changename:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% fast exit
if(nargin==1); return; end

% default options
option.NAME=[];
option.PREPEND=[];
option.APPEND=[];
option.DELETE=[];
option.CHANGE=[];

% get options from SEIZMO global
global SEIZMO
try
    fields=fieldnames(SEIZMO.CHANGENAME);
    for i=1:numel(fields)
        option.(fields{i})=SEIZMO.CHANGENAME.(fields{i});
    end
catch
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:changename:badInput',...
            'Options must be specified as a strings!');
    end
    if(strcmpi(varargin{i},'filename') || strcmpi(varargin{i},'file'))
        varargin{i}='name';
    end
    option.(upper(varargin{i}))=varargin{i+1};
end

% check options
nrecs=numel(data);
fields=fieldnames(option);
for i=1:numel(fields)
    % specific checks
    switch lower(fields{i})
        case {'name' 'append' 'prepend'}
            if(isempty(option.(fields{i}))); continue; end
            if(ischar(option.(fields{i})))
                option.(fields{i})=cellstr(option.(fields{i}));
            elseif(iscellstr(option.(fields{i})))
                % flatten multi-line names
                option.(fields{i})=cellstr(char(option.(fields{i})));
            end
            if(~iscellstr(option.(fields{i})) || ...
                    numel(size(option.(fields{i})))~=2 || ...
                    ~any(size(option.(fields{i}),1)==[1 nrecs]))
                error('seizmo:changename:badInput',...
                    ['%s option must be a cellstr/char array\n'...
                    'with a single string or one for each record!'],...
                    fields{i});
            end
        case 'delete'
            if(isempty(option.(fields{i}))); continue; end
            if(ischar(option.(fields{i})))
                option.(fields{i})=cellstr(option.(fields{i}));
            end
            if(~iscellstr(option.(fields{i})) || ...
                    numel(size(option.(fields{i})))~=2 || ...
                    ~any(size(option.(fields{i}),1)==[1 nrecs]))
                error('seizmo:changename:badInput',...
                    ['%s option must be a cellstr/char array\n'...
                    'with a single string or one for each record!'],...
                    fields{i});
            end
            if(any(cellfun('size',option.(fields{i}),1)>1))
                error('seizmo:changename:badInput',...
                    ['%s option does not allow cellstr arrays\n'...
                    'with multi-line cells!'],fields{i});
            end
            if(size(option.(fields{i}),1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
        case 'change'
            if(isempty(option.(fields{i}))); continue; end
            if(~iscellstr(option.(fields{i})) || ...
                    numel(size(option.(fields{i})))~=2 || ...
                    ~any(size(option.(fields{i}),1)==[1 nrecs]) || ...
                    mod(size(option.(fields{i}),2),2))
                error('seizmo:changename:badInput',...
                    ['%s option must be a cellstr array\n'...
                    'with input/output string pairs!'],fields{i});
            end
            if(any(cellfun('size',option.(fields{i}),1)>1))
                error('seizmo:changename:badInput',...
                    ['%s option does not allow cellstr arrays\n'...
                    'with multi-line cells!'],fields{i});
            end
            if(size(option.(fields{i}),1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
        otherwise
            error('seizmo:changename:badField',...
                'Unknown Option: %s',fields{i});
    end
end

% implement options
if(~isempty(option.NAME))
    [data.name]=deal(option.NAME{:});
end
if(~isempty(option.PREPEND))
    option.PREPEND=strcat(option.PREPEND,{data.name}');
    [data.name]=deal(option.PREPEND{:});
end
if(~isempty(option.APPEND))
    option.APPEND=strcat({data.name}',option.APPEND);
    [data.name]=deal(option.APPEND{:});
end
if(~isempty(option.DELETE))
    for i=1:nrecs
        for j=1:size(option.DELETE,2)
            data(i).name=strrep(data(i).name,option.DELETE{i,j},'');
        end
    end
end
if(~isempty(option.CHANGE))
    for i=1:nrecs
        for j=1:2:size(option.CHANGE,2)
            data(i).name=strrep(data(i).name,...
                option.CHANGE{i,j},option.CHANGE{i,j+1});
        end
    end
end

end
