function [data]=changepath(data,varargin)
%CHANGEPATH    Change the filepath of SEIZMO records
%
%    Usage:    data=changepath(data,...,'path',path,...)
%              data=changepath(data,...,'prepend',string,...)
%              data=changepath(data,...,'append',string,...)
%              data=changepath(data,...,'delete',string,...)
%              data=changepath(data,...,'delete',{string1 ... stringN},...)
%              data=changepath(data,...,'change',{orig replacement},...)
%
%    Description: CHANGEPATH(DATA,...,'PATH',PATH,...) sets the path field
%     of DATA to PATH.  The path field is the filepath associated with the
%     record.  If DATA contains more than 1 record and PATH is just a
%     single string, the path field for all records are changed to PATH.
%     Otherwise PATH must be a char/cellstr array with 1 row/element for
%     each record in DATA.
%
%     CHANGEPATH(DATA,...,'PREPEND',STRING,...) prepends STRING to the
%     path field of DATA.  If DATA contains more than 1 record and STRING
%     is just a single string, STRING is prepended to the path field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.
%
%     CHANGEPATH(DATA,...,'APPEND',STRING,...) appends STRING to the path
%     field of DATA.  If DATA contains more than 1 record and STRING is
%     just a single string, STRING is appended to the path field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.
%
%     CHANGEPATH(DATA,...,'DELETE',STRING,...) deletes STRING from the path
%     field of DATA.  If DATA contains more than 1 record and STRING is
%     just a single string, STRING is deleted from the path field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.  To delete multiple parts from
%     the path field, see the next DELETE option usage.
%
%     CHANGEPATH(DATA,...,'DELETE',{STRING1 ... STRINGN},...) allows
%     deleting multiple strings from the path field of DATA.  Multi-line
%     entries are not allowed (use one string per cell if working on
%     multiple records where each row corresponds to a separate record).  
%
%     CHANGEPATH(DATA,...,'CHANGE',{ORIGINAL REPLACEMENT},...) replaces all
%     occurances of the string ORIGINAL in the path field of DATA with the
%     string REPLACEMENT.  If DATA contains more than 1 record and ORIGINAL
%     and REPLACEMENT are just single strings, then the same string set is
%     applied to every record.  Otherwise a string pair must be given for
%     each record (requires a cellstr array).  To change multiple portions
%     of the path field for records in DATA, just expand the number of
%     columns in the cellstr array.
%
%    Notes:
%     - CHANGEPATH does NOT check the validity of input strings! 
%     - The order in which the options are implemented:
%        path, prepend, append, delete, change
%       This may come in handy for some more complicated cases.
%     - Using an option more than once replaces the previous option value
%       with the later one.
%
%    Examples:
%     Prepend and append to all filepaths:
%      data=changepath(data,'prepend','/some/dir/','append','/two/more')
%
%     Change the path of a record:
%      data(3)=changepath(data(3),'path','my/new/file/path')
%
%     Delete multiple parts from all filepaths:
%      data=changepath(data,'delete',{'/some/abs/path' '../pointless/..'})
%
%     Replace certain portions of all filepaths:
%      data=changepath(data,'change',{'//' '/' 'dir1' 'dir2'})
%
%    See also: changebyteorder, changename, writeheader, writeseizmo

%     Version History:
%        May  28, 2009 - initial version
%        May  29, 2009 - fix empty option handling, allow scalar expansion
%                        for path
%        June 12, 2009 - add testing table
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2009 at 17:30 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:changepath:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% fast exit
if(nargin==1); return; end

% default options
option.PATH=[];
option.PREPEND=[];
option.APPEND=[];
option.DELETE=[];
option.CHANGE=[];

% get options from SEIZMO global
global SEIZMO
try
    fields=fieldnames(SEIZMO.CHANGEPATH);
    for i=1:numel(fields)
        option.(fields{i})=SEIZMO.CHANGEPATH.(fields{i});
    end
catch
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:changepath:badInput',...
            'Options must be specified as a strings!');
    end
    option.(upper(varargin{i}))=varargin{i+1};
end

% check options
nrecs=numel(data);
fields=fieldnames(option);
for i=1:numel(fields)
    % skip empty fields
    if(isempty(option.(fields{i}))); continue; end
    % specific checks
    switch lower(fields{i})
        case {'path' 'append' 'prepend'}
            if(ischar(option.(fields{i})))
                option.(fields{i})=cellstr(option.(fields{i}));
            elseif(iscellstr(option.(fields{i})))
                % flatten multi-line paths
                option.(fields{i})=cellstr(char(option.(fields{i})));
            end
            if(~iscellstr(option.(fields{i})) || ...
                    numel(size(option.(fields{i})))~=2 || ...
                    ~any(size(option.(fields{i}),1)==[1 nrecs]))
                error('seizmo:changepath:badInput',...
                    ['%s option must be a cellstr/char array\n'...
                    'with a single string or one for each record!'],...
                    fields{i});
            end
        case 'delete'
            if(ischar(option.(fields{i})))
                option.(fields{i})=cellstr(option.(fields{i}));
            end
            if(~iscellstr(option.(fields{i})) || ...
                    numel(size(option.(fields{i})))~=2 || ...
                    ~any(size(option.(fields{i}),1)==[1 nrecs]))
                error('seizmo:changepath:badInput',...
                    ['%s option must be a cellstr/char array\n'...
                    'with a single string or one for each record!'],...
                    fields{i});
            end
            if(any(cellfun('size',option.(fields{i}),1)>1))
                error('seizmo:changepath:badInput',...
                    ['%s option does not allow cellstr arrays\n'...
                    'with multi-line cells!'],fields{i});
            end
            if(size(option.(fields{i}),1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
        case 'change'
            if(~iscellstr(option.(fields{i})) || ...
                    numel(size(option.(fields{i})))~=2 || ...
                    ~any(size(option.(fields{i}),1)==[1 nrecs]) || ...
                    mod(size(option.(fields{i}),2),2))
                error('seizmo:changepath:badInput',...
                    ['%s option must be a cellstr array\n'...
                    'with input/output string pairs!'],fields{i});
            end
            if(any(cellfun('size',option.(fields{i}),1)>1))
                error('seizmo:changepath:badInput',...
                    ['%s option does not allow cellstr arrays\n'...
                    'with multi-line cells!'],fields{i});
            end
            if(size(option.(fields{i}),1)==1)
                option.(fields{i})=option.(fields{i})(ones(nrecs,1),:);
            end
    end
end

% implement options
if(~isempty(option.PATH))
    [data.path]=deal(option.PATH{:});
end
if(~isempty(option.PREPEND))
    option.PREPEND=strcat(option.PREPEND,{data.path}');
    [data.path]=deal(option.PREPEND{:});
end
if(~isempty(option.APPEND))
    option.APPEND=strcat({data.path}',option.APPEND);
    [data.path]=deal(option.APPEND{:});
end
if(~isempty(option.DELETE))
    for i=1:nrecs
        for j=1:size(option.DELETE,2)
            data(i).path=strrep(data(i).path,option.DELETE{i,j},'');
        end
    end
end
if(~isempty(option.CHANGE))
    for i=1:nrecs
        for j=1:2:size(option.CHANGE,2)
            data(i).path=strrep(data(i).path,...
                option.CHANGE{i,j},option.CHANGE{i,j+1});
        end
    end
end

end
