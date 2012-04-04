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
%    Description:
%     DATA=CHANGEPATH(DATA,...,'PATH',PATH,...) sets the path field of DATA
%     to PATH.  The path field is the filepath associated with the record.
%     If DATA contains more than 1 record and PATH is just a single string,
%     the path field for all records are changed to PATH.  Otherwise PATH
%     must be a char/cellstr array with 1 row/element for each record in
%     DATA.
%
%     DATA=CHANGEPATH(DATA,...,'PREPEND',STRING,...) prepends STRING to the
%     path field of DATA.  If DATA contains more than 1 record and STRING
%     is just a single string, STRING is prepended to the path field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.
%
%     DATA=CHANGEPATH(DATA,...,'APPEND',STRING,...) appends STRING to the
%     path field of DATA.  If DATA contains more than 1 record and STRING
%     is just a single string, STRING is appended to the path field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.
%
%     DATA=CHANGEPATH(DATA,...,'DELETE',STRING,...) deletes STRING from the
%     path field of DATA.  If DATA contains more than 1 record and STRING
%     is just a single string, STRING is deleted from the path field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.  To delete multiple parts from
%     the path field, see the next DELETE option usage.
%
%     DATA=CHANGEPATH(DATA,...,'DELETE',{STRING1 ... STRINGN},...) allows
%     deleting multiple strings from the path field of DATA.  Multi-line
%     entries are not allowed (use one string per cell if working on
%     multiple records where each row corresponds to a separate record).  
%
%     DATA=CHANGEPATH(DATA,...,'CHANGE',{ORIGINAL REPLACEMENT},...) alters
%     all occurances of the string ORIGINAL in the path field of DATA with
%     the string REPLACEMENT.  If DATA contains more than 1 record and
%     ORIGINAL and REPLACEMENT are just single strings, then the same
%     string set is applied to every record.  Otherwise a string pair must
%     be given for each record (requires a cellstr array).  To change
%     multiple portions of the path field for records in DATA, just expand
%     the number of columns in the cellstr array.
%
%    Notes:
%     - CHANGEPATH does NOT check the validity of the file paths created!
%
%    Examples:
%     % Prepend and append to all filepaths:
%     data=changepath(data,'prepend','/some/dir/','append','/two/more')
%
%     % Change the path of a record:
%     data(3)=changepath(data(3),'path','my/new/file/path')
%
%     % Delete multiple parts from all filepaths:
%     data=changepath(data,'delete',{'/some/abs/path' '../pointless/..'})
%
%     % Replace certain portions of all filepaths:
%     data=changepath(data,'change',{'//' '/' 'dir1' 'dir2'})
%
%    See also: CHANGEBYTEORDER, CHANGENAME, WRITEHEADER, WRITESEIZMO

%     Version History:
%        May  28, 2009 - initial version
%        May  29, 2009 - fix empty option handling, allow scalar expansion
%                        for path
%        Sep.  5, 2009 - PATHNAME and NAME also point to PATH option, error
%                        on unknown option
%        Jan. 28, 2010 - drop global options, options are implemented in
%                        the order given, multiple calls for the same
%                        option are all done, seizmoverbose support
%        Apr. 10, 2010 - fix bug where path was written as a cellstr
%        Apr. 25, 2010 - allow several more strings to specify options
%        Feb. 11, 2011 - mass seizmocheck fix
%        Apr.  2, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 15:05 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:changepath:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data));

% fast exit
if(nargin==1); return; end

% number of records
nrecs=numel(data);

% detail message
if(seizmoverbose)
    disp('Changing Paths of Record(s)')
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:changepath:badInput',...
            'Options must be specified as strings!');
    end
    
    % check
    switch lower(varargin{i})
        case {'pathname' 'path' 'filepath' 'append' 'prepend'...
                'pathappend' 'appendpath' 'prependpath' 'pathprepend'}
            if(isempty(varargin{i+1})); continue; end
            if(ischar(varargin{i+1}))
                varargin{i+1}=cellstr(varargin{i+1});
            elseif(iscellstr(varargin{i+1}))
                % flatten multi-line names
                varargin{i+1}=cellstr(char(varargin{i+1}));
            end
            if(~iscellstr(varargin{i+1}) || ...
                    numel(size(varargin{i+1}))~=2 || ...
                    ~any(size(varargin{i+1},1)==[1 nrecs]))
                error('seizmo:changepath:badInput',...
                    ['%s option must be a cellstr/char array\n'...
                    'with a single string or one for each record!'],...
                    upper(varargin{i}));
            end
        case {'delete' 'pathdelete' 'deletepath'}
            if(isempty(varargin{i+1})); continue; end
            if(ischar(varargin{i+1}))
                varargin{i+1}=cellstr(varargin{i+1});
            end
            if(~iscellstr(varargin{i+1}) || ...
                    numel(size(varargin{i+1}))~=2 || ...
                    ~any(size(varargin{i+1},1)==[1 nrecs]))
                error('seizmo:changepath:badInput',...
                    ['%s option must be a cellstr/char array\n'...
                    'with a single string or one for each record!'],...
                    upper(varargin{i}));
            end
            if(any(cellfun('size',varargin{i+1},1)>1))
                error('seizmo:changepath:badInput',...
                    ['%s option does not allow cellstr arrays\n'...
                    'with multi-line cells!'],upper(varargin{i}));
            end
            if(size(varargin{i+1},1)==1)
                varargin{i+1}=varargin{i+1}(ones(nrecs,1),:);
            end
        case {'change' 'changepath' 'pathchange'}
            if(isempty(varargin{i+1})); continue; end
            if(~iscellstr(varargin{i+1}) || ...
                    numel(size(varargin{i+1}))~=2 || ...
                    ~any(size(varargin{i+1},1)==[1 nrecs]) || ...
                    mod(size(varargin{i+1},2),2))
                error('seizmo:changepath:badInput',...
                    ['%s option must be a cellstr array\n'...
                    'with input/output string pairs!'],upper(varargin{i}));
            end
            if(any(cellfun('size',varargin{i+1},1)>1))
                error('seizmo:changepath:badInput',...
                    ['%s option does not allow cellstr arrays\n'...
                    'with multi-line cells!'],upper(varargin{i}));
            end
            if(size(varargin{i+1},1)==1)
                varargin{i+1}=varargin{i+1}(ones(nrecs,1),:);
            end
        otherwise
            error('seizmo:changepath:badField',...
                'Unknown Option: %s',upper(varargin{i}));
    end
    
    % implement
    switch lower(varargin{i})
        case {'pathname' 'path' 'filepath'}
            [data.path]=deal(varargin{i+1}{:});
        case {'prepend' 'prependpath' 'pathprepend'}
            varargin{i+1}=strcat(varargin{i+1},{data.path}');
            [data.path]=deal(varargin{i+1}{:});
        case {'append' 'pathappend' 'appendpath'}
            varargin{i+1}=strcat({data.path}',varargin{i+1});
            [data.path]=deal(varargin{i+1}{:});
        case {'delete' 'pathdelete' 'deletepath'}
            for j=1:nrecs
                for k=1:size(varargin{i+1},2)
                    data(j).path=...
                        strrep(data(j).path,varargin{i+1}{j,k},'');
                end
            end
        case {'change' 'changepath' 'pathchange'}
            for j=1:nrecs
                for k=1:2:size(varargin{i+1},2)
                    data(j).path=strrep(data(j).path,...
                        varargin{i+1}{j,k},varargin{i+1}{j,k+1});
                end
            end
    end
end

end
