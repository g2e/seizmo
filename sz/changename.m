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
%    Description:
%     DATA=CHANGENAME(DATA,...,'NAME',NAME,...) sets the name field of DATA
%     to NAME.  The name field is the filename associated with the record.
%     If DATA contains more than 1 record and NAME is just a single string,
%     the name field for all records are changed to NAME.  Otherwise NAME
%     must be a char/cellstr array with 1 row/element for each record in
%     DATA.
%
%     DATA=CHANGENAME(DATA,...,'PREPEND',STRING,...) prepends STRING to the
%     name field of DATA.  If DATA contains more than 1 record and STRING
%     is just a single string, STRING is prepended to the name field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.
%
%     DATA=CHANGENAME(DATA,...,'APPEND',STRING,...) appends STRING to the
%     name field of DATA.  If DATA contains more than 1 record and STRING
%     is just a single string, STRING is appended to the name field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.
%
%     DATA=CHANGENAME(DATA,...,'DELETE',STRING,...) deletes STRING from the
%     name field of DATA.  If DATA contains more than 1 record and STRING
%     is just a single string, STRING is deleted from the name field of all
%     records in DATA.  Otherwise STRING must be a char/cellstr array with
%     1 row/element for each record in DATA.  To delete multiple parts from
%     the name field, see the next DELETE option usage.
%
%     DATA=CHANGENAME(DATA,...,'DELETE',{STRING1 ... STRINGN},...) allows
%     deleting multiple strings from the name field of DATA.  Multi-line
%     entries are not allowed (use one string per cell if working on
%     multiple records where each row corresponds to a separate record).  
%
%     DATA=CHANGENAME(DATA,...,'CHANGE',{ORIGINAL REPLACEMENT},...) sets
%     occurances of the string ORIGINAL in the name field of DATA to the
%     string REPLACEMENT.  If DATA contains more than 1 record and ORIGINAL
%     and REPLACEMENT are just single strings, then the same string set is
%     applied to every record.  Otherwise a string pair must be given for
%     each record (requires a cellstr array).  To change multiple portions
%     of the name field for records in DATA, just expand the number of
%     columns in the cellstr array.
%
%    Notes:
%     - CHANGENAME does NOT check the validity of the filenames created!
%
%    Examples:
%     % Prepend and append to all filenames:
%     data=changename(data,'prepend','zzz','append','.new')
%
%     % Change the name of a record:
%     data(3)=changename(data(3),'name','mynewfilename')
%
%     % Delete multiple parts from all filenames:
%     data=changename(data,'delete',{'__' '.merged' 'SAC'})
%
%     % Replace certain portions of all filenames:
%     data=changename(data,'change',{'..' '.__.' 'part1' 'part2'})
%
%    See also: CHANGEBYTEORDER, CHANGEPATH, WRITEHEADER, WRITESEIZMO,
%              WRITEPARAMETERS

%     Version History:
%        May  28, 2009 - initial version
%        May  29, 2009 - fix empty option handling, allow scalar expansion
%                        for name
%        Sep.  5, 2009 - FILENAME and FILE also point to NAME option, error
%                        on unknown option
%        Jan. 28, 2010 - drop global options, options are implemented in
%                        the order given, multiple calls for the same
%                        option are all done, seizmoverbose support
%        Apr. 10, 2010 - fixed bug where name was written as a cellstr
%        Apr. 25, 2010 - allow several more strings to specify options
%        Feb. 11, 2011 - mass seizmocheck fix
%        Feb.  7, 2012 - doc update
%        Aug. 30, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  7, 2012 at 15:05 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:changename:badNumInputs',...
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
    disp('Changing Names of Record(s)')
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error('seizmo:changename:badInput',...
            'Options must be specified as strings!');
    end
    
    % check
    switch lower(varargin{i})
        case {'filename' 'file' 'name' 'namefile' 'append' 'prepend'...
                'nameappend' 'appendname' 'prependname' 'nameprepend'}
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
                error('seizmo:changename:badInput',...
                    ['%s option must be a cellstr/char array\n'...
                    'with a single string or one for each record!'],...
                    upper(varargin{i}));
            end
        case {'delete' 'namedelete' 'deletename'}
            if(isempty(varargin{i+1})); continue; end
            if(ischar(varargin{i+1}))
                varargin{i+1}=cellstr(varargin{i+1});
            end
            if(~iscellstr(varargin{i+1}) || ...
                    numel(size(varargin{i+1}))~=2 || ...
                    ~any(size(varargin{i+1},1)==[1 nrecs]))
                error('seizmo:changename:badInput',...
                    ['%s option must be a cellstr/char array\n'...
                    'with a single string or one for each record!'],...
                    upper(varargin{i}));
            end
            if(any(cellfun('size',varargin{i+1},1)>1))
                error('seizmo:changename:badInput',...
                    ['%s option does not allow cellstr arrays\n'...
                    'with multi-line cells!'],upper(varargin{i}));
            end
            if(size(varargin{i+1},1)==1)
                varargin{i+1}=varargin{i+1}(ones(nrecs,1),:);
            end
        case {'change' 'namechange' 'changename'}
            if(isempty(varargin{i+1})); continue; end
            if(~iscellstr(varargin{i+1}) || ...
                    numel(size(varargin{i+1}))~=2 || ...
                    ~any(size(varargin{i+1},1)==[1 nrecs]) || ...
                    mod(size(varargin{i+1},2),2))
                error('seizmo:changename:badInput',...
                    ['%s option must be a cellstr array\n'...
                    'with input/output string pairs!'],upper(varargin{i}));
            end
            if(any(cellfun('size',varargin{i+1},1)>1))
                error('seizmo:changename:badInput',...
                    ['%s option does not allow cellstr arrays\n'...
                    'with multi-line cells!'],upper(varargin{i}));
            end
            if(size(varargin{i+1},1)==1)
                varargin{i+1}=varargin{i+1}(ones(nrecs,1),:);
            end
        otherwise
            error('seizmo:changename:badField',...
                'Unknown Option: %s',upper(varargin{i}));
    end
    
    % implement
    switch lower(varargin{i})
        case {'filename' 'file' 'name' 'namefile'}
            [data.name]=deal(varargin{i+1}{:});
        case {'prepend' 'prependname' 'nameprepend'}
            varargin{i+1}=strcat(varargin{i+1},{data.name}');
            [data.name]=deal(varargin{i+1}{:});
        case {'append' 'appendname' 'nameappend'}
            varargin{i+1}=strcat({data.name}',varargin{i+1});
            [data.name]=deal(varargin{i+1}{:});
        case {'delete' 'deletename' 'namedelete'}
            for j=1:nrecs
                for k=1:size(varargin{i+1},2)
                    data(j).name=...
                        strrep(data(j).name,varargin{i+1}{j,k},'');
                end
            end
        case {'change' 'changename' 'namechange'}
            for j=1:nrecs
                for k=1:2:size(varargin{i+1},2)
                    data(j).name=strrep(data(j).name,...
                        varargin{i+1}{j,k},varargin{i+1}{j,k+1});
                end
            end
    end
end

end
