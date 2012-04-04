function [option]=cutparameters(varargin)
%CUTPARAMETERS    Parses inputs defining the data window(s)
%
%    Usage:    options=cutparameters(nrecs,args)
%
%    Description:
%     CUTPARAMETERS(NRECS,VARARGIN) parses parameters VARARGIN passed to
%     CUT and READDATAWINDOW, passing results back as a struct.  Parameters
%     include those that define the window as well as other options (fill,
%     cmplist, etc).  Parameters go through basic type/size checks (which
%     is why NRECS is an argument).
%
%    Notes:
%     - Empty input arguments may have unexpected results!
%
%    Examples:
%     % CUTPARAMETERS is what allows CUT and READDATAWINDOW
%     % to have very flexible input lists.
%
%    See also: CUT, READDATAWINDOW

%     Version History:
%        Apr. 17, 2008 - initial version
%        Apr. 18, 2008 - bugfix
%        June 24, 2008 - more input checks, major doc update
%        Sep. 14, 2008 - minor doc update, minor code cleaning
%        Oct.  2, 2008 - output now a struct, added CMPLIST, SEIZMO global
%                        access, code cleaning
%        Oct.  7, 2008 - improved CMPLIST checking, input change to check
%                        sizes
%        Oct. 16, 2008 - name changed from CUTPARAM to CUTPARAMETERS
%        Nov. 15, 2008 - update for new name schema
%        Apr. 23, 2009 - move usage up
%        May  29, 2009 - minor doc update, add nargin check
%        Oct.  6, 2009 - drop use of LOGICAL function
%        Aug. 21, 2010 - nargchk fix
%        Apr. 12, 2011 - allow pad/padding
%        Nov.  2, 2011 - doc update, allow 1x2/3/5/6 absolute time entries,
%                        fix ref with no offset bug
%        Dec.  1, 2011 - IZTYPE option added
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  1, 2011 at 20:45 GMT

% todo:

% check number of inputs
error(nargchk(1,inf,nargin));

% option defaults
option.TRIM=true;
option.FILL=false;
option.FILLER=0;
option.REF1='b';
option.REF2='e';
option.OFFSET1=0;
option.OFFSET2=0;
option.CMPLIST={':'};
option.IZTYPE=[];

% get options set by SEIZMO global
global SEIZMO; fields=fieldnames(option).';
if(isfield(SEIZMO,'CUTPARAMETERS'))
    for i=fields
        if(isfield(SEIZMO.CUTPARAMETERS,i{:})); 
            option.(i{:})=SEIZMO.CUTPARAMETERS.(i{:}); 
        end
    end
end

% number of records
nrecs=varargin{1};
varargin=varargin(2:end);

% read in and check inline cut parameters
skip=false; % skip next argument? (already read in)
ref1=false; % is REF1/OFFSET1 set?
ref2=false; % is REF2/OFFSET2 set?
for i=1:nargin-1
    % skip flag
    if(skip); skip=false; continue; end
    
    % skip empty
    if(isempty(varargin{i})); continue; end
    
    % offset without reference or an absolute time
    % NOTE: column vectors for offsets, row vectors for abs times
    if(isnumeric(varargin{i}))
        % offset or abs time?
        sz=size(varargin{i});
        if(sz(2)==1)
            % first offset (set reference to 'z')
            if(~ref1)
                option.REF1='z';
                option.OFFSET1=varargin{i};
                ref1=true;
            % second offset (use first reference)
            elseif(~ref2)
                option.REF2=option.REF1;
                option.OFFSET2=varargin{i};
                ref2=true;
            % only two offsets allowed
            else
                error('seizmo:cutparameters:badInput',...
                    'Too many window offsets!');
            end
        elseif(sz(1)==1) % abs time?
            if(~any(sz(2)==[2 3 5 6]))
                error('seizmo:cutparameters:badInput',...
                    'Malformed absolute time argument!');
            else % abs time
                % window reference
                if(~ref1)
                    option.REF1=varargin{i};
                    ref1=true;
                    
                    % quit if last arg
                    if(i==nargin-1); continue; end
                    
                    % first offset
                    sz=size(varargin{i+1});
                    if(isempty(varargin{i+1}))
                        % default offset (0)
                        skip=true;
                    elseif(isnumeric(varargin{i+1}))
                        if(sz(2)==1)
                            % offsets to abs time
                            option.OFFSET1=varargin{i+1};
                            skip=true;
                        elseif(sz(1)==1)
                            % abs time 2
                            skip=false;
                        else
                            error('seizmo:cutparameters:badInput',...
                                'Numeric inputs must be scalars/vectors!');
                        end
                    else
                        % string, offset not given (use default of 0)
                        skip=false;
                    end
                elseif(~ref2)
                    option.REF2=varargin{i};
                    ref2=true;
                    
                    % quit if last arg
                    if(i==nargin-1); continue; end
                    
                    % second offset
                    if(isempty(varargin{i+1}))
                        % default offset (0)
                        skip=true;
                    elseif(isnumeric(varargin{i+1}))
                        if(sz(2)==1)
                            % offsets to abs time
                            option.OFFSET2=varargin{i+1};
                            skip=true;
                        elseif(sz(1)==1)
                            % abs time 3?? only two references allowed!
                            error('seizmo:cutparameters:badInput',...
                                'Too many window references!');
                        else
                            error('seizmo:cutparameters:badInput',...
                                'Numeric inputs must be scalars/vectors!');
                        end
                    else
                        % string, offset not given (use default of 0)
                        skip=false;
                    end
                else
                    % only two references allowed
                    error('seizmo:cutparameters:badInput',...
                        'Too many window references!');
                end
            end
        else
            error('seizmo:cutparameters:badInput',...
                    'Numeric arguments must be scalars or vectors!');
        end
    % string (could be reference or other option)
    elseif(ischar(varargin{i}))
        skip=true;
        
        % checks
        if(any(strcmpi({'fill' 'filler' 'trim'},varargin{i})))
            if(isempty(varargin{i+1}))
                continue; % leave as default/current by skipping
            elseif(~isnumeric(varargin{i+1}) && ...
                    ~islogical(varargin{i+1}))
                error('seizmo:cutparameters:badInput',...
                    'OPTION value must be numeric or logical!');
            end
        end
        
        % options
        switch lower(varargin{i})
            case {'fill' 'pad'}
                option.FILL=varargin{i+1}~=0;
            case {'filler' 'padding'}
                option.FILLER=double(varargin{i+1});
            case 'trim'
                if(~isscalar(varargin{i+1}))
                    error('seizmo:cutparameters:badInput',...
                        'TRIM option must be scalar!');
                end
                option.TRIM=varargin{i+1}~=0;
            case 'iztype'
                if(isempty(varargin{i+1})); continue; end
                if(ischar(varargin{i+1}))
                    varargin{i+1}=cellstr(varargin{i+1});
                end
                if(~iscellstr(varargin{i+1}))
                    error('seizmo:cutparameters:badInput',...
                        'IZTYPE option must be a string!');
                end
                option.IZTYPE=varargin{i+1};
            case 'cmplist'
                if(iscell(varargin{i+1}))
                    for j=1:numel(varargin{i+1})
                        if(~isnumeric(varargin{i+1}{j}) ...
                            && ~strcmp(varargin{i+1}{j},':'))
                            error('seizmo:cutparameters:badInput',...
                                'CMPLIST incomprehensible!');
                        end
                    end
                    option.CMPLIST=varargin{i+1};
                elseif(isnumeric(varargin{i+1}) ...
                        && all(varargin{i+1}==fix(varargin{i+1})) ...
                        && all(varargin{i+1}>0))
                    option.CMPLIST=num2cell(varargin{i+1},2);
                elseif(ischar(varargin{i+1}) && strcmp(varargin{i+1},':'))
                    option.CMPLIST=cellstr(varargin{i+1});
                else
                    error('seizmo:cutparameters:badInput',...
                        'CMPLIST incomprehensible!');
                end
            otherwise
                % window reference
                if(~ref1)
                    option.REF1=varargin{i};
                    ref1=true;
                    
                    % quit if last arg
                    if(i==nargin-1); continue; end
                    
                    % first offset
                    sz=size(varargin{i+1});
                    if(isempty(varargin{i+1}))
                        continue;
                    elseif(isnumeric(varargin{i+1}))
                        if(sz(2)==1)
                            % offsets to abs time
                            option.OFFSET1=varargin{i+1};
                            skip=true;
                        elseif(sz(1)==1)
                            % abs time 2
                            skip=false;
                        else
                            error('seizmo:cutparameters:badInput',...
                                'Numeric inputs must be scalars/vectors!');
                        end
                    else
                        % offset not given
                        skip=false;
                    end
                elseif(~ref2)
                    option.REF2=varargin{i};
                    ref2=true;
                    
                    % quit if last arg
                    if(i==nargin-1); continue; end
                    
                    % second offset
                    sz=size(varargin{i+1});
                    if(isempty(varargin{i+1}))
                        continue;
                    elseif(isnumeric(varargin{i+1}))
                        if(sz(2)==1)
                            % offsets to abs time
                            option.OFFSET2=varargin{i+1};
                            skip=true;
                        elseif(sz(1)==1)
                            % abs time 3?? only two references allowed!
                            error('seizmo:cutparameters:badInput',...
                                'Too many window references!');
                        else
                            error('seizmo:cutparameters:badInput',...
                                'Numeric inputs must be scalars/vectors!');
                        end
                    else
                        % offset not given
                        skip=false;
                    end
                else
                    % only two references allowed
                    error('seizmo:cutparameters:badInput',...
                        'Too many window references!');
                end
        end
    else
        error('seizmo:cutparameters:badInput',...
            'Window parameters must be strings or numeric!');
    end
end

% expand scalars
if(numel(option.OFFSET1)==1); option.OFFSET1(1:nrecs,1)=option.OFFSET1; end
if(numel(option.OFFSET2)==1); option.OFFSET2(1:nrecs,1)=option.OFFSET2; end
if(numel(option.CMPLIST)==1); option.CMPLIST(1:nrecs,1)=option.CMPLIST; end
if(numel(option.FILLER)==1); option.FILLER(1:nrecs,1)=option.FILLER; end
if(numel(option.FILL)==1); option.FILL(1:nrecs,1)=option.FILL; end
if(numel(option.IZTYPE)==1); option.IZTYPE(1:nrecs,1)=option.IZTYPE; end

% assure column vectors
option.OFFSET1=option.OFFSET1(:);
option.OFFSET2=option.OFFSET2(:);
option.CMPLIST=option.CMPLIST(:);
option.FILLER=option.FILLER(:);
option.FILL=option.FILL(:);
option.IZTYPE=option.IZTYPE(:);

% cut parameter checks
if(numel(option.OFFSET1)~=nrecs || numel(option.OFFSET2)~=nrecs)
    error('seizmo:cutparameters:badInputSize',...
        'Number of elements for OFFSET not correct!');
elseif(numel(option.CMPLIST)~=nrecs)
    error('seizmo:cutparameters:badInputSize',...
        'Number of elements for CMPLIST option not correct!');
elseif(numel(option.FILLER)~=nrecs)
    error('seizmo:cutparameters:badInputSize',...
        'Number of elements for FILLER option not correct!');
elseif(numel(option.FILL)~=nrecs)
    error('seizmo:cutparameters:badInputSize',...
        'Number of elements for FILL option not correct!');
elseif(~isempty(option.IZTYPE) && numel(option.IZTYPE)~=nrecs)
    error('seizmo:cutparameters:badInputSize',...
        'Number of elements for IZTYPE option not correct!');
end

end
