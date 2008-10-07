function [option]=cutparam(varargin)
%CUTPARAM    Parses inputs defining the cut window for CUTIM/RPDW
%
%    Description: CUTPARAM(NRECS,VARARGIN) parses parameters VARARGIN
%     passed to CUTIM and RPDW, passing the results back as a struct.  
%     Parameters include those that define the window as well as other
%     options (fill, cmplist, etc).  Parameters go through basic type/size
%     checks (which is why NRECS is an argument).
%
%    Notes:
%     - empty input arguments may have unexpected results
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage:    options=cutparam(nrecs,args)
%
%    Examples: NONE
%
%    See also: cutim, rpdw

%     Version History:
%        Apr. 17, 2008 - initial version
%        Apr. 18, 2008 - bugfix
%        June 24, 2008 - more input checks, major doc update
%        Sep. 14, 2008 - minor doc update, minor code cleaning
%        Oct.  2, 2008 - output now a struct, added CMPLIST, SACLAB global
%                        access, code cleaning
%        Oct.  7, 2008 - improved CMPLIST checking, input change to check
%                        sizes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  7, 2008 at 06:50 GMT

% todo:

% option defaults
option.TRIM=true;
option.FILL=false;
option.FILLER=0;
option.REF1='b';
option.REF2='e';
option.OFFSET1=0;
option.OFFSET2=0;
option.CMPLIST=':';

% get options set by SACLAB global
global SACLAB; fields=fieldnames(option).';
if(isfield(SACLAB,'CUTPARAM'))
    for i=fields
        if(isfield(SACLAB.CUTPARAM,i{:})); 
            option.(i{:})=SACLAB.CUTPARAM.(i{:}); 
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
for i=1:nargin
    % skip flag
    if(skip); skip=false; continue; end
    
    % skip empty
    if(isempty(varargin{i})); continue; end
    
    % offset without reference
    if(isnumeric(varargin{i}))
        % first offset (set reference to 'z')
        if(~ref1)
            option.REF1='z';
            option.OFFSET1=varargin{i};
            ref1=true;
        % second offset (use first reference)
        elseif(~ref2)
            option.REF2=ref1;
            option.OFFSET2=varargin{i};
            ref2=true;
        % only two offsets allowed
        else
            error('SAClab:cutparam:badInput','Too many window offsets!')
        end
    % string (could be reference or other option)
    elseif(ischar(varargin{i}))
        skip=true;
        
        % checks
        if(any(strcmpi({'fill' 'filler' 'trim'},varargin{i})))
            if(isempty(varargin{i+1}))
                continue; % leave as default by skipping
            elseif(~isnumeric(varargin{i+1}) && ...
                    ~islogical(varargin{i+1}))
                error('SAClab:cutparam:badInput',...
                    'OPTION value must be numeric or logical!')
            end
        end
        
        % options
        switch lower(varargin{i})
            case 'fill'
                option.FILL=logical(varargin{i+1});
            case 'filler'
                option.FILLER=double(varargin{i+1});
            case 'trim'
                if(~isscalar(varargin{i+1}))
                    error('SAClab:cutparam:badInput',...
                        'TRIM option must be scalar!');
                end
                option.TRIM=logical(varargin{i+1});
            case 'cmplist'
                if(iscell(varargin{i+1}))
                    for j=1:numel(varargin{i+1})
                        if(~isnumeric(varargin{i+1}{j}) ...
                            && ~strcmp(varargin{i+1}{j},':'))
                            error('SAClab:cutparam:badInput',...
                                'CMPLIST incomprehensible!');
                        end
                    end
                elseif(~isnumeric(varargin{i+1}) ...
                    && ~strcmp(varargin{i+1},':'))
                    error('SAClab:cutparam:badInput',...
                        'CMPLIST incomprehensible!');
                end
                option.CMPLIST=varargin{i+1};
            otherwise
                % window reference
                if(~ref1)
                    option.REF1=varargin{i};
                    ref1=true;
                    
                    % first offset
                    if(isempty(varargin{i+1}))
                        continue;
                    elseif(isnumeric(varargin{i+1}))
                        option.OFFSET1=varargin{i+1};
                    else
                        % offset not given
                        option.OFFSET1=0;
                        skip=false;
                    end
                elseif(~ref2)
                    option.REF2=varargin{i};
                    ref2=true;
                    
                    % second offset
                    if(isempty(varargin{i+1}))
                        continue;
                    elseif(isnumeric(varargin{i+1}))
                        option.OFFSET2=varargin{i+1};
                    else
                        % offset not given
                        option.OFFSET2=0;
                        skip=false;
                    end
                else
                    % only two references allowed
                    error('SAClab:cutparam:badInput',...
                        'Too many window references!')
                end
        end
    else
        error('SAClab:cutparam:badInput',...
            'Window parameters must be strings or numeric!')
    end
end

% expand scalars
if(numel(option.OFFSET1)==1); option.OFFSET1(1:nrecs,1)=option.OFFSET1; end
if(numel(option.OFFSET2)==1); option.OFFSET2(1:nrecs,1)=option.OFFSET2; end
if(numel(option.CMPLIST)==1); option.CMPLIST(1:nrecs,1)=option.CMPLIST; end
if(numel(option.FILLER)==1); option.FILLER(1:nrecs,1)=option.FILLER; end
if(numel(option.FILL)==1); option.FILL(1:nrecs,1)=option.FILL; end

% cut parameter checks
if(numel(option.OFFSET1)~=nrecs || numel(option.OFFSET2)~=nrecs)
    error('SAClab:cutparam:badInputSize',...
        'Number of elements in OFFSET not correct!')
elseif(numel(option.CMPLIST)~=nrecs)
    error('SAClab:cutparam:badInputSize',...
        'Number of elements in CMPLIST not correct!')
elseif(numel(option.FILLER)~=nrecs)
    error('SAClab:cutparam:badInputSize',...
        'Number of elements in FILLER not correct!')
elseif(numel(option.FILL)~=nrecs)
    error('SAClab:cutparam:badInputSize',...
        'Number of elements in FILL not correct!')
end

end
