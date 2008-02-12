function []=lh(data,varargin)
%LH    List SAClab headers
%
%    Description: List specified header field(s) and their value(s) of
%     SAClab data structure in a manner similar to SAC's lh display format.
%     Fields must be strings corresponding to a valid header field or a
%     valid group field (ie. t,kt,resp,user,kuser).
%
%    Usage:    lh(saclab_struct,'field1','field2',...)
%
%    Examples:
%     lh(data)          % lists all header variables
%     lh(data,'t')      % lists t group
%     lh(data,'dEltA')
%     lh(data,'StLA','stLo')
%
%    by Garrett Euler (2/2008)  ggeuler@wustl.edu
%
%    See also:  ch, gh, rh, wh, rpdw, rdata, rsac, bsac, wsac, sachp, gv

% do nothing on no input
if(nargin<1); return; end

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head'))
    error('data structure does not have proper fields')
end

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% check header field types match up to what lh currently supports
lhtypes={'real' 'int' 'enum' 'lgc' 'char'};
for i=1:nver
    for j=1:length(h(i).types)
        if(~any(strcmp(h(i).types{j},lhtypes)))
            warning('SAClab:missingTypes',...
                'Header definition has unsupported field type, %s',...
                h(i).types{j});
        end
    end
end

% loop over files
for i=1:length(data)
    % get header info individually
    v=data(i).version==vers;
    
    % check for list all
    if (nargin==1)
        clear varargin
        % all available field type sets
        for j=1:length(h(v).types)
            for k=1:length(h(v).(h(v).types{j}))
                varargin{k,j}=fieldnames(h(v).(h(v).types{j})(k).pos)';
            end
        end
        varargin=[varargin{:}]';
    end
    
    % get filename if available
    if(isfield(data,'name'))
        name=data(i).name;
    else
        name='';
    end
    
    % formatted header
    disp(' ')
    disp(sprintf(' FILE: %s - %d',name,i))
    disp('+--------------+--------------+')
    disp('| FIELD        | STORAGE      | VALUE               ')
    disp('+--------------+--------------+---------------------')
    
    % loop over fields
    for j=1:length(varargin)
        % force lowercase
        gf=lower(varargin{j});
        
        % check for group fields
        group=0; glen=1;
        if(isfield(h(v).grp,gf))
            g=h(v).grp.(gf).min:h(v).grp.(gf).max;
            group=1; glen=length(g);
        end
        
        % expansion loop
        for k=1:glen
            % modify field if group
            if(group); f=[gf num2str(g(k))];
            else f=gf;
            end
            
            % valid field flag
            set=0;
            
            % format by field type
            for m=1:length(h(v).real)
                if(isfield(h(v).real(m).pos,f))
                    if(data(i).head(h(v).real(m).pos.(f))==h(v).real(m).undef)
                        disp(sprintf('| %-12s | %-12s | UNDEFINED (%d)',...
                            upper(f),upper(h(v).real(m).store),...
                            data(i).head(h(v).real(m).pos.(f))));
                    else
                        disp(sprintf('| %-12s | %-12s | %d',upper(f),...
                            upper(h(v).real(m).store),...
                            data(i).head(h(v).real(m).pos.(f))));
                    end
                    set=1;
                end
            end
            for m=1:length(h(v).int)
                if(isfield(h(v).int(m).pos,f))
                    if(data(i).head(h(v).int(m).pos.(f))==h(v).int(m).undef)
                        disp(sprintf('| %-12s | %-12s | UNDEFINED (%d)',...
                            upper(f),upper(h(v).int(m).store),...
                            data(i).head(h(v).int(m).pos.(f))));
                    else
                        disp(sprintf('| %-12s | %-12s | %d',upper(f),...
                            upper(h(v).int(m).store),...
                            data(i).head(h(v).int(m).pos.(f))));
                    end
                    set=1;
                end
            end
            for m=1:length(h(v).enum)
                if(isfield(h(v).enum(m).pos,f))
                    ival=data(i).head(h(v).enum(m).pos.(f));
                    % searches only enum(1) for now (prefer one big
                    % set of enums to share vs separate sets)
                    if(ival-round(ival)==0 && ival>=h(v).enum(1).minval && ...
                            ival<=h(v).enum(1).maxval)
                        disp(sprintf('| %-12s | %-12s | %s',upper(f),...
                            upper(h(v).enum(m).store),...
                            h(v).enum(1).desc{ival+1}));
                    elseif(ival==h(v).enum(m).undef)
                        disp(sprintf('| %-12s | %-12s | UNDEFINED (%d)',...
                            upper(f),upper(h(v).enum(m).store),ival));
                    else
                        disp(sprintf('| %-12s | %-12s | UNKNOWN (%d)',...
                            upper(f),upper(h(v).enum(m).store),ival));
                    end
                    set=1;
                end
            end
            for m=1:length(h(v).lgc)
                if(isfield(h(v).lgc(m).pos,f))
                    if(data(i).head(h(v).lgc(m).pos.(f))==h(v).false)
                        disp(sprintf('| %-12s | %-12s | FALSE',upper(f),...
                            upper(h(v).lgc(m).store)));
                    elseif(data(i).head(h(v).lgc(m).pos.(f))==h(v).true)
                        disp(sprintf('| %-12s | %-12s | TRUE',upper(f),...
                            upper(h(v).lgc(m).store)));  
                    elseif(data(i).head(h(v).lgc(m).pos.(f))==h(v).lgc(m).undef)
                        disp(sprintf('| %-12s | %-12s | UNDEFINED (%d)',...
                            upper(f),upper(h(v).lgc(m).store),...
                            data(i).head(h(v).lgc(m).pos.(f))));
                    else
                        disp(sprintf('| %-12s | %-12s | INVALID (%d)',...
                            upper(f),upper(h(v).lgc(m).store),...
                            data(i).head(h(v).lgc(m).pos.(f))));
                    end
                    set=1;
                end
            end
            for m=1:length(h(v).char)
                if(isfield(h(v).char(m).pos,f))
                    n=h(v).char(m).pos.(f);
                    o=n(2)-n(1)+1; p=length(h(v).char(m).undef);
                    if(p>o) % only happens if header definition is messed
                        warning('SAClab:undefOutOfRange',...
                            ['Length of UNDEFINED string value exceeds'...
                            'maximum string length for field %s'],f)
                    elseif(strcmp([h(v).char(m).undef ones(1,o-p)*32],...
                            char(data(i).head(n(1):n(2))')))
                        disp(sprintf('| %-12s | %-12s | UNDEFINED (%s)',...
                            upper(f),upper(h(v).char(m).store),...
                            char(data(i).head(n(1):n(2)))));
                    else
                        disp(sprintf('| %-12s | %-12s | %s',...
                            upper(f),upper(h(v).char(m).store),...
                            char(data(i).head(n(1):n(2)))));
                    end
                    set=1;
                end
            end
            
            % check field was found
            if(~set)
                warning('SAClab:fieldInvalid','Invalid field: %s',f);
            end
        end
    end
    
    % formatted footer    
    disp('+--------------+--------------+---------------------')
    disp(' ')
end

end

