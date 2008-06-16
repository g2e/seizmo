function []=lh(data,varargin)
%LH    List SAClab data headers
%
%    Description: LH(DATA,FIELD1,FIELD2,...) lists the specified header 
%     field(s) FIELD1, FIELD2, ... and their value(s) in DATA in a manner 
%     similar to SAC's lh formating.  FIELDS must be strings corresponding 
%     to a valid header field or a valid group field (ie. t, kt, resp, 
%     user, kuser).
%
%    Usage:    lh(data,'field1','field2',...)
%
%    Examples:
%      lh(data)          % lists all header variables
%      lh(data,'t')      % lists t group
%
%     Fields are case independent:
%      lh(data,'dEltA')
%      lh(data,'StLA','stLo')
%
%    See also:  ch, gh, glgc, genum, genumdesc

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - Documentation Update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 16:25 GMT

% require at least one input
if(nargin<1)
    error('MATLAB:nargchk:notEnoughInputs',...
        'Not enough input arguments.')
end

% check data structure
error(seischk(data))

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=seisdef(vers(nver));
for i=1:nver-1
    h(i)=seisdef(vers(i));
end

% extra padding
disp(' ')

% loop over files
for i=1:length(data)
    % get header info individually
    v=data(i).version==vers;
    
    % list all case
    if (nargin==1)
        clear varargin
        % all available field type sets
        for j=1:length(h(v).types)
            for k=1:length(h(v).(h(v).types{j}))
                varargin{k,j}=fieldnames(h(v).(h(v).types{j})(k).pos)';
            end
        end
        varargin=[varargin{:}].';
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
    disp('---------------------------')
    
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
        
        % group loop
        for k=1:glen
            % modify field if group
            if(group); f=[gf num2str(g(k))];
            else f=gf;
            end
            
            % display field/value
            lh_disp(h(v),f,data(i));
        end
    end
end

end

function []=lh_disp(h,f,data)
%LH_DISP    Finds and displays field/value pairs for lh

% handle enum/lgc fields specially
for m=1:length(h.enum)
    if(isfield(h.enum(m).pos,f))
        ival=data.head(h.enum(m).pos.(f));
        % searches only enum(1) for now (prefer one big
        % set of enums to share vs separate sets)
        if(ival==fix(ival) && ival>=h.enum(1).minval && ...
                ival<=h.enum(1).maxval)
            disp(sprintf('%12s = %s',upper(f),h.enum(1).desc{ival+1}));
        elseif(ival==h.undef.ntype)
            disp(sprintf('%12s = UNDEFINED (%d)',upper(f),ival));
        else
            disp(sprintf('%12s = UNKNOWN (%d)',upper(f),ival));
        end
        return;
    end
end
for m=1:length(h.lgc)
    if(isfield(h.lgc(m).pos,f))
        if(data.head(h.lgc(m).pos.(f))==h.false)
            disp(sprintf('%12s = FALSE',upper(f)));
        elseif(data.head(h.lgc(m).pos.(f))==h.true)
            disp(sprintf('%12s = TRUE',upper(f)));  
        elseif(data.head(h.lgc(m).pos.(f))==h.undef.ntype)
            disp(sprintf('%12s = UNDEFINED (%d)',upper(f),...
                data.head(h.lgc(m).pos.(f))));
        else
            disp(sprintf('%12s = INVALID (%d)',upper(f),...
                data.head(h.lgc(m).pos.(f))));
        end
        return;
    end
end

% string types
for m=1:length(h.stype)
    for n=1:length(h.(h.stype{m}))
        if(isfield(h.(h.stype{m})(n).pos,f))
            p=h.(h.stype{m})(n).pos.(f);
            q=p(2)-p(1)+1; u=length(h.undef.stype);
            if(strcmp([h.undef.stype ones(1,q-u)*32],...
                char(data.head(p(1):p(2)).')))
                disp(sprintf('%12s = UNDEFINED (%s)',upper(f),...
                    char(data.head(p(1):p(2)))));
            else
                disp(sprintf('%12s = %s',upper(f),...
                    char(data.head(p(1):p(2)))));
            end
            return;
        end
    end
end

% remaining numeric types
for m=1:length(h.ntype)
    for n=1:length(h.(h.ntype{m}))
        if(isfield(h.(h.ntype{m})(n).pos,f))
            if(data.head(h.(h.ntype{m})(n).pos.(f))==h.undef.ntype)
                disp(sprintf('%12s = UNDEFINED (%d)',upper(f),...
                    data.head(h.(h.ntype{m})(n).pos.(f))));
            else
                disp(sprintf('%12s = %d',upper(f),...
                    data.head(h.(h.ntype{m})(n).pos.(f))));
            end
            return;
        end
    end
end

% field not found
warning('SAClab:lh:fieldInvalid',...
    '\nHeader Version: %d\nInvalid field: %s',h.version,f);

end
