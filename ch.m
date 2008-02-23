function [data]=ch(data,varargin)
%CH    Change SAClab headers
%
%    Description: Changes the specified SAClab header field(s) to the
%     specified value(s).  The field variable must be a string
%     corresponding to a valid header field.  Values may be scalar
%     (assigns same value to all) or vectors of length equal to the number
%     of records.  Values may be contained in a numeric, char, or cell
%     array.  Char arrays must be arranged as a column vector with one
%     string value per row.  Cell arrays must have one value per cell.
%
%    Usage:    new_struct=ch(saclab_struct,'group_field1',value1,...
%                               'field2',{'value2a' 'value2b' ...},...
%                               'field3',[value3a; value3b; ...])
%
%    Examples:
%     data=ch(data,'DELTA',gh(data,'delta')*2);
%     data=ch(data,'STLA',lats,'STLO',lons)
%     data=ch(data,'KT0','sSKS');
%
%    See also:  lh, gh, rh, wh, glgc, genum, genumdesc, 
%               rpdw, rdata, rsac, bsac, wsac, sachi, gv

% throw error if unpaired fields
if (mod(nargin-1,2)); error('Unpaired Field/Value!'); end

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head'))
    error('data structure does not have proper fields')
end

% number of records
nrecs=length(data);

% recursive section (breaks up data so ch only deals with 1 header version)
v=[data.version];
vers=unique(v(:));
nver=length(vers);
if(nver>1)
    for i=1:nver
        % need a way to parse varargin quickly
        temp=cell(1,nargin-1);
        for j=1:2:nargin-2
            temp{j}=varargin{j};
        end
        for j=2:2:nargin-1
            if(isscalar(varargin{j}) || isempty(varargin{j}))
                temp{j}=varargin{j};
            elseif(isvector(varargin{j}))
                if(length(varargin{j})~=nrecs)
                    error('Value vector for field %s not correct size',...
                        varargin{j-1});
                end
                temp{j}=varargin{j}(vers(i)==v);
            else
                if(size(varargin{j},1)~=nrecs)
                    error('Value array for field %s not correct size',...
                        varargin{j-1});
                end
                temp{j}=varargin{j}(vers(i)==v,:);
            end
        end
        data(vers(i)==v)=ch(data(vers(i)==v),temp{:});
    end
    return;
end

% headers setup
h=sachi(vers);

% check header field types match up to what ch currently supports
chtypes={'real' 'int' 'enum' 'lgc' 'char'};
for j=1:length(chtypes)
    if(~any(strcmp(h.types,chtypes{j})))
        error('ch field type %s unsupported by header vers %d',...
            chtypes{j},vers);
    end
end
for j=1:length(h.types)
    if(~any(strcmp(h.types{j},chtypes)))
        error('Header version %d has unsupported field type %s',...
            vers,h.types{j});
    end
end

% pull entire header
head=[data.head];

% loop over field/value pairs
for i=1:2:(nargin-2)
    % force values into cell array
    if(isnumeric(varargin{i+1}))
        varargin{i+1}=num2cell(varargin{i+1});
    elseif(ischar(varargin{i+1}))
        varargin{i+1}=cellstr(varargin{i+1});
    end
    
    % force field name to be lowercase
    gf=lower(varargin{i});
    
    % check for group fields
    group=0; glen=1;
    if(isfield(h.grp,gf))
        g=h.grp.(gf).min:h.grp.(gf).max;
        group=1; glen=length(g);
    end
    
    % check & expand values
    if(isscalar(varargin{i+1}))
        varargin{i+1}=varargin{i+1}(ones(nrecs,glen));
    elseif(isvector(varargin{i+1}))
        len=length(varargin{i+1});
        if(len==0)
            continue;
        elseif(len~=nrecs)
            error('Value vector for field %s not correct size: %d ~= %d',...
                varargin{i},len,nrecs);
        end
        varargin{i+1}=varargin{i+1}(:);
        if(group)
            varargin{i+1}=varargin{i+1}(:,ones(glen,1));
        end
    else
        dim=size(varargin{i+1});
        if(dim(1)~=nrecs || dim(2)<glen)
            error('Value array for field %s not correct size',varargin{i})
        end
    end
    
    % group field loop
    for j=1:glen
        % modify field if group
        if(group); f=[gf num2str(g(j))];
        else f=gf; end
        
        % assign values to field
        head=a2h(head,h,f,varargin{i+1}(:,j));
    end
end

% put head back in
head=mat2cell(head,h.size,ones(1,nrecs));
[data.head]=deal(head{:});

end

function [head]=a2h(head,h,f,v)
%A2H    Assign values to header field

% input based on field type (any new types will have to added
% to the header definition must be added here too)
for m=1:length(h.real)
    if(isfield(h.real(m).pos,f))
        head(h.real(m).pos.(f),:)=cell2mat(v);
        return;
    end
end
for m=1:length(h.int)
    if(isfield(h.int(m).pos,f))
        head(h.int(m).pos.(f),:)=round(cell2mat(v));
        return;
    end
end
for m=1:length(h.enum)
    if(isfield(h.enum(m).pos,f))
        % string enum values (lookup enum)
        if(iscellstr(v))
            % lookup only in enum(1)
            if(isfield(h.enum(1).val,v{1}) && length(unique(v))==1)
                head(h.enum(m).pos.(f),:)=h.enum(1).val.(v{1});
            elseif(all(isfield(h.enum(1).val,v)))
                % non-vectorized
                for n=1:length(v)
                    head(h.enum(m).pos.(f),n)=h.enum(1).val.(v{n});
                end
            else
                error('Some Enum IDs not valid for field %s',f)
            end
        % assume numeric cell array of enum values
        % (no check - breaks if mixed num/char cells)
        else head(h.enum(m).pos.(f),:)=round(cell2mat(v));
        end
        return;
    end
end
for m=1:length(h.lgc)
    if(isfield(h.lgc(m).pos,f))
        v=cell2mat(v);
        true=(v==h.true | v==116);
        false=(v==h.false | v==102);
        undef=(v==h.undef.ntype | v==117);
        if(~all(true | false | undef))
            error('Illogical values for logic field %s',f)
        end
        head(h.lgc(m).pos.(f),true)=h.true;
        head(h.lgc(m).pos.(f),false)=h.false;
        head(h.lgc(m).pos.(f),undef)=h.undef.ntype;
        return;
    end
end
for m=1:length(h.char)
    if(isfield(h.char(m).pos,f))
        n=h.char(m).pos.(f); 
        o=n(2)-n(1)+1;
        v=strnlen(v,o);
        head(n(1):n(2),:)=char(v).';
        return;
    end
end
            
% field not found
warning('SAClab:fieldInvalid','Invalid field: %s',f);

end
