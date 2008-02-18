function [data]=ch(data,varargin)
%CH    Change SAClab headers
%
%    Description: Changes the specified header field(s) to the specified 
%     value(s) for data in a SAClab structure.  The field variable must be
%     a string corresponding to a valid header field.  Values may be scalar
%     (assign same value to all) or a vector of length equal to the number
%     of records.  Values may be contained in a numeric, char, or cell
%     array.  Char arrays must be arranged as a column vector with one
%     string value per row.  Cell arrays must have one value per cell.
%
%    Usage:    new_struct=ch(saclab_struct,'group_field1',value1,...
%                               'field2',{'value2a' 'value2b' ...},...
%                               'field3',[value3a; value3b; ...])
%
%    Examples:
%     data=ch(data,'DELTA',dt);
%     data=ch(data,'STLA',lats,'STLO',lons)
%     data=ch(data,'KT0','sSKS');
%
%    by Garrett Euler (2/2008)  ggeuler@wustl.edu
%
%    See also:  lh, gh, rh, wh, rpdw, rdata, rsac, bsac, wsac, sachp, gv

% do nothing on no input
if(nargin<=1); return; end

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

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% check header field types match up to what ch currently supports
chtypes={'real' 'int' 'enum' 'lgc' 'char'};
for i=1:nver
    for j=1:length(h(i).types)
        if(~any(strcmp(h(i).types{j},chtypes)))
            warning('SAClab:missingTypes',...
                'Header definition has unsupported field type, %s',...
                h(i).types{j});
        end
    end
end

% warn if mixed versions
if(nver>1)
    warning('SAClab:mixedVersions','Dataset has varying header versions!');
end

% loop over field/value pairs
for i=1:2:(nargin-1)
    % force values into cell array
    if(isnumeric(varargin{i+1}))
        varargin{i+1}=num2cell(varargin{i+1});
    elseif(ischar(varargin{i+1}))
        varargin{i+1}=cellstr(varargin{i+1});
    end
    
    % input checks - values in vector/scalar form
    if(isvector(varargin{i+1}))
        % get vector length
        len=length(varargin{i+1});
        
        % act based on length
        if(len==0)
            % no values - skip field
            warning('SAClab:noValues','Empty value array for field %s',...
                varargin{i});
            continue;
        elseif(len==1)
            % 1 value - replicate
            varargin{i+1}=repmat(varargin{i+1},nrecs,1);
        elseif(len~=nrecs)
            % bad size - skip field
            warning('SAClab:valuesBadSize',...
                'Value vector for field %s not correct size: %d ~= %d',...
                varargin{i},len,nrecs);
            continue;
        end
        
        % force values to be in column vector
        varargin{i+1}=varargin{i+1}(:);
    % input checks - values in an array
    else
        % get array dimensions
        dim=size(varargin{i+1});
        
        % check if a dimension matches number of records
        match=find(dim==nrecs,1);
        if(min(dim)==0)
            % empty array (0x0) - skip field
            warning('SAClab:noValues','Empty value array for field %s',...
                varargin{i});
            continue;
        elseif(isempty(match))
            % bad size - skip field
            error('value array not correct size for field %s',varargin{i})
        elseif(match==2)
            % array on its side - transpose
            varargin{i+1}=varargin{i+1}.';
        end
    end
    
    % force field name to be lowercase
    gf=lower(varargin{i});
    
    % loop over records
    for j=1:nrecs
        % logical index of header info
        v=data(j).version==vers;
        
        % check for group fields
        group=0; glen=1;
        if(isfield(h(v).grp,gf))
            g=h(v).grp.(gf).min:h(v).grp.(gf).max;
            group=1; glen=length(g);
        end
        
        % checking group field size vs value array size
        if(size(varargin{i+1},2)<glen)
            % not enough columns - skip
            warning('SAClab:notEnoughValues',...
                'Not enough values given for group field %s, record %d',gf,j);
            continue;
        end
        
        % group field loop
        for k=1:glen
            % modify field if group
            if(group)
                f=[gf num2str(g(k))];
                
                % if values are vector for group field (sets group fields
                % to same value), reindex rather than replicate to save
                % memory and time
                if(isvector(varargin{i+1})); k=1; end %#ok<FXSET>
            else
                f=gf;
            end
            
            % valid field flag
            set=0;
            
            % input based on field type (any new types will have to added
            % to the header definition must be added here too)
            for m=1:length(h(v).real)
                if(isfield(h(v).real(m).pos,f))
                    data(j).head(h(v).real(m).pos.(f))=varargin{i+1}{j,k};
                    set=1;
                end
            end
            for m=1:length(h(v).int)
                if(isfield(h(v).int(m).pos,f))
                    data(j).head(h(v).int(m).pos.(f))=round(varargin{i+1}{j,k});
                    set=1;
                end
            end
            for m=1:length(h(v).enum)
                if(isfield(h(v).enum(m).pos,f))
                    % string enum value (lookup enum)
                    if(ischar(varargin{i+1}{j,k}))
                        % searches only enum(1) for now (prefer one big
                        % set of enums to share vs separate sets)
                        if(isfield(h(v).enum(1).val,varargin{i+1}{j,k}))
                            data(j).head(h(v).enum(m).pos.(f))=...
                                h(v).enum(1).val.(varargin{i+1}{j,k});
                        else
                            warning('SAClab:enumIdInvalid',...
                                'Enum ID %s not valid - no change to header made.',...
                                varargin{i+1}{j,k});
                        end
                    % numeric enum value (no check)
                    elseif(isnumeric(varargin{i+1}{j,k}))
                        data(j).head(h(v).enum(m).pos.(f))=...
                            round(varargin{i+1}{j,k});
                    end
                    set=1;
                end
            end
            for m=1:length(h(v).lgc)
                if(isfield(h(v).lgc(m).pos,f))
                    % false
                    if(varargin{i+1}{j,k}(1)==h(v).false || ...
                            varargin{i+1}{j,k}(1)==102) % always allow string 'f...'
                        data(j).head(h(v).lgc(m).pos.(f))=h(v).false;
                    % true
                    elseif(varargin{i+1}{j,k}(1)==h(v).true || ...
                            varargin{i+1}{j,k}(1)==116) % always allow string 't...'
                        data(j).head(h(v).lgc(m).pos.(f))=h(v).true;
                    % illogical
                    else
                        if(ischar(varargin{i+1}{j,k}))
                            warning('SAClab:logicInvalid',...
                                '%s illogical - no change to header made.',...
                                varargin{i+1}{j,k});
                        else
                            warning('SAClab:logicInvalid',...
                                '%d illogical - no change to header made.',...
                                varargin{i+1}{j,k});
                        end
                    end
                    set=1;
                end
            end
            for m=1:length(h(v).char)
                if(isfield(h(v).char(m).pos,f))
                    % get string lengths
                    n=h(v).char(m).pos.(f); 
                    o=n(2)-n(1)+1; p=length(varargin{i+1}{j,k});
                    if (p>o)
                        % truncate string
                        data(j).head(n(1):n(2))=varargin{i+1}{j,k}(1:o);
                    else
                        % pad string with spaces
                        data(j).head(n(1):n(2))=...
                            [varargin{i+1}{j,k}(:).' ones(1,o-p)*32];
                    end
                    set=1;
                end
            end
            
            % check field was found
            if(~set)
                warning('SAClab:fieldInvalid','Invalid field: %s',field);
            end
        end
    end
end

end

