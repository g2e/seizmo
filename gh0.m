function [varargout]=gh(data,varargin)
%GH    Get SAClab header values
%
%    Description: Gets specified header values from SAClab structure.  
%     Fields must be strings corresponding to a valid header field or a 
%     valid group field (ie. t,kt,resp,user,kuser).  Values are returned 
%     in numeric arrays or cell string arrays oriented as a column vectors
%     with one value per cell.  One value array per field or group field.
%     Calling with no fields will attempt to return all headers in a
%     single numeric array (will crash for datasets with headers of
%     mixed sizes).
%    
%    Usage:    [value_array1,...]=gh(saclab_struct,'field1',...
%                                       'group_field1','field2')
%
%    Examples:
%     head=gh(data)       % put all header variables in one numeric array
%     times=gh(data,'t')  % put all t values in one array
%     dt=gh(data,'DeLtA')
%     [stla,stlo]=gh(data,'stla','STLO')
%
%    by Garrett Euler (2/2008)  ggeuler@wustl.edu
%
%    See also:  ch, lh, rh, wh, rpdw, rdata, rsac, bsac, wsac, sachp, gv

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

% number of files
nrecs=length(data);

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% warn if mixed versions
if(nver>1)
    warning('SAClab:mixedVersions','Dataset has varying header versions!');
end

% pull all header case (may fail if mixed versions!)
if(nargin==1); varargout{1}=[data.head]; end

% loop over given fields
for i=1:length(varargin)
    % initialize output
    varargout{i}=[];
    
    % force lowercase
    f=lower(varargin{i});
    
    % determine field properties in each header version
    valid=zeros(nver,1);
    string=zeros(nver,1);
    glen=ones(nver,1);
    g=cell(nver,1);
    group=zeros(nver,1);
    for j=1:nver
        % group?
        if(isfield(h(j).grp,f))
            g{j}=h(j).grp.(f).min:h(j).grp.(f).max;
            group(j)=1; glen(j)=length(g{j});
            ft=[f num2str(g{j}(1))]; % check if first is valid
        else
            ft=f;
        end
        
        % string?
        for k=1:length(h(j).stype)
            for m=1:length(h(j).(h(j).stype{k}))
                if(isfield(h(j).(h(j).stype{k})(m).pos,ft))
                    string(j)=1; valid(j)=1;
                end
            end
        end
        
        % numeric?
        if(~valid(j))
            for k=1:length(h(j).ntype)
                for m=1:length(h(j).(h(j).ntype{k}))
                    if(isfield(h(j).(h(j).ntype{k})(m).pos,ft))
                        valid(j)=1;
                    end
                end
            end
        end
        
        % invalid
        if(~valid(j))
            warning('SAClab:fieldInvalid','Invalid field: %s',ft);
        end
        clear ft;
    end
    
    % mixing check
    types=unique(string);
    valids=unique(valid);
    if(length(types)>1)
        error('mixed output types for a field not allowed: %s',f);
    end
    if(length(valids)>1)
        warning('SAClab:mixedValidity','Field has mixed validity: %s',f);
    end
    
    % if totally invalid move to next otherwise get column size
    if(length(valids)==1 && valids==0); continue;
    else ncols=max(glen);
    end
    
    % preallocate with zeros
    if(types); varargout{i}=repmat({'0'},nrecs,ncols);
    else varargout{i}=zeros(nrecs,ncols); % output is double precision
    end
    
    % loop over records
    for j=1:nrecs
        % logical index of header info
        v=data(j).version==vers;
        
        % replace zeros with undefined
        if(types); varargout{i}(j,:)=repmat({h(v).undef.stype},1,ncols);
        else varargout{i}(j,:)=ones(1,ncols)*h(v).undef.ntype;
        end
        
        % loop over group fields
        for k=1:glen(v)
            % group field name
            if(group(v)); gf=[f num2str(g{v}(k))];
            else gf=f;
            end
            
            % output by type
            if(types)
                for n=1:length(h(v).stype)
                    for m=1:length(h(v).(h(v).stype{n}))
                        if(isfield(h(v).(h(v).stype{n})(m).pos,gf))
                            p=h(v).(h(v).stype{n})(m).pos.(gf);
                            varargout{i}{j,k}=char(data(j).head(p(1):p(2)))';
                        end
                    end
                end
            else
                for n=1:length(h(v).ntype)
                    for m=1:length(h(v).(h(v).ntype{n}))
                        if(isfield(h(v).(h(v).ntype{n})(m).pos,gf))
                            varargout{i}(j,k)=data(j).head(h(v).(h(v).ntype{n})(m).pos.(gf));
                        end
                    end
                end
            end
        end
    end
end

end
