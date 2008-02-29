function [varargout]=gh(data,varargin)
%GH    Get seislab header values
%
%    Description: Gets specified header values from seislab data structure.  
%     Fields must be strings corresponding to a valid header field or a 
%     valid group field (ie. t,kt,resp,user,kuser).  Values are returned 
%     in numeric arrays or cell string arrays oriented as a column vectors
%     with one value per cell.  One value array per field or group field.
%     Calling with no fields will attempt to return all headers in a
%     single numeric array.
%    
%    Usage:    [value_array1,...]=gh(seislab_struct,'field1',...
%                                       'group_field1','field2',...)
%
%    Examples:
%     head=gh(data)       % put all header variables in one numeric array
%     times=gh(data,'t')  % put all t values in one array
%     dt=gh(data,'DeLtA')
%     [stla,stlo]=gh(data,'stla','STLO')
%
%    See also:  ch, lh, rh, wh, glgc, genum, genumdesc

% require at least one input
if(nargin<1)
    error('MATLAB:nargchk:notEnoughInputs',...
        'Not enough input arguments.')
end

% check data structure
error(seischk(data))

% number of records
nrecs=length(data);

% recursive section
%   breaks up dataset with multiple calls so that the
%   rest of gh deals with only one header version
v=[data.version];
vers=unique(v(:));
nver=length(vers);
sfill=repmat({'NaN'},nrecs,1);
nfill=nan(nrecs,1);
if(nver>1)
    for i=1:nver
        varargoutn=cell(1,nargin-1);
        [varargoutn{:}]=gh(data(vers(i)==v),varargin{:});
        % assign to varargout
        for j=1:nargin-1
            % preallocate by type
            if(i==1)
                if(iscellstr(varargoutn{j}))
                    varargout{j}=sfill; type=1;
                else
                    varargout{j}=nfill; type=0;
                end
            end
            % expand
            in=size(varargoutn{j},2);
            out=size(varargout{j},2);
            if(in>out)
                if(type)
                    varargout{j}(:,out+1:in)=sfill(:,ones(1,in-out));
                else
                    varargout{j}(:,out+1:in)=nfill(:,ones(1,in-out));
                end
            end
            % assign
            varargout{j}(vers(i)==v,1:in)=varargoutn{j};
        end
    end
    return;
end

% headers setup
h=seishi(vers);

% pull entire header
head=[data.head];

% push out entire header
if(nargin==1); varargout{1}=head; return; end

% loop over fields
for i=1:nargin-1
    % force field to be lowercase
    gf=lower(varargin{i});
    
    % check for group fields
    group=0; glen=1;
    if(isfield(h.grp,gf))
        g=h.grp.(gf).min:h.grp.(gf).max;
        group=1; glen=length(g);
    end
    
    % group field loop
    for j=1:glen
        % modify field if group
        if(group); f=[gf num2str(g(j))];
        else f=gf; end
        
        % pull header values
        [val,type]=ph(head,h,f);
        
        % preallocate (did not know the type until now)
        if(j==1)
            if(type)
                varargout{i}=repmat({'NaN'},nrecs,glen);
            else
                varargout{i}=nan(nrecs,glen);
            end
        end
        
        % assign to varargout
        varargout{i}(:,j)=val;
    end
end

end

function [head,type]=ph(head,h,f)
%PH    Pull header values

% output by type
for n=1:length(h.stype)
    for m=1:length(h.(h.stype{n}))
        if(isfield(h.(h.stype{n})(m).pos,f))
            p=h.(h.stype{n})(m).pos.(f);
            head=cellstr(char(head(p(1):p(2),:).'));
            type=1; return;
        end
    end
end
for n=1:length(h.ntype)
    for m=1:length(h.(h.ntype{n}))
        if(isfield(h.(h.ntype{n})(m).pos,f))
            head=head(h.(h.ntype{n})(m).pos.(f),:).';
            type=0; return;
        end
    end
end

% field not found
warning('seislab:gh:fieldInvalid',...
    '\nHeader Version: %d\nInvalid field: %s',h.version,f);
head=nan;
type=0;

end

