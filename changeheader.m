function [data]=changeheader(data,varargin)
%CHANGEHEADER    Change SEIZMO data header values
%
%    Description: CHANGEHEADER(DATA,FIELD,VALUE) changes the header field 
%     FIELD to the value(s) in VALUE for each record in DATA and returns
%     an updated SEIZMO data structure.  FIELD must be a string 
%     corresponding to a valid header field.  VALUE may be a scalar 
%     (assigns same value to all) or a vector of length equal to the number
%     of records and can be of type numeric, logical, char, or cell.  If
%     FIELD is a group field ('t','resp','user','kt','kuser') then VALUE
%     can be a scalar (assigns same value to all fields for all records), a
%     vector (separate values for each record, but same value for each
%     field), or an array (separate values for all fields in all records).
%     For group field assignment, VALUE should be arranged so that columns
%     delimit values for each field in the group and rows delimit records 
%     (first row = first record, etc).  See the examples below to see how 
%     to replicate a set of values across several records.
%
%     CHANGEHEADER(DATA,FIELD1,VALUE1,...,FIELDN,VALUEN) allows changing
%     multiple fields in a single call.
%
%    Notes:
%     - Assigning a vector of values to a group field can behave 
%       inconsistently when working with multiple records.  Basically when
%       the number of records equals the number of fields, replication 
%       across records becomes replication across fields when there are 2 
%       or more header versions.  Avoid all this headache by not relying on
%       CHANGEHEADER to do the replication for vectors.  Expand the vector
%       to an array beforehand.  You have been warned!
%     - Passing a nan, inf, -inf value will change a numeric field to its
%       undefined value.  Use 'nan', 'undef' or 'undefined' to do the same
%       for a character field.  This is useful for not having to remember
%       what the field's actual undefined value is.
%
%    Tested on: Matlab r2007b
%
%    Header changes: Determined by input list.
%
%    Usage: data=changeheader(data,'field1',values1)
%           data=changeheader(data,'field1,'values1,...,'fieldN',valuesN)
%
%    Examples:
%     Some simple examples:
%      data=changeheader(data,'DELTA',getheader(data,'delta')*2);
%      data=changeheader(data,'STLA',lats,'STLO',lons)
%      data=changeheader(data,'KT0','sSKS');
%     
%     Copying the values of a group field from one record to other records
%     can have troubles with expansion if the number of records being
%     assigned to is the same size as the group field (ie 10 fields and 10
%     records) AND there are multiple data filetypes/versions (ie SAC
%     binary vs SEIZMO binary -- NOT timeseries vs spectral records). For
%     example, the following has behavior dependent on number of record
%     types in records 2-11.  If all records have the same filetype/version
%     then the values in record 1 copy to records 2-11 so their 't' fields
%     match.  If there are multiple filetypes though, this will cause each
%     record to get a separate value from record 1 and replicate that value
%     across it's 't' fields.
%      data=changeheader(data(2:11),'t',getheader(data(1),'t'))
%
%     This will always make the 't' fields of records 2-11 match record 1:
%      t=getheader(data(1),'t');
%      data=changeheader(data(2:11),'t',t(ones(10,1),:));
%     or
%      data=changeheader(data(2:11),'t',getheader(data(ones(10,1)),'t'));
%
%     Replication across fields will not have this problem.  Copy values in
%     field 't0' to rest of 't' fields for each record:
%         data=changeheader(data(1:10),'t',getheader(data(1:10),'t0'))
%
%    See also:  listheader, getheader, readheader, writeheader, getlgc,
%               getenumid, getenumdesc

%     Version History:
%        Oct. 29, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 28, 2008 - new sachp support
%        Feb. 18, 2008 - rewrite - parses by version before updating
%        Feb. 23, 2008 - major code cleanup
%        Feb. 28, 2008 - major code cleanup
%        Mar.  4, 2008 - cleanup errors and warnings
%        June 16, 2008 - doc update
%        June 20, 2008 - enum ids & logicals now support uppercase strings
%        June 28, 2008 - doc update
%        June 30, 2008 - undefine field supported using
%                        nan, inf, -inf, 'nan', 'undef', 'undefined'
%        Oct. 17, 2008 - added VINFO support, supports new struct layout
%        Nov. 16, 2008 - rename from CH to CHANGEHEADER, doc update, error
%                        msg update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 16, 2008 at 07:20 GMT

% todo:

% throw error if unpaired fields
if (mod(nargin-1,2))
    error('seizmo:changeheader:badNargs','Unpaired Field/Value!');
end

% check data structure
error(seizmocheck(data))

% quick exit
if(nargin==1); return; end

% number of records
nrecs=numel(data);

% recursive section
%   breaks up dataset with multiple calls so that the
%   rest of changeheader deals with only one header version
[h,idx]=versioninfo(data);
nver=numel(h);
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
                % dice up input if size is the same as the number of
                % records, otherwise replicate for each ch call. This
                % breaks group field replication when nrecs==ngfields.
                if(length(varargin{j})==nrecs)
                    temp{j}=varargin{j}(idx==i);
                else
                    temp{j}=varargin{j};
                end
            else
                % check and dice up input array
                dim=size(varargin{j},1);
                if(dim==nrecs)
                    temp{j}=varargin{j}(idx==i,:);
                else
                    error('seizmo:changeheader:invalidInputSize',...
                        'Value array for field %s incorrect size!',...
                        varargin{j-1});
                end
            end
        end
        data(idx==i)=changeheader(data(idx==i),temp{:});
    end
    return;
end

% pull entire header
head=[data.head];

% loop over field/value pairs
for i=1:2:(nargin-2)
    % force values into cell array
    if(isnumeric(varargin{i+1}) || islogical(varargin{i+1}))
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
        varargin{i+1}=varargin{i+1}(:);
        len=numel(varargin{i+1});
        if(len==0); continue; end
        % figure out group expansion
        if(group)
            % expand
            if(len==nrecs)
                varargin{i+1}=varargin{i+1}(:,ones(glen,1));
            elseif(len==glen)
                varargin{i+1}=varargin{i+1}(:,ones(nrecs,1)).';
            else
                error('seizmo:changeheader:invalidInputSize',...
                    ['\nFiletype: %s, Version: %d\n'...
                    'Group value vector for field %s incorrect size!'],...
                    h.filetype,h.version,varargin{i});
            end
        % check for correct length
        elseif(len~=nrecs)
            error('seizmo:changeheader:invalidInputSize',...
                ['\nFiletype: %s, Version: %d\n'...
                'Value vector for field %s incorrect size!'],...
                h.filetype,h.version,varargin{i});
        end
    % array
    else
        dim=size(varargin{i+1});
        if(dim(1)~=nrecs || dim(2)<glen)
            error('seizmo:changeheader:invalidInputSize',...
                ['\nFiletype: %s, Version: %d\n'...
                'Group value array for field %s not correct size'],...
                h.filetype,h.version,varargin{i})
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
%A2H    Subfunction that assigns values to header field

% special treatment for enums
for m=1:length(h.enum)
    if(isfield(h.enum(m).pos,f))
        % string id/desc enum values
        if(iscellstr(v))
            % lowercase values
            v=lower(v);
            % replace special undefined strings
            v(strcmpi(v,'undef'))={h.undef.stype};
            v(strcmpi(v,'undefined'))={h.undef.stype};
            v(strcmpi(v,'nan'))={h.undef.stype};
            % all are the same (quick)
            if(length(unique(v))==1)
                % all are ids
                if(isfield(h.enum(1).val,v{1}))
                    head(h.enum(m).pos.(f),:)=h.enum(1).val.(v{1});
                % all are undefined
                elseif(strcmpi(v{1},h.undef.stype))
                    head(h.enum(m).pos.(f),:)=h.undef.ntype;
                else
                    % check if description
                    desc=strcmpi(v{1},h.enum(1).desc);
                    if(any(desc))
                        head(h.enum(m).pos.(f),:)=find(desc)-1;
                    else
                        warning('seizmo:changeheader:enumBad',...
                            ['\nFiletype: %s, Version: %d\n'...
                            'Enum ID/Desc Invalid for field %s !'],...
                            h.filetype,h.version,f);
                    end
                end
            % not all are the same (slow b/c looped)
            else
                for n=1:length(v)
                    % check if string id
                    if(isfield(h.enum(1).val,v{n}))
                        head(h.enum(m).pos.(f),n)=h.enum(1).val.(v{n});
                    elseif(strcmpi(v{n},h.undef.stype))
                        head(h.enum(m).pos.(f),n)=h.undef.ntype;
                    else
                        % check if description
                        desc=strcmpi(v{n},h.enum(1).desc);
                        if(any(desc))
                            head(h.enum(m).pos.(f),n)=find(desc)-1;
                        else
                            warning('seizmo:changeheader:enumBad',...
                                ['\nFiletype: %s, Version: %d\n'...
                                'Enum ID/Desc Invalid for field %s !'],...
                                h.filetype,h.version,f);
                        end
                    end
                end
            end
        % assume numeric cell array of enum values
        % (no check - breaks if mixed num/char cells)
        else
            v=cell2mat(v);
            v(isnan(v))=h.undef.ntype;
            v(isinf(v))=h.undef.ntype;
            head(h.enum(m).pos.(f),:)=v;
        end
        return;
    end
end

% special treatment for logical words
for m=1:length(h.lgc)
    if(isfield(h.lgc(m).pos,f))
        if(iscellstr(v))
            % logic words (unknown => undefined here)
            true=strncmpi(v,'t',1);
            false=strncmpi(v,'f',1);
            undef=strncmpi(v,'u',1) | strcmpi(v,'nan');
            head(h.lgc(m).pos.(f),true)=h.true;
            head(h.lgc(m).pos.(f),false)=h.false;
            head(h.lgc(m).pos.(f),undef)=h.undef.ntype;
        else
            % logic numbers (unknown => unknown here)
            v=cell2mat(v);
            if(islogical(v))
                head(h.lgc(m).pos.(f),v)=h.true;
                head(h.lgc(m).pos.(f),~v)=h.false;
            else
                v(isnan(v))=h.undef.ntype;
                v(isinf(v))=h.undef.ntype;
                head(h.lgc(m).pos.(f),:)=v;
            end
        end
        return;
    end
end

% string types
for n=1:length(h.stype)
    for m=1:length(h.(h.stype{n}))
        if(isfield(h.(h.stype{n})(m).pos,f))
            p=h.(h.stype{n})(m).pos.(f);
            o=p(2)-p(1)+1;
            v(strcmpi(v,'undef'))={h.undef.stype};
            v(strcmpi(v,'undefined'))={h.undef.stype};
            v(strcmpi(v,'nan'))={h.undef.stype};
            v=strnlen(v,o);
            head(p(1):p(2),:)=char(v).';
            return;
        end
    end
end

% remainder of numeric types
for n=1:length(h.ntype)
    for m=1:length(h.(h.ntype{n}))
        if(isfield(h.(h.ntype{n})(m).pos,f))
            v=cell2mat(v);
            v(isnan(v))=h.undef.ntype;
            v(isinf(v))=h.undef.ntype;
            head(h.(h.ntype{n})(m).pos.(f),:)=v;
            return;
        end
    end
end
            
% field not found
warning('seizmo:changeheader:fieldInvalid',...
    'Filetype: %s, Version: %d\nInvalid field: %s !',...
    h.filetype,h.version,f);

end
