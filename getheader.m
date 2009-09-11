function [varargout]=getheader(data,varargin)
%GETHEADER    Get SEIZMO data header values
%
%    Usage:  headers=getheader(data)
%            values=getheader(data,'field1')
%            [values1,...,valuesN]=getheader(data,'field1',...,'fieldN')
%
%    Description: GETHEADER(DATA) will return header values for all records
%     in DATA as a single numeric array.  Rows in the output array
%     correspond to the header values of individual records.  The order of
%     the fields follows that of how they are stored in the SEIZMO struct
%     (see SEIZMO's function SEISMODEF for details).  Character fields are 
%     returned as a series of their ascii number equivalents (0-255).
%
%     GETHEADER(DATA,FIELD) returns the specified header field FIELD's
%     values for each record stored in the SEIZMO data structure DATA.
%     FIELD must be a string corresponding to a valid header field or a
%     valid group field (ie. t,kt,resp,user,kuser).  Values are returned in
%     numeric arrays or cellstring arrays oriented such that each column
%     corresponds to an individual header field and each row to an
%     individual record.  So the group field 't' would return a numeric
%     array with 10 columns and as many rows as records in DATA while group
%     field 'kuser' would return a cellstring array with 3 columns and as
%     many rows as records in DATA.
%     
%     GETHEADER(DATA,FIELD1,...,FIELDN) returns one array of values per
%     field or group field.
%
%    Notes:
%     - Enumerated fields return the value actually stored, an integer used
%       to look up the enum string id and description in a table.  To 
%       retrieve the associated string id or description use the functions 
%       GETENUM or GETENUMDESC.
%     - Logical fields return the value actually stored, not a logical.  To
%       get a more useful value use GETLGC.
%
%    Examples:
%     Put all t series values in one array:
%      times=getheader(data,'t')
%
%     Pull just the sample rates for records (fields are case insensitive):
%      dt=getheader(data,'DeLtA')
%
%     Get the station lat and lon for records:
%      [stla,stlo]=getheader(data,'stla','STLO')
%
%     Enumerated fields return the table lookup index which is
%     the value stored in the header:
%      getheader(data,'iftype')
%
%    See also:  changeheader, listheader, readheader, writeheader, getlgc,
%               getenumid, getenumdesc, getncmp, getarrival, compareheader

%     Version History:
%        Oct. 29, 2008 - initial version
%        Nov.  7, 2008 - doc update
%        Jan. 28, 2008 - new sachp support
%        Feb. 18, 2008 - rewrite - parses by version first
%        Feb. 23, 2008 - minor doc update
%        Feb. 28, 2008 - code cleanup
%        Feb. 29, 2008 - minor code cleaning
%        Mar.  4, 2008 - doc update
%        June 12, 2008 - doc update, full header dump fixes
%        Oct. 17, 2008 - added VINFO support, supports new struct layout
%        Nov. 16, 2008 - update for new name schema (now GETHEADER)
%        Apr. 23, 2009 - fix seizmocheck for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:05 GMT

% todo:

% require at least one input
if(nargin<1)
    error('seizmo:getheader:notEnoughInputs',...
        'Not enough input arguments.')
end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% number of records
nrecs=numel(data);

% recursive section
%   breaks up dataset with multiple calls so that the
%   rest of gh deals with only one header version
[h,idx]=versioninfo(data);
nver=numel(h);
sfill=repmat({'NaN'},nrecs,1);
nfill=nan(nrecs,1);
if(nver>1)
    nout=max([1 nargin-1]);
    for i=1:nver
        varargoutn=cell(1,nout);
        [varargoutn{:}]=getheader(data(idx==i),varargin{:});
        % assign to varargout
        for j=1:nout
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
            varargout{j}(idx==i,1:in)=varargoutn{j};
        end
    end
    return;
end

% pull entire header
head=[data.head];

% push out entire header
if(nargin==1); varargout{1}=head.'; return; end

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
warning('seizmo:getheader:fieldInvalid',...
    'Filetype: %s, Version: %d\nInvalid field: %s',h.filetype,h.version,f);
head=nan;
type=0;

end
