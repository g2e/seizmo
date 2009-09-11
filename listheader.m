function []=listheader(data,varargin)
%LISTHEADER    List SEIZMO data headers
%
%    Usage:    listheader(data)
%              listheader(data,'field1','field2',...)
%
%    Description: LISTHEADER(DATA) lists the entire header of all records
%     in DATA with a format similar to SAC's lh command.
%
%     LISTHEADER(DATA,FIELD1,...,FIELDN) lists the header field(s) FIELD1
%     to FIELDN and their value(s) from records in DATA in a manner similar
%     to SAC's lh command.  FIELDS must be a string corresponding to a
%     valid header field, group field (ie. t, kt, resp, user, kuser), or
%     wildcards ('nz*' 'dep*' etc).  Only * and ? are valid wildcards.
%
%    Notes:
%
%    Examples:
%      listheader(data)          % lists all header variables
%      listheader(data,'t')      % lists t group
%
%     Fields are case independent:
%      listheader(data,'dEltA')
%      listheader(data,'StLA','stLo')
%
%    See also:  compareheader, changeheader, getheader

%     Version History:
%        Oct. 29, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 28, 2008 - new sachp support
%        Feb. 23, 2008 - minor doc update
%        Feb. 29, 2008 - major code overhaul
%        Mar.  4, 2008 - minor doc update
%        June 12, 2008 - doc update, added history
%        Oct. 17, 2008 - added VINFO support, supports new struct layout
%        Oct. 27, 2008 - update for struct change, doc update, history fix
%        Nov. 16, 2008 - update for new name schema (now listheader)
%        Apr. 23, 2009 - move usage up
%        Aug.  4, 2009 - changed numeric format from %d to %g so Octave
%                        output is not integer
%        Sep.  4, 2009 - changed numeric format again from %g to %-.10g
%        Sep. 11, 2009 - add support for wildcards * and ?, skip on empty
%                        fix, better nargin checking
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 11, 2009 at 06:15 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% headers setup
[h,idx]=versioninfo(data);

% extra padding
disp(' ')

% loop over files
for i=1:numel(data)
    % all available fields
    fields=cell(1,1);
    for j=1:length(h(idx(i)).types)
        for k=1:length(h(idx(i)).(h(idx(i)).types{j}))
            fields{k,j}=...
                fieldnames(h(idx(i)).(h(idx(i)).types{j})(k).pos).';
        end
    end
    fields=[fields{:}].';
    
    % list all case
    if (nargin==1)
        varargin=fields;
    end
    
    % get filename
    name=fullfile(data(i).path,data(i).name);
    
    % formatted header
    disp(' ')
    disp(sprintf(' FILE: %s - %d',name,i))
    disp('---------------------------')
    
    % loop over fields
    for j=1:numel(varargin)
        % force lowercase
        gf=lower(varargin{j});
        
        % skip empty
        if(isempty(gf)); continue; end
        
        % check for group fields
        group=0; glen=1;
        if(isfield(h(idx(i)).grp,gf))
            g=h(idx(i)).grp.(gf).min:h(idx(i)).grp.(gf).max;
            group=1; glen=length(g);
        end
        
        % wildcard case (?==63,*==42) - pass to regexp
        wild=false;
        if(~group && (any(gf==42) || any(gf==63)))
            % declare as wildcard
            wild=true;
            % no partial matches
            ngf=['^' gf '$'];
            % replace ? with .
            ngf(ngf==63)=46;
            % replace * with .*
            tmp=find(ngf==42);
            for k=1:numel(tmp)
                ngf=[ngf(1:tmp(k)-2+k) '.' ngf(tmp(k)+k-1:end)];
            end
            % now find matches in fields
            ngf=fields(~cellfun('isempty',regexp(fields,ngf)));
            glen=numel(ngf);
        end
        
        % group loop
        for k=1:glen
            % modify field if in a group
            if(group); f=[gf num2str(g(k))];
            elseif(wild); f=ngf{k};
            else f=gf;
            end
            
            % display field/value
            lh_disp(h(idx(i)),f,data(i));
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
            disp(sprintf('%12s = UNDEFINED (%g)',upper(f),ival));
        else
            disp(sprintf('%12s = UNKNOWN (%g)',upper(f),ival));
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
            disp(sprintf('%12s = UNDEFINED (%g)',upper(f),...
                data.head(h.lgc(m).pos.(f))));
        else
            disp(sprintf('%12s = INVALID (%g)',upper(f),...
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
                disp(sprintf('%12s = UNDEFINED (%g)',upper(f),...
                    data.head(h.(h.ntype{m})(n).pos.(f))));
            else
                disp(sprintf('%12s = %-.10g',upper(f),...
                    data.head(h.(h.ntype{m})(n).pos.(f))));
            end
            return;
        end
    end
end

% field not found
warning('seizmo:listheader:fieldInvalid',...
    'Filetype: %s, Version: %d\nInvalid field: %s',h.filetype,h.version,f);

end
