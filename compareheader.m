function []=compareheader(data,varargin)
%COMPAREHEADER    List SEIZMO data headers side-by-side for easy comparison
%
%    Usage:    compareheader(data)
%              compareheader(data,'field1',...,'fieldN')
%
%    Description: COMPAREHEADER(DATA) prints out a table comparing all of
%     the header fields of records in DATA.  Rows in the table correspond
%     to a specific header field (see the first column for the field name)
%     and the columns correspond to different records (there is a record
%     path and filename listing above the table as an aid).
%
%     COMPAREHEADER(DATA,'FIELD1',...,'FIELDN') prints out the header
%     fields FIELD1 to FIELDN for records in DATA as a table.  FIELDS may
%     be normal fields ('b' 'kt1' 'xmaximum' etc), group fields ('t' 'kt'
%     etc), or wildcards ('*t1' '?' etc).  Only * and ? are valid
%     wildcards.
%
%    Notes:
%
%    Examples:
%     Some simple cases:
%      compareheader(data)          % compare all header variables
%      compareheader(data,'t')      % compare t group
%
%     Fields are case independent:
%      compareheader(data,'dEltA')
%      compareheader(data,'StLA','stLo')
%
%    See also: listheader, getheader, changeheader

%     Version History:
%        Sep. 11, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 11, 2009 at 03:40 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% headers setup
[h,idx]=versioninfo(data);

% number of records
nrecs=numel(data);

% define column width
cn=25; cs=num2str(cn);

% gather all possible fields
fields=cell(1,1,1);
for i=1:numel(h)
    % all available field type sets
    for j=1:length(h(i).types)
        for k=1:length(h(i).(h(i).types{j}))
            fields{j,k,i}=...
                fieldnames(h(i).(h(i).types{j})(k).pos)';
        end
    end
end
% only list each field once
% - this forces an alphabetical listing
fields=unique([fields{:}].');

% list all if no fields given
if(nargin==1)
    varargin=fields;
end

% list file-lookup list
disp(' ')
disp(' RECORDS:')
disp('---------------------------')
for i=1:nrecs
    disp(sprintf('%d - %s%s',i,data(i).path,data(i).name))
end
disp('---------------------------')
disp(' ')

% table header
disp(sprintf('%15s','\    '))
disp([sprintf('%15s','FIELD    \   ') '   RECORD NUMBER'])
disp([sprintf('%15s','\  ') sprintf(['%' cs 'd'],1:nrecs)])
disp(char(45*ones(1,15+cn*nrecs)))

% loop over fields
nvarg=numel(varargin);
for i=1:nvarg
    % force lowercase
    gf=lower(varargin{i});
    
    % check for group fields (similar to list all case)
    group=false; glen=1; g=cell(nrecs,1);
    for j=1:nrecs
        if(isfield(h(idx(j)).grp,gf))
            g{j}=h(idx(j)).grp.(gf).min:h(idx(j)).grp.(gf).max;
            group=true;
        end
    end
    
    % clean up group list (only list each field once)
    % - this forces a sorted listing
    if(group)
        g=unique([g{:}]);
        glen=numel(g);
    end
    
    % wildcard case (?==63,*==42) - pass to regexpi
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
        for j=1:numel(tmp)
            ngf=[ngf(1:tmp(j)-2+j) '.' ngf(tmp(j)+j-1:end)];
        end
        % now find matches in fields
        ngf=fields(~cellfun('isempty',regexp(fields,ngf)));
        glen=numel(ngf);
    end
    
    % loop over fields in group
    % - nongroup fields are treated as a group of 1
    for j=1:glen
        % modify field name if in a group
        if(group); f=[gf num2str(g(j))];
        elseif(wild); f=ngf{j};
        else f=gf;
        end
        
        % get value formatted as string (20 chars, right justified)
        values=nan(nrecs,cn);
        for k=1:nrecs
            values(k,:)=cmph_disp(h(idx(k)),f,data(k),cs,cn);
        end
        
        % display
        values=values.';
        disp([sprintf('%12s | ',upper(f)) char(values(:).')])
    end
end

end

function [string]=cmph_disp(h,f,data,cs,cn)
%CMPH_DISP    Returns the value of a header field as a string

% handle enum/lgc fields specially
for m=1:length(h.enum)
    if(isfield(h.enum(m).pos,f))
        ival=data.head(h.enum(m).pos.(f));
        % searches only enum(1) for now (prefer one big
        % set of enums to share vs separate sets)
        if(ival==fix(ival) && ival>=h.enum(1).minval && ...
                ival<=h.enum(1).maxval)
            string=sprintf(['%' cs 's'],h.enum(1).id{ival+1});
        elseif(ival==h.undef.ntype)
            string=sprintf(['%' cs 's'],sprintf('UNDEFINED (%g)',ival));
        else
            string=sprintf(['%' cs 's'],sprintf('UNKNOWN (%g)',ival));
        end
        return;
    end
end
for m=1:length(h.lgc)
    if(isfield(h.lgc(m).pos,f))
        if(data.head(h.lgc(m).pos.(f))==h.false)
            string=sprintf(['%' cs 's'],'FALSE');
        elseif(data.head(h.lgc(m).pos.(f))==h.true)
            string=sprintf(['%' cs 's'],'TRUE');  
        elseif(data.head(h.lgc(m).pos.(f))==h.undef.ntype)
            string=sprintf(['%' cs 's'],sprintf('UNDEFINED (%g)',...
                data.head(h.lgc(m).pos.(f))));
        else
            string=sprintf(['%' cs 's'],sprintf('INVALID (%g)',...
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
                % give 2 char padding
                p(2)=p(1)+min(cn-15,p(2)-p(1));
                string=sprintf(['%' cs 's'],sprintf('UNDEFINED (%s)',...
                    char(data.head(p(1):p(2)))));
            else
                % give 2 char padding
                p(2)=p(1)+min(cn-3,p(2)-p(1));
                string=sprintf(['%' cs 's'],char(data.head(p(1):p(2))));
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
                string=sprintf(['%' cs 's'],sprintf('UNDEFINED (%g)',...
                    data.head(h.(h.ntype{m})(n).pos.(f))));
            else
                string=sprintf(['%' cs 's'],sprintf('%-.10g',...
                    data.head(h.(h.ntype{m})(n).pos.(f))));
            end
            return;
        end
    end
end

% field not found
string=sprintf(['%' cs 's'],'NOT A FIELD!');

end

