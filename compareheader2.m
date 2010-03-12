function []=compareheader2(data,varargin)
%COMPAREHEADER2    List SEIZMO headers side-by-side for easy comparison
%
%    Usage:    compareheader2(data)
%              compareheader2(data,'field1',...,'fieldN')
%
%    Description: COMPAREHEADER2(DATA) prints out a table comparing all of
%     the header fields of records in DATA.  Rows in the table correspond
%     to a specific record (see the first column for the record index and
%     the lookup table above the header value table for the record's path
%     & filename) while the columns correspond to different header fields.
%     This is basically like COMPAREHEADER but transposed.
%
%     COMPAREHEADER2(DATA,'FIELD1',...,'FIELDN') prints out the header
%     fields FIELD1 to FIELDN for records in DATA as a table.  FIELDS may
%     be normal fields ('b' 'kt1' 'xmaximum' etc), group fields ('t' 'kt'
%     etc), absolute fields ('t9 utc' 'user3 tai' 'resp utc' etc), or
%     wildcards ('*t1' '?' etc).  Only * and ? are valid wildcards.
%
%    Notes:
%     - group fields:    t, kt, user, kuser, resp, dep, st, ev, nz, nzdttm,
%                         kname, {real group} utc, {real group} tai
%     - list fields:     picks, all, full
%     - virtual fields:  nzmonth, nzcday, kzdttm, kzdate, kztime, z, ztai
%     - abs time fields: {real field} utc, {real field} tai
%
%    Examples:
%     Some simple cases:
%      compareheader2(data)          % compare all header variables
%      compareheader2(data,'t')      % compare t group
%
%     Fields are case independent:
%      compareheader2(data,'dEltA')
%      compareheader2(data,'StLA','stLo')
%
%     Compare picks:
%      compareheader2(data,'picks')
%
%    See also: COMPAREHEADER, LISTHEADER, GETHEADER, CHANGEHEADER

%     Version History:
%        Feb. 24, 2010 - initial version
%        Mar. 12, 2010 - fixes for Octave (explicit celling)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 12, 2010 at 01:40 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% headers setup
[h,idx]=versioninfo(data);
nh=numel(h);

% number of records
nrecs=numel(data);

% gather all possible fields and reftimes
fields=cell(2,5,nh); % just a guess
utcf=fields; taif=fields; utc6f=fields; tai6f=fields;
ref=nan(nrecs,5); good=false(nrecs,1);
vf=cell(nh,1); 
for i=1:nh
    % vf in wildcard search
    vf{i}=fieldnames(h(i).vf).';
    
    % all available field type sets
    for j=1:numel(h(i).types)
        for k=1:numel(h(i).(h(i).types{j}))
            fields(k,j,i)=...
                {fieldnames(h(i).(h(i).types{j})(k).pos).'};
            % special absolute time alternate fields for all real
            if(strcmp(h(i).types{j},'real'))
                utcf(k,j,i)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' utc')};
                taif(k,j,i)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' tai')};
                utc6f(k,j,i)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' 6utc')};
                tai6f(k,j,i)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' 6tai')};
            end
        end
    end
    
    % get reference times hack
    vidx=find(idx==i);
    [ref(vidx,:),good(vidx,1)]=vf_gh_z(h(i),[data(vidx).head]);
end

% only list each field once
% - this forces an alphabetical listing
fields=unique([[fields{:}].'; [vf{:}].']);
utcf=unique([utcf{:}].');
taif=unique([taif{:}].');
utc6f=unique([utc6f{:}].');
tai6f=unique([tai6f{:}].');

% list all if no fields given
if(nargin==1)
    varargin=fields;
end

% list file-lookup list
disp(' ')
disp(' RECORDS:')
disp('---------------------------')
for i=1:nrecs
    disp(sprintf('%d - %s',i,fullfile(data(i).path,data(i).name)))
end
disp('---------------------------')
disp(' ')

% loop over fields
v=cell(1,0); hf=cell(1,0); cw=nan(1,0); cnt=0;
nvarg=numel(varargin);
for i=1:nvarg
    % force lowercase
    gf=strtrim(lower(varargin{i}));
    wf=getwords(gf);
    
    % skip empty
    if(isempty(gf)); continue; end
    
    % check for vlist/vgrp
    group=false; glen=1; ngf=cell(nh,1);
    if(strcmp(gf,'all') || strcmp(gf,'full'))
        % list all
        ngf=fields; group=true; glen=numel(ngf);
    elseif(strcmp(gf,'picks'))
        % list picks - hardcoded...
        ngf={'kt0' 't0' 'kt1' 't1' 'kt2' 't2' 'kt3' 't3'...
            'kt4' 't4' 'kt5' 't5' 'kt6' 't6' 'kt7' 't7'...
            'kt8' 't8' 'kt9' 't9'};
        group=true; glen=numel(ngf);
    else
        % check for group fields (similar to list all case)
        for j=1:nh
            if(isfield(h(j).vgrp,wf{1}))
                ngf{j}=strtrim(strcat(h(j).vgrp.(wf{1}),{' '},...
                    joinwords(wf(2:end))));
                group=true;
            end
        end
        
        % clean up group list (only list each field once)
        % - this forces a sorted listing
        if(group)
            ngf=unique([ngf{:}]);
            glen=numel(ngf);
        end
    end
    
    % wildcard case (?==63,*==42) - pass to regexp
    wild=false;
    if(~group && (any(wf{1}==42) || any(wf{1}==63)))
        % declare as wildcard
        wild=true;
        % only show absolute time fields if explicitly sought
        if(strcmp(joinwords(wf(2:end)),'utc'))
            % no partial matches, trim excess blanks between field and utc
            ngf=['^' regexptranslate('wildcard',...
                strtrim(wf{1})) ' utc$'];
            % find matches in utc set
            ngf=utcf(~cellfun('isempty',regexp(utcf,ngf)));
        elseif(strcmp(joinwords(wf(2:end)),'tai'))
            % no partial matches, trim excess blanks between field and tai
            ngf=['^' regexptranslate('wildcard',...
                strtrim(wf{1})) ' tai$'];
            % find matches in tai set
            ngf=taif(~cellfun('isempty',regexp(taif,ngf)));
        elseif(strcmp(joinwords(wf(2:end)),'6utc'))
            % no partial matches, trim excess blanks between field and utc
            ngf=['^' regexptranslate('wildcard',...
                strtrim(wf{1})) ' 6utc$'];
            % find matches in utc set
            ngf=utc6f(~cellfun('isempty',regexp(utc6f,ngf)));
        elseif(strcmp(joinwords(wf(2:end)),'6tai'))
            % no partial matches, trim excess blanks between field and tai
            ngf=['^' regexptranslate('wildcard',...
                strtrim(wf{1})) ' 6tai$'];
            % find matches in tai set
            ngf=tai6f(~cellfun('isempty',regexp(tai6f,ngf)));
        elseif(numel(wf)==1)
            % no partial matches
            ngf=['^' regexptranslate('wildcard',gf) '$'];
            % find matches in normal fields
            ngf=fields(~cellfun('isempty',regexp(fields,ngf)));
        end
        glen=numel(ngf);
    end
    
    % loop over fields in group
    % - nongroup fields are treated as a group of 1
    for j=1:glen
        % modify field name if in a group
        if(group); f=ngf{j};
        elseif(wild); f=ngf{j};
        else f=gf;
        end
        
        % get minimum column width required by field name
        cn=numel(f)+2; cs=num2str(cn);
        
        % get value formatted as string (right justified)
        values=cell(nrecs,1);
        for k=1:nrecs
            values{k}=cmph_disp(h(idx(k)),f,data(k),cs,ref(k,:),good(k));
        end
        
        % justify
        values=strtrim(strjust(char(values),'right'));
        
        % now get minimum column width
        cn=max(cn,size(values,2)+2); cs=num2str(cn);
        
        % fill new cell
        cnt=cnt+1;
        cw(cnt)=cn;
        hf{cnt}=sprintf(['%' cs 's'],upper(f));
        v{cnt}=char(strcat({32*ones(1,cn-size(values,2))},values));
    end
end

% table header
disp( sprintf('%14s','      \ HEADER'))
disp( sprintf('%14s','RECORD \ FIELD'))
disp([sprintf('%11s',' NUMBER \  ') hf{:}])
disp(char(45*ones(1,11+sum(cw))))

% make record number column
rnc=num2str((1:nrecs)');
rnc=strcat({char(32*ones(1,11-size(rnc,2)-3))},rnc,{' | '});

% display table
disp([char(rnc),v{:}]);

end

function [string]=cmph_disp(h,f,data,cs,reftime,good)
%CMPH_DISP    Returns the value of a header field as a string

% virtual fields
if(isfield(h.vf,f))
    switch h.vf.(f).type
        case 'enum'
            % searches only enum(1) for now (prefer one big
            % set of enums to share vs separate sets)
            ival=h.vf.(f).gh(h,data.head);
            if(ival==fix(ival) && ival>=h.enum(1).minval && ...
                    ival<=h.enum(1).maxval)
                string=sprintf(['%' cs 's'],h.enum(1).id{ival+1});
            elseif(ival==h.undef.ntype)
                string=sprintf(['%' cs 's'],...
                    sprintf('UNDEFINED (%g)',ival));
            else
                string=sprintf(['%' cs 's'],sprintf('UNKNOWN (%g)',ival));
            end
        case 'lgc'
            ival=h.vf.(f).gh(h,data.head);
            switch ival
                case h.false
                    string=sprintf(['%' cs 's'],'FALSE');
                case h.true
                    string=sprintf(['%' cs 's'],'TRUE');
                case h.undef.ntype
                    string=sprintf(['%' cs 's'],...
                        sprintf('UNDEFINED (%g)',ival));
                otherwise
                    string=sprintf(['%' cs 's'],...
                        sprintf('INVALID (%g)',ival));
            end
        case 'char'
            ival=h.vf.(f).gh(h,data.head);
            string=sprintf(['%' cs 's'],ival{:});
        case 'abs'
            ival=h.vf.(f).lh(h,data.head);
            string=sprintf(['%' cs 's'],ival{:});
        otherwise
            ival=h.vf.(f).gh(h,data.head);
            if(ival==h.undef.ntype)
                string=sprintf(['%' cs 's'],sprintf('UNDEFINED (%g)',...
                    ival));
            else
                string=sprintf(['%' cs 's'],sprintf('%-.10g',...
                    ival));
            end
    end
    return;
end

% handle enum/lgc fields
for m=1:numel(h.enum)
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
for m=1:numel(h.lgc)
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
for m=1:numel(h.stype)
    for n=1:numel(h.(h.stype{m}))
        if(isfield(h.(h.stype{m})(n).pos,f))
            p=h.(h.stype{m})(n).pos.(f);
            q=p(2)-p(1)+1; u=numel(h.undef.stype);
            if(strcmp([h.undef.stype ones(1,q-u)*32],...
                    char(data.head(p(1):p(2)).')))
                string=sprintf(['%' cs 's'],sprintf('UNDEFINED (%s)',...
                    char(data.head(p(1):p(2)))));
            else
                string=sprintf(['%' cs 's'],char(data.head(p(1):p(2))));
            end
            return;
        end
    end
end

% remaining numeric types
for m=1:numel(h.ntype)
    for n=1:numel(h.(h.ntype{m}))
        if(isfield(h.(h.ntype{m})(n).pos,f))
            if(data.head(h.(h.ntype{m})(n).pos.(f))==h.undef.ntype)
                string=sprintf(['%' cs 's'],sprintf('UNDEFINED (%g)',...
                    data.head(h.(h.ntype{m})(n).pos.(f))));
            else
                string=sprintf(['%' cs 's'],sprintf('%-.10g',...
                    data.head(h.(h.ntype{m})(n).pos.(f))));
            end
            return;
        elseif(strcmpi(h.ntype{n},'real'))
            % absolute time fields section
            wf=getwords(f);
            if(isfield(h.real(n).pos,wf{1}))
                v=data.head(h.real(n).pos.(wf{1}));
                if(any(strcmpi(joinwords(wf(2:end)),{'utc' '6utc'})))
                    if(v==h.undef.ntype || isnan(v) || isinf(v))
                        string=sprintf(['%' cs 's'],...
                            sprintf('UNDEFINED (%g)',v));
                    elseif(~good)
                        string=sprintf(['%' cs 's'],...
                            sprintf('NO REFTIME (%g)',v));
                    else
                        % get values for output
                        utc=tai2utc(reftime+[0 0 0 0 v]);
                        cal=doy2cal(utc(1:2));
                        utc(5)=round(1000*utc(5));
                        string=sprintf(['%' cs 's'],...
                            sprintf(['%04d-%02d-%02d (%03d) '...
                            '%02d:%02d:%02d.%03d'],utc(1),cal(2),cal(3),...
                            utc(2),utc(3),utc(4),fix(utc(5)/1000),...
                            mod(utc(5),1000)));
                    end
                    return;
                elseif(any(strcmpi(joinwords(wf(2:end)),{'tai' '6tai'})))
                    if(v==h.undef.ntype || isnan(v) || isinf(v))
                        string=sprintf(['%' cs 's'],...
                            sprintf('UNDEFINED (%g)',v));
                    elseif(~good)
                        string=sprintf(['%' cs 's'],...
                            sprintf('NO REFTIME (%g)',v));
                    else
                        tai=fixtimes(reftime+[0 0 0 0 v]);
                        cal=doy2cal(tai(1:2));
                        tai(5)=round(1000*tai(5));
                        string=sprintf(['%' cs 's'],...
                            sprintf(['%04d-%02d-%02d (%03d) '...
                            '%02d:%02d:%02d.%03d'],tai(1),cal(2),cal(3),...
                            tai(2),tai(3),tai(4),fix(tai(5)/1000),...
                            mod(tai(5),1000)));
                    end
                    return;
                end
            end
        end
    end
end

% field not found
string=sprintf(['%' cs 's'],'NOT A FIELD!');

end
