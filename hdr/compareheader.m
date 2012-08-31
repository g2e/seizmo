function []=compareheader(data,varargin)
%COMPAREHEADER    List SEIZMO headers in table form for easy comparison
%
%    Usage:    compareheader(data)
%              compareheader(data,'field1',...,'fieldN')
%
%    Description:
%     COMPAREHEADER(DATA) prints out a table comparing all of the header
%     fields of records in DATA.  Rows in the table correspond to a
%     specific header field (see the first column for the field names) and
%     the columns correspond to different records (the record index is
%     given).  See QUERYHEADER for a transposed version of the header value
%     table.
%
%     COMPAREHEADER(DATA,'FIELD1',...,'FIELDN') prints out the header
%     fields FIELD1 to FIELDN for records in DATA as a table.  FIELDS may
%     be normal fields ('b' 'kt1' 'xmaximum' etc), group fields ('t' 'kt'
%     etc), absolute fields ('t9 utc' 'user3 tai' 'resp utc' etc), list
%     fields ('picks', 'all', 'full' -- these are not available to the
%     functions CHANGEHEADER or GETHEADER) or wildcards ('*t1' '?' etc).
%     Only * and ? are valid wildcard characters.
%
%    Notes:
%     - Group fields:    T, KT, USER, KUSER, RESP, DEP, ST, EV, NZ, NZDTTM,
%                        KNAME, CMP, DELAZ, REAL, INT, ENUM, LGC, CHAR
%     - List fields:     PICKS, ALL, FULL
%     - Virtual fields:  NZMONTH, NZCDAY, KZDTTM, KZDATE, KZTIME, Z, Z6,
%                        ZTAI, ZTAI6, GCP, NCMP
%     - AbsTime fields:  <FIELD> UTC, <FIELD> TAI
%
%    Examples:
%     % Some simple cases:
%     compareheader(data)          % compare all header variables
%     compareheader(data,'t')      % compare t group
%
%     % Fields are case independent:
%     compareheader(data,'dEltA')
%     compareheader(data,'StLA','stLo')
%
%     % Compare picks:
%     compareheader(data,'picks')
%
%    See also: QUERYHEADER, LISTHEADER, GETHEADER, CHANGEHEADER, VGRP

%     Version History:
%        Sep. 11, 2009 - initial version
%        Sep. 12, 2009 - support vgrp addition, use regexptranslate
%        Sep. 13, 2009 - added utc/tai abs time fields, added abs time
%                        vgrp, added vf support, global option to set
%                        column width, vf show up in wildcards
%        Sep. 15, 2009 - doc update
%        Sep. 18, 2009 - 2nd pass at abs time support
%        Oct.  6, 2009 - dropped use of LOGICAL function
%        Jan. 30, 2010 - use VF_GH_Z to get reference time
%        Feb. 24, 2010 - doc update for compareheader2
%        Mar. 12, 2010 - fixes for Octave (explicit celling)
%        July 30, 2010 - update for compareheader2=>queryheader change
%        Aug. 21, 2010 - minor printing fixes
%        Jan.  5, 2011 - improved H1 line
%        Jan. 30, 2012 - doc update, drop record listing, utc6/tai6 changes
%        Feb.  1, 2012 - utc/tai bugfix
%        Feb.  7, 2012 - edit h1 line
%        Feb. 11, 2012 - divergent octave behavior workaround (cellcat)
%        Feb. 20, 2012 - properly allocate arrays causing above issue
%        Aug. 30, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 30, 2012 at 01:40 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% headers setup
[h,idx]=versioninfo(data);
nh=numel(h);

% number of records
nrecs=numel(data);

% allow changing column width via global
global SEIZMO
try
    cn=SEIZMO.COMPAREHEADER.COLUMNWIDTH;
    cs=num2str(cn);
catch
    % define default column width
    cn=35; cs=num2str(cn);
end

% check column width
if(~isscalar(cn) || cn~=round(cn) || cn<=0)
    error('seizmo:compareheader:badColumnWidth',...
        'COLUMNWIDTH must be a scalar integer >0!');
end

% gather all possible fields and reftimes
fields=cell(2,5,nh); [fields{:}]=deal({}); % size is just a guess
absf=fields; lgcf=fields; enumf=fields;
ref=nan(nrecs,5); good=false(nrecs,1); vf=cell(nh,1); 
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
                absf(k,j,i,1)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' utc')};
                absf(k,j,i,2)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' tai')};
                absf(k,j,i,3)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' utc6')};
                absf(k,j,i,4)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' tai6')};
            elseif(strcmp(h(i).types{j},'lgc'))
                lgcf(k,j,i)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' lgc')};
            elseif(strcmp(h(i).types{j},'enum'))
                enumf(k,j,i,1)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' id')};
                enumf(k,j,i,2)=...
                    {strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' desc')};
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
absf=unique([absf{:}].');
lgcf=unique([lgcf{:}].');
enumf=unique([enumf{:}].');

% list all if no fields given
if(nargin==1); varargin=fields; end

% table header
fprintf('%16s\n','\    ')
disp([sprintf('%16s','FIELD    \   ') '   RECORD NUMBER'])
disp([sprintf('%16s','\  ') sprintf(['%' cs 'd'],1:nrecs)])
disp(char(45*ones(1,16+cn*nrecs))) % 45=dash

% loop over fields
nvarg=numel(varargin);
for i=1:nvarg
    % force lowercase
    gf=strtrim(lower(varargin{i}));
    wf=getwords(gf);
    
    % skip empty
    if(isempty(gf)); continue; end
    
    % require certain modifiers
    if(numel(wf)>1 && ~any(strcmp(wf{2},...
            {'utc' 'tai' 'utc6' 'tai6' 'lgc' 'id' 'desc'})))
        error('seizmo:compareheader:badInput',...
            'Unknown field modifier: %s !',wf{2});
    end
    
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
        % only show modified fields if explicitly sought
        if(numel(wf)==1)
            % no partial matches
            ngf=['^' regexptranslate('wildcard',gf) '$'];
            % find matches in normal fields
            ngf=fields(~cellfun('isempty',regexp(fields,ngf)));
        elseif(any(strcmp(wf(2),{'utc' 'tai' 'utc6' 'tai6'})))
            % no partial matches, trim excess blanks between field/modifier
            ngf=['^' regexptranslate('wildcard',...
                strtrim(wf{1})) ' ' wf{2} '$'];
            % find matches
            ngf=absf(~cellfun('isempty',regexp(absf,ngf)));
        elseif(any(strcmp(wf(2),'lgc')))
            % no partial matches, trim excess blanks between field/modifier
            ngf=['^' regexptranslate('wildcard',...
                wf{1}) ' ' wf{2} '$'];
            % find matches
            ngf=lgcf(~cellfun('isempty',regexp(lgcf,ngf)));
        elseif(any(strcmp(wf(2),{'id' 'desc'})))
            % no partial matches, trim excess blanks between field/modifier
            ngf=['^' regexptranslate('wildcard',...
                wf{1}) ' ' wf{2} '$'];
            % find matches
            ngf=enumf(~cellfun('isempty',regexp(enumf,ngf)));
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
        
        % get value formatted as string (right justified)
        values=nan(nrecs,cn);
        for k=1:nrecs
            tmp=cmph_disp(h(idx(k)),f,data(k),cs,cn,ref(k,:),good(k));
            
            % check for oversized output
            if(numel(tmp)>cn)
                % replace oversized output with *
                tmp=[];
                tmp(1)=' ';
                tmp(2:cn)='*';
            end
            values(k,:)=tmp;
        end
        
        % display
        values=values.';
        disp([sprintf('%13s | ',upper(f)) char(values(:).')])
    end
end

end

function [string]=cmph_disp(h,f,data,cs,cn,reftime,good)
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
wf=getwords(f); solo=numel(wf)==1;
for m=1:numel(h.enum)
    if(isfield(h.enum(m).pos,wf{1}))
        ival=data.head(h.enum(m).pos.(wf{1}));
        % searches only enum(1) for now (prefer one big
        % set of enums to share vs separate sets)
        if(ival==fix(ival) && ival>=h.enum(1).minval && ...
                ival<=h.enum(1).maxval)
            if(solo)
                string=sprintf(['%' cs 's'],h.enum(1).id{ival+1});
            else
                if(strcmp(wf{2},'id'))
                    string=sprintf(['%' cs 's'],h.enum(1).id{ival+1});
                else % desc
                    string=sprintf(['%' cs 's'],h.enum(1).desc{ival+1});
                end
            end
        elseif(ival==h.undef.ntype)
            string=sprintf(['%' cs 's'],sprintf('UNDEFINED (%g)',ival));
        else
            string=sprintf(['%' cs 's'],sprintf('UNKNOWN (%g)',ival));
        end
        return;
    end
end
for m=1:numel(h.lgc)
    if(isfield(h.lgc(m).pos,wf{1}))
        if(data.head(h.lgc(m).pos.(wf{1}))==h.false)
            string=sprintf(['%' cs 's'],'FALSE');
        elseif(data.head(h.lgc(m).pos.(wf{1}))==h.true)
            string=sprintf(['%' cs 's'],'TRUE');  
        elseif(data.head(h.lgc(m).pos.(wf{1}))==h.undef.ntype)
            string=sprintf(['%' cs 's'],sprintf('UNDEFINED (%g)',...
                data.head(h.lgc(m).pos.(wf{1}))));
        else
            string=sprintf(['%' cs 's'],sprintf('INVALID (%g)',...
                data.head(h.lgc(m).pos.(wf{1}))));
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
if(solo)
    for m=1:numel(h.ntype)
        for n=1:numel(h.(h.ntype{m}))
            if(isfield(h.(h.ntype{m})(n).pos,f))
                if(data.head(h.(h.ntype{m})(n).pos.(f))==h.undef.ntype)
                    string=sprintf(['%' cs 's'],...
                        sprintf('UNDEFINED (%g)',...
                        data.head(h.(h.ntype{m})(n).pos.(f))));
                else
                    string=sprintf(['%' cs 's'],sprintf('%-.10g',...
                        data.head(h.(h.ntype{m})(n).pos.(f))));
                end
                return;
            end
        end
    end
else
    for n=1:numel(h.real)
        % absolute time fields section
        if(isfield(h.real(n).pos,wf{1}))
            v=data.head(h.real(n).pos.(wf{1}));
            if(any(strcmp(wf{2},{'utc' 'utc6' 'tai' 'tai6'})))
                if(v==h.undef.ntype || isnan(v) || isinf(v))
                    string=sprintf(['%' cs 's'],...
                        sprintf('UNDEFINED (%g)',v));
                elseif(~good)
                    string=sprintf(['%' cs 's'],...
                        sprintf('NO REFTIME (%g)',v));
                else
                    % get values for output
                    switch wf{2}
                        case {'utc' 'utc6'}
                            utc=fixtimes(reftime+[0 0 0 0 v],'utc');
                        case {'tai' 'tai6'}
                            utc=utc2tai(reftime+[0 0 0 0 v]);
                    end
                    cal=doy2cal(utc(1:2));
                    utc(5)=round(1000*utc(5));
                    string=sprintf(['%' cs 's'],...
                        sprintf(['%04d-%02d-%02d (%03d) '...
                        '%02d:%02d:%02d.%03d'],utc(1),cal(2),cal(3),...
                        utc(2),utc(3),utc(4),fix(utc(5)/1000),...
                        mod(utc(5),1000)));
                end
                return;
            end
        end
    end
end

% field not found
string=sprintf(['%' cs 's'],'NOT A FIELD!');

end
