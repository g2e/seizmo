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
%        Sep. 12, 2009 - support vgrp addition, use regexptranslate
%        Sep. 13, 2009 - added utc/tai abs time fields, added abs time
%                        vgrp, added vf support, global option to set
%                        column width, vf show up in wildcards
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 13, 2009 at 23:10 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

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
fields=cell(2,5,nh); utcf=fields; taif=fields;
reftime=nan(nrecs,6); vf=cell(nh,1);
for i=1:nh
    % vf in wildcard search
    vf{i}=fieldnames(h(i).vf).';
    
    % all available field type sets
    for j=1:numel(h(i).types)
        for k=1:numel(h(i).(h(i).types{j}))
            fields{k,j,i}=...
                fieldnames(h(i).(h(i).types{j})(k).pos).';
            % special absolute time alternate fields for all real
            if(strcmp(h(i).types{j},'real'))
                utcf{k,j,i}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' utc');
                taif{k,j,i}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' tai');
            end
        end
    end
    
    % get reference times hack
    % (skips getheader)
    vidx=find(idx==i);
    head=[data(vidx).head];
    reftime(vidx,:)=head(h.reftime,:).';
    reftime(vidx(reftime(vidx,:)==h.undef.ntype),:)=nan;
    reftime(:,5)=reftime(:,5)+reftime(:,6)/1000;
    reftime=utc2tai(reftime(:,1:5));
end
% only list each field once
% - this forces an alphabetical listing
fields=unique([[fields{:}].'; [vf{:}].']);
utcf=unique([utcf{:}].');
taif=unique([taif{:}].');

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

% table header
disp(sprintf('%15s','\    '))
disp([sprintf('%15s','FIELD    \   ') '   RECORD NUMBER'])
disp([sprintf('%15s','\  ') sprintf(['%' cs 'd'],1:nrecs)])
disp(char(45*ones(1,15+cn*nrecs)))

% loop over fields
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
        
        % get value formatted as string (right justified)
        values=nan(nrecs,cn);
        for k=1:nrecs
            tmp=cmph_disp(h(idx(k)),f,data(k),cs,cn,reftime(k,:));
            
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
        disp([sprintf('%12s | ',upper(f)) char(values(:).')])
    end
end

end

function [string]=cmph_disp(h,f,data,cs,cn,reftime)
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
        % UTC absolute times section
        elseif(strcmp(f(max(1,end-3):end),' utc'))
            f=strtrim(f(1:end-4));
            if(isfield(h.(h.ntype{m})(n).pos,f))
                if(data.head(h.(h.ntype{m})(n).pos.(f))==h.undef.ntype)
                    string=sprintf(['%' cs 's'],...
                        sprintf('UNDEFINED (%g)',...
                        data.head(h.(h.ntype{m})(n).pos.(f))));
                elseif(any(isnan(reftime) | isinf(reftime)))
                    string=sprintf(['%' cs 's'],...
                        sprintf('NO REFTIME (%g)',...
                        data.head(h.(h.ntype{m})(n).pos.(f))));
                else
                    % get values for output
                    utc=tai2utc(reftime+...
                        [0 0 0 0 data.head(h.(h.ntype{m})(n).pos.(f))]);
                    cal=doy2cal(utc(1:2));
                    utc(5)=round(1000*utc(5));
                    string=sprintf(['%' cs 's'],...
                        sprintf(['%04d-%02d-%02d (%03d) '...
                        '%02d:%02d:%02d.%03d'],utc(1),cal(2),cal(3),...
                        utc(2),utc(3),utc(4),floor(utc(5)/1000),...
                        mod(utc(5),1000)));
                end
                return;
            end
        % TAI absolute times section
        elseif(strcmp(f(max(1,end-3):end),' tai'))
            f=strtrim(f(1:end-4));
            if(isfield(h.(h.ntype{m})(n).pos,f))
                if(data.head(h.(h.ntype{m})(n).pos.(f))==h.undef.ntype)
                    string=sprintf(['%' cs 's'],...
                        sprintf('UNDEFINED (%g)',...
                        data.head(h.(h.ntype{m})(n).pos.(f))));
                elseif(any(isnan(reftime) | isinf(reftime)))
                    string=sprintf(['%' cs 's'],...
                        sprintf('NO REFTIME (%g)',...
                        data.head(h.(h.ntype{m})(n).pos.(f))));
                else
                    % get values for output
                    tai=fixtimes(reftime+...
                        [0 0 0 0 data.head(h.(h.ntype{m})(n).pos.(f))]);
                    cal=doy2cal(tai(1:2));
                    tai(5)=round(1000*tai(5));
                    string=sprintf(['%' cs 's'],...
                        sprintf(['%04d-%02d-%02d (%03d) '...
                        '%02d:%02d:%02d.%03d'],tai(1),cal(2),cal(3),...
                        tai(2),tai(3),tai(4),floor(tai(5)/1000),...
                        mod(tai(5),1000)));
                end
                return;
            end
        end
    end
end

% field not found
string=sprintf(['%' cs 's'],'NOT A FIELD!');

end
