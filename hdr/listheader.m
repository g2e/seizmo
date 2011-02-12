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
%     valid header field, group field (ie. 't' 'kt' 'resp' 'user' 'kuser'
%     etc), absolute fields ('t9 utc' 'user3 tai' 'resp utc' etc), or
%     wildcards ('nz*' 'dep*' etc).  Only * and ? are valid wildcards.
%
%    Notes:
%     - group fields:    t, kt, user, kuser, resp, dep, st, ev, nz, nzdttm,
%                         kname, {real group} utc, {real group} tai
%     - list fields:     picks, all, full
%     - virtual fields:  nzmonth, nzcday, kzdttm, kzdate, kztime, z, ztai
%     - abs time fields: {real field} utc, {real field} tai
%
%    Examples:
%      listheader(data)          % lists all header variables
%      listheader(data,'t')      % lists t group
%
%     Fields are case independent:
%      listheader(data,'dEltA')
%      listheader(data,'StLA','stLo')
%
%     List only single char fields
%      listheader(data,'?')
%
%    See also:  COMPAREHEADER, CHANGEHEADER, GETHEADER

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
%        Sep. 12, 2009 - added vgrp support, use regexptranslate
%        Sep. 13, 2009 - added utc/tai abs time fields
%        Sep. 14, 2009 - vlists, abs time vgrp, vf, vf via wildcards
%        Sep. 15, 2009 - doc updated
%        Sep. 18, 2009 - 2nd pass at abs time support
%        Jan. 30, 2010 - use VF_GH_Z to get reference time
%        Feb. 11, 2011 - mass nargchk fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:
% - skip undefined (set via global)

% check nargin
error(nargchk(1,inf,nargin));

% headers setup
[h,idx]=versioninfo(data);
nh=numel(h);

% number of records
nrecs=numel(data);

% extra padding
disp(' ')

% gather all possible fields and reftimes
nfields=cell(nh,1);
nutcf=nfields; ntaif=nfields; nutc6f=nfields; ntai6f=nfields;
tmpfields=cell(2,5); % just a guess
tmputcf=tmpfields; tmptaif=tmpfields;
tmputc6f=tmpfields; tmptai6f=tmpfields;
ref=nan(nrecs,5); good=false(nrecs,1);
for i=1:nh
    % add virtual fields to wildcard search
    vf=fieldnames(h(i).vf);
    
    % all available field type sets
    for j=1:numel(h(i).types)
        for k=1:numel(h(i).(h(i).types{j}))
            tmpfields{k,j}=...
                fieldnames(h(i).(h(i).types{j})(k).pos).';
            % special absolute time alternate fields for all real
            if(strcmp(h(i).types{j},'real'))
                tmputcf{k,j}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' utc');
                tmptaif{k,j}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' tai');
                tmputc6f{k,j}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' 6utc');
                tmptai6f{k,j}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' 6tai');
            end
        end
    end
    
    % combine
    nfields{i}=[[tmpfields{:}].'; vf];
    nutcf{i}=[tmputcf{:}].';
    ntaif{i}=[tmptaif{:}].';
    nutc6f{i}=[tmputc6f{:}].';
    ntai6f{i}=[tmptai6f{:}].';
    tmpfields=cell(2,5);
    tmputcf=tmpfields;
    tmptaif=tmpfields;
    tmputc6f=tmpfields;
    tmptai6f=tmpfields;
    
    % get reference times hack
    vidx=find(idx==i);
    [ref(vidx,:),good(vidx,1)]=vf_gh_z(h(i),[data(vidx).head]);
end

% loop over files
for i=1:nrecs
    % all available fields
    fields=nfields{idx(i)};
    utcf=nutcf{idx(i)};
    taif=ntaif{idx(i)};
    utc6f=nutc6f{idx(i)};
    tai6f=ntai6f{idx(i)};
    
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
        gf=strtrim(lower(varargin{j}));
        wf=getwords(gf);
        
        % skip empty
        if(isempty(gf)); continue; end
        
        % check for vlist/vgrp
        group=false; glen=1;
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
            if(isfield(h(idx(i)).vgrp,wf{1}))
                ngf=strtrim(strcat(h(idx(i)).vgrp.(wf{1}),{' '},...
                    joinwords(wf(2:end))));
                group=true; glen=numel(ngf);
            end
        end
        
        % wildcard case (?==63,*==42) - pass to regexp
        wild=false;
        if(~group && (any(gf==42) || any(gf==63)))
            % declare as wildcard
            wild=true;
            % only show absolute time fields if explicitly sought
            if(strcmp(joinwords(wf(2:end)),'utc'))
                % no partial matches, trim excess blanks between field, utc
                ngf=['^' regexptranslate('wildcard',...
                    strtrim(wf{1})) ' utc$'];
                % find matches in utc set
                ngf=utcf(~cellfun('isempty',regexp(utcf,ngf)));
            elseif(strcmp(joinwords(wf(2:end)),'tai'))
                % no partial matches, trim excess blanks between field, tai
                ngf=['^' regexptranslate('wildcard',...
                    strtrim(wf{1})) ' tai$'];
                % find matches in tai set
                ngf=taif(~cellfun('isempty',regexp(taif,ngf)));
            elseif(strcmp(joinwords(wf(2:end)),'6utc'))
                % no partial matches, trim excess blanks between field, utc
                ngf=['^' regexptranslate('wildcard',...
                    strtrim(wf{1})) ' 6utc$'];
                % find matches in utc set
                ngf=utc6f(~cellfun('isempty',regexp(utc6f,ngf)));
            elseif(strcmp(joinwords(wf(2:end)),'6tai'))
                % no partial matches, trim excess blanks between field, tai
                ngf=['^' regexptranslate('wildcard',...
                    strtrim(wf{1})) ' 6tai$'];
                % find matches in tai set
                ngf=tai6f(~cellfun('isempty',regexp(tai6f,ngf)));
            else
                % no partial matches
                ngf=['^' regexptranslate('wildcard',gf) '$'];
                % find matches in normal fields
                ngf=fields(~cellfun('isempty',regexp(fields,ngf)));
            end
            glen=numel(ngf);
        end
        
        % group loop
        for k=1:glen
            % modify field if in a group
            if(group); f=ngf{k};
            elseif(wild); f=ngf{k};
            else f=gf;
            end
            
            % display field/value
            lh_disp(h(idx(i)),f,data(i),ref(i,:),good(i));
        end
    end
end

end

function []=lh_disp(h,f,data,reftime,good)
%LH_DISP    Finds and displays field/value pairs for lh

% virtual fields
if(isfield(h.vf,f))
    switch h.vf.(f).type
        case 'enum'
            % searches only enum(1) for now (prefer one big
            % set of enums to share vs separate sets)
            ival=h.vf.(f).gh(h,data.head);
            if(ival==fix(ival) && ival>=h.enum(1).minval && ...
                    ival<=h.enum(1).maxval)
                disp(sprintf('%17s = %s',upper(f),h.enum(1).id{ival+1}));
            elseif(ival==h.undef.ntype)
                disp(sprintf('%17s = UNDEFINED (%g)',upper(f),ival));
            else
                disp(sprintf('%17s = UNKNOWN (%g)',upper(f),ival));
            end
        case 'lgc'
            ival=h.vf.(f).gh(h,data.head);
            switch ival
                case h.false
                    disp(sprintf('%17s = FALSE',upper(f)));
                case h.true
                    disp(sprintf('%17s = TRUE',upper(f)));
                case h.undef.ntype
                    disp(sprintf('%17s = UNDEFINED (%g)',upper(f),ival));
                otherwise
                    disp(sprintf('%17s = INVALID (%g)',upper(f),ival));
            end
        case 'char'
            ival=h.vf.(f).gh(h,data.head);
            disp(sprintf('%17s = %s',upper(f),ival{:}));
        case 'abs'
            ival=h.vf.(f).lh(h,data.head);
            disp(sprintf('%17s = %s',upper(f),ival{:}));
        otherwise
            ival=h.vf.(f).gh(h,data.head);
            if(ival==h.undef.ntype)
                disp(sprintf('%17s = UNDEFINED (%g)',upper(f),ival));
            else
                disp(sprintf('%17s = %-.10g',upper(f),ival));
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
            disp(sprintf('%17s = %s',upper(f),h.enum(1).desc{ival+1}));
        elseif(ival==h.undef.ntype)
            disp(sprintf('%17s = UNDEFINED (%g)',upper(f),ival));
        else
            disp(sprintf('%17s = UNKNOWN (%g)',upper(f),ival));
        end
        return;
    end
end
for m=1:numel(h.lgc)
    if(isfield(h.lgc(m).pos,f))
        if(data.head(h.lgc(m).pos.(f))==h.false)
            disp(sprintf('%17s = FALSE',upper(f)));
        elseif(data.head(h.lgc(m).pos.(f))==h.true)
            disp(sprintf('%17s = TRUE',upper(f)));  
        elseif(data.head(h.lgc(m).pos.(f))==h.undef.ntype)
            disp(sprintf('%17s = UNDEFINED (%g)',upper(f),...
                data.head(h.lgc(m).pos.(f))));
        else
            disp(sprintf('%17s = INVALID (%g)',upper(f),...
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
                disp(sprintf('%17s = UNDEFINED (%s)',upper(f),...
                    char(data.head(p(1):p(2)))));
            else
                disp(sprintf('%17s = %s',upper(f),...
                    char(data.head(p(1):p(2)))));
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
                disp(sprintf('%17s = UNDEFINED (%g)',upper(f),...
                    data.head(h.(h.ntype{m})(n).pos.(f))));
            else
                disp(sprintf('%17s = %-.10g',upper(f),...
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
                        disp(sprintf('%17s = UNDEFINED (%g)',upper(f),v));
                    elseif(~good)
                        disp(sprintf('%17s = NO REFTIME (%g)',upper(f),v));
                    else
                        % get values for output
                        utc=tai2utc(reftime+[0 0 0 0 v]);
                        cal=doy2cal(utc(1:2));
                        utc(5)=round(1000*utc(5));
                        disp(sprintf(['%17s = %04d-%02d-%02d (%03d) '...
                            '%02d:%02d:%02d.%03d'],upper(f),utc(1),...
                            cal(2),cal(3),utc(2),utc(3),utc(4),...
                            fix(utc(5)/1000),mod(utc(5),1000)));
                    end
                    return;
                elseif(any(strcmpi(joinwords(wf(2:end)),{'tai' '6tai'})))
                    if(v==h.undef.ntype || isnan(v) || isinf(v))
                        disp(sprintf('%17s = UNDEFINED (%g)',upper(f),v));
                    elseif(~good)
                        disp(sprintf('%17s = NO REFTIME (%g)',upper(f),v));
                    else
                        tai=fixtimes(reftime+[0 0 0 0 v]);
                        cal=doy2cal(tai(1:2));
                        tai(5)=round(1000*tai(5));
                        disp(sprintf(['%17s = %04d-%02d-%02d (%03d) '...
                            '%02d:%02d:%02d.%03d'],upper(f),tai(1),...
                            cal(2),cal(3),tai(2),tai(3),tai(4),...
                            fix(tai(5)/1000),mod(tai(5),1000)));
                    end
                    return;
                end
            end
        end
    end
end

% field not found
%warning('seizmo:listheader:fieldInvalid',...
%   'Filetype: %s, Version: %d\nInvalid field: %s',h.filetype,h.version,f);

end
