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
%     valid header field, group field (ie. 't' 'kt' 'resp' 'user' etc), or
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
%        Sep. 12, 2009 - added vgrp support, use regexptranslate
%        Sep. 13, 2009 - added utc/tai abs time fields
%        Sep. 14, 2009 - vlists, abs time vgrp, vf, vf via wildcards
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2009 at 08:20 GMT

% todo:
% - skip undefined (set via global)
% - multiple columns? (just like cmph)

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% headers setup
[h,idx]=versioninfo(data);
nh=numel(h);

% number of records
nrecs=numel(data);

% extra padding
disp(' ')

% gather all possible fields and reftimes
nfields=cell(nh,1); nutcf=nfields; ntaif=nfields;
tmpfields=cell(2,5); tmputcf=tmpfields; tmptaif=tmpfields;
reftime=nan(nrecs,6);
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
            end
        end
    end
    
    % combine
    nfields{i}=[[tmpfields{:}].'; vf];
    nutcf{i}=[tmputcf{:}].';
    ntaif{i}=[tmptaif{:}].';
    tmpfields=cell(2,5);
    tmputcf=tmpfields;
    tmptaif=tmpfields;
    
    % get reference times hack
    % (skips getheader)
    vidx=find(idx==i);
    head=[data(vidx).head];
    reftime(vidx,:)=head(h.reftime,:).';
    reftime(vidx(reftime(vidx,:)==h.undef.ntype),:)=nan;
    reftime(:,5)=reftime(:,5)+reftime(:,6)/1000;
    reftime=utc2tai(reftime(:,1:5));
end

% loop over files
for i=1:nrecs
    % all available fields
    fields=nfields{idx(i)};
    utcf=nutcf{idx(i)};
    taif=ntaif{idx(i)};
    
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
            lh_disp(h(idx(i)),f,data(i),reftime(i,:));
        end
    end
end

end

function []=lh_disp(h,f,data,reftime)
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
                disp(sprintf('%12s = %s',upper(f),h.enum(1).id{ival+1}));
            elseif(ival==h.undef.ntype)
                disp(sprintf('%12s = UNDEFINED (%g)',upper(f),ival));
            else
                disp(sprintf('%12s = UNKNOWN (%g)',upper(f),ival));
            end
        case 'lgc'
            ival=h.vf.(f).gh(h,data.head);
            switch ival
                case h.false
                    disp(sprintf('%12s = FALSE',upper(f)));
                case h.true
                    disp(sprintf('%12s = TRUE',upper(f)));
                case h.undef.ntype
                    disp(sprintf('%12s = UNDEFINED (%g)',upper(f),ival));
                otherwise
                    disp(sprintf('%12s = INVALID (%g)',upper(f),ival));
            end
        case 'char'
            ival=h.vf.(f).gh(h,data.head);
            disp(sprintf('%12s = %s',upper(f),ival{:}));
        otherwise
            ival=h.vf.(f).gh(h,data.head);
            if(ival==h.undef.ntype)
                disp(sprintf('%12s = UNDEFINED (%g)',upper(f),ival));
            else
                disp(sprintf('%12s = %-.10g',upper(f),ival));
            end
    end
    return;
end

% handle enum/lgc fields
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
        % UTC absolute times section
        elseif(strcmp(f(max(1,end-3):end),' utc'))
            f2=strtrim(f(1:end-4));
            if(isfield(h.(h.ntype{m})(n).pos,f2))
                if(data.head(h.(h.ntype{m})(n).pos.(f2))==h.undef.ntype)
                    disp(sprintf('%12s = UNDEFINED (%g)',upper(f),...
                        data.head(h.(h.ntype{m})(n).pos.(f2))));
                elseif(any(isnan(reftime) | isinf(reftime)))
                    disp(sprintf('%12s = NO REFTIME (%g)',upper(f),...
                        data.head(h.(h.ntype{m})(n).pos.(f2))));
                else
                    % get values for output
                    utc=tai2utc(reftime+...
                        [0 0 0 0 data.head(h.(h.ntype{m})(n).pos.(f2))]);
                    cal=doy2cal(utc(1:2));
                    utc(5)=round(1000*utc(5));
                    disp(sprintf(['%12s = %04d-%02d-%02d (%03d) '...
                        '%02d:%02d:%02d.%03d'],upper(f),utc(1),cal(2),...
                        cal(3),utc(2),utc(3),utc(4),floor(utc(5)/1000),...
                        mod(utc(5),1000)));
                end
                return;
            end
        % TAI absolute times section
        elseif(strcmp(f(max(1,end-3):end),' tai'))
            f2=strtrim(f(1:end-4));
            if(isfield(h.(h.ntype{m})(n).pos,f2))
                if(data.head(h.(h.ntype{m})(n).pos.(f2))==h.undef.ntype)
                    disp(sprintf('%12s = UNDEFINED (%g)',upper(f),...
                        data.head(h.(h.ntype{m})(n).pos.(f2))));
                elseif(any(isnan(reftime) | isinf(reftime)))
                    disp(sprintf('%12s = NO REFTIME (%g)',upper(f),...
                        data.head(h.(h.ntype{m})(n).pos.(f2))));
                else
                    % get values for output
                    tai=fixtimes(reftime+...
                        [0 0 0 0 data.head(h.(h.ntype{m})(n).pos.(f2))]);
                    cal=doy2cal(tai(1:2));
                    tai(5)=round(1000*tai(5));
                    disp(sprintf(['%12s = %04d-%02d-%02d (%03d) '...
                        '%02d:%02d:%02d.%03d'],upper(f),tai(1),cal(2),...
                        cal(3),tai(2),tai(3),tai(4),floor(tai(5)/1000),...
                        mod(tai(5),1000)));
                end
                return;
            end
        end
    end
end

% field not found
warning('seizmo:listheader:fieldInvalid',...
    'Filetype: %s, Version: %d\nInvalid field: %s',h.filetype,h.version,f);

end
