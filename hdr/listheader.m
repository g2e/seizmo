function []=listheader(data,varargin)
%LISTHEADER    List SEIZMO data headers
%
%    Usage:    listheader(data)
%              listheader(data,'field1',...,'fieldN')
%
%    Description:
%     LISTHEADER(DATA) lists the header fields of all records in DATA in a
%     format similar to SAC's lh command.  Undefined fields are skipped.
%
%     LISTHEADER(DATA,'FIELD1',...,'FIELDN') prints a list of the header
%     field(s) FIELD1 to FIELDN and their value(s) from records in DATA
%     like SAC's lh command does.  FIELDS may be normal fields ('b' 'kt1'
%     'xmaximum' etc), group fields ('t' 'kt' etc), absolute fields
%     ('t9 utc' 'user3 tai' 'resp utc' etc), list fields ('picks' 'all'
%     'full' -- these are not available to the functions CHANGEHEADER or
%     GETHEADER) or wildcards ('*t1' '?' etc).  Only * and ? are valid
%     wildcard characters.
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
%     % List some timing info:
%     listheader(data,'b','e','delta','npts');
%
%     % List the t group field values:
%     listheader(data,'t')
%
%     % Fields are case independent:
%     listheader(data,'dEltA')
%     listheader(data,'StLA','stLo')
%
%     % List only single char fields
%     listheader(data,'?')
%
%     % List only the logical fields:
%     listheader(data,'lgc')
%
%    See also: QUERYHEADER, COMPAREHEADER, CHANGEHEADER, GETHEADER, VGRP

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
%        Jan. 30, 2012 - doc update, fix disp+sprintf, utc6/tai6 changes,
%                        skip undefined values, sac-like default list,
%                        multi-column support
%        Feb.  1, 2012 - fix tai/tai6 output, adaptive column spacing,
%                        split path/file output, utc/tai bugfix
%        Feb. 20, 2012 - clean up cell array initialization issues
%        Aug. 30, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 30, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% hide undefined fields?
global SEIZMO
try hide=SEIZMO.LISTHEADER.HIDE; catch; hide=true; end
try cols=SEIZMO.LISTHEADER.COLUMNS; catch; cols=2; end

% SAC's default field list
defaults={'npts' 'b' 'e' 'iftype' 'leven' 'delta' 'odelta' 'idep' 'dep' ...
    'ko' 'o' 'ka' 'a' 'picks' 'kf' 'f' 'kzdate' 'kztime' 'iztype' ...
    'kinst' 'resp' 'kdatrd' 'kstnm' 'cmp' 'istreg' 'st' 'kevnm' ...
    'ievreg' 'ev' 'ievtyp' 'khole' 'dist' 'az' 'baz' 'gcarc' 'lovrok' ...
    'iqual' 'isynth' 'user' 'kuser' 'nxsize' 'xminimum' 'xmaximum' ...
    'nysize' 'yminimum' 'ymaximum' 'nvhdr' 'scale' 'norid' 'nevid' ...
    'nwfid' 'iinst' 'lpspol' 'lcalda' 'kcmpnm' 'knetwk' 'mag' 'imagtyp' ...
    'imagsrc'};

% headers setup
[h,idx]=versioninfo(data);
nh=numel(h);

% number of records
nrecs=numel(data);

% gather all possible fields and reftimes
nfields=cell(nh,1);
nabsf=nfields; nlgcf=nfields; nenumf=nfields;
tmpfields=cell(2,5); [tmpfields{:}]=deal({}); % size is just a guess
tmpabsf=tmpfields; tmplgcf=tmpfields; tmpenumf=tmpfields;
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
                tmpabsf{k,j,1}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' utc');
                tmpabsf{k,j,2}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' tai');
                tmpabsf{k,j,3}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' utc6');
                tmpabsf{k,j,4}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' tai6');
            elseif(strcmp(h(i).types{j},'lgc'))
                tmplgcf{k,j,1}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' lgc');
            elseif(strcmp(h(i).types{j},'enum'))
                tmpenumf{k,j,1}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' id');
                tmpenumf{k,j,2}=...
                    strcat(fieldnames(h(i).(h(i).types{j})(k).pos).',...
                    ' desc');
            end
        end
    end
    
    % combine
    nfields{i}=[[tmpfields{:}].'; vf];
    nabsf{i}=[tmpabsf{:}].';
    nlgcf{i}=[tmplgcf{:}].';
    nenumf{i}=[tmpenumf{:}].';
    tmpfields=cell(2,5);
    tmpabsf=tmpfields;
    tmplgcf=tmpfields;
    tmpenumf=tmpfields;
    
    % get reference times hack
    vidx=find(idx==i);
    [ref(vidx,:),good(vidx,1)]=vf_gh_z(h(i),[data(vidx).head]);
end

% loop over files
for i=1:nrecs
    % all available fields
    fields=nfields{idx(i)};
    absf=nabsf{idx(i)};
    lgcf=nlgcf{idx(i)};
    enumf=nenumf{idx(i)};
    
    % list all case
    if(nargin==1); varargin=defaults; end
    
    % formatted header
    fprintf('\n PATH: %s\n',data(i).path)
    fprintf(' FILE: (%d) %s\n',i,data(i).name)
    
    % setup out
    cnt=1; out=cell(cnt,1);
    
    % loop over fields
    for j=1:numel(varargin)
        % force lowercase
        gf=strtrim(lower(varargin{j}));
        wf=getwords(gf);
        
        % skip empty
        if(isempty(gf)); continue; end
        
        % require certain modifiers
        if(numel(wf)>1 && ~any(strcmp(wf{2},...
                {'utc' 'tai' 'utc6' 'tai6' 'lgc' 'id' 'desc'})))
            error('seizmo:listheader:badInput',...
                'Unknown field modifier: %s !',wf{2});
        end
        
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
        
        % group loop
        for k=1:glen
            % modify field if in a group
            if(group); f=ngf{k};
            elseif(wild); f=ngf{k};
            else f=gf;
            end
            
            % get formatted field/value
            out{cnt}=lh_format(h(idx(i)),f,data(i),ref(i,:),good(i),hide);
            if(isempty(out{cnt})); continue; end
            cnt=cnt+1;
        end
    end
    
    % display field/value
    if(rem(numel(out),cols)>0)
        out(end+1:end+cols-rem(numel(out),cols))={''};
    end
    rows=numel(out)/cols;
    out=mat2cell(reshape(out,cols,[])',rows,ones(1,cols));
    for j=1:cols; out{j}=[char(32*ones(rows,2)) strtrim(char(out{j}))]; end
    out=[out{:}];
    disp(char(45*ones(1,size(out,2))));
    disp(out);
end

end

function [out]=lh_format(h,f,data,reftime,good,hide)
%LH_FORMAT    Finds and formats field/value pairs for lh display

% preset out to empty
out='';

% virtual fields
if(isfield(h.vf,f))
    switch h.vf.(f).type
        case 'enum'
            % searches only enum(1) for now (prefer one big
            % set of enums to share vs separate sets)
            ival=h.vf.(f).gh(h,data.head);
            if(ival==fix(ival) && ival>=h.enum(1).minval && ...
                    ival<=h.enum(1).maxval)
                out=sprintf('%17s = %s',upper(f),h.enum(1).id{ival+1});
            elseif(ival==h.undef.ntype)
                if(hide); return; end
                out=sprintf('%17s = UNDEFINED (%g)',upper(f),ival);
            else
                out=sprintf('%17s = UNKNOWN (%g)',upper(f),ival);
            end
        case 'lgc'
            ival=h.vf.(f).gh(h,data.head);
            switch ival
                case h.false
                    out=sprintf('%17s = FALSE',upper(f));
                case h.true
                    out=sprintf('%17s = TRUE',upper(f));
                case h.undef.ntype
                    if(hide); return; end
                    out=sprintf('%17s = UNDEFINED (%g)',upper(f),ival);
                otherwise
                    out=sprintf('%17s = INVALID (%g)',upper(f),ival);
            end
        case 'char'
            ival=h.vf.(f).gh(h,data.head);
            if(hide && strcmp(ival,h.undef.stype)); return; end
            out=sprintf('%17s = %s',upper(f),ival{:});
        case 'abs'
            ival=h.vf.(f).lh(h,data.head);
            if(hide && strcmp(ival,h.undef.stype)); return; end
            out=sprintf('%17s = %s',upper(f),ival{:});
        otherwise
            ival=h.vf.(f).gh(h,data.head);
            if(ival==h.undef.ntype)
                if(hide); return; end
                out=sprintf('%17s = UNDEFINED (%g)',upper(f),ival);
            else
                out=sprintf('%17s = %-.10g',upper(f),ival);
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
                out=sprintf('%17s = %s',upper(f),h.enum(1).desc{ival+1});
            else
                if(strcmp(wf{2},'id'))
                    out=sprintf('%17s = %s',...
                        upper(f),h.enum(1).id{ival+1});
                else % desc
                    out=sprintf('%17s = %s',...
                        upper(f),h.enum(1).desc{ival+1});
                end
            end
        elseif(ival==h.undef.ntype)
            if(hide); return; end
            out=sprintf('%17s = UNDEFINED (%g)',upper(f),ival);
        else
            out=sprintf('%17s = UNKNOWN (%g)',upper(f),ival);
        end
        return;
    end
end
for m=1:numel(h.lgc)
    if(isfield(h.lgc(m).pos,wf{1}))
        if(data.head(h.lgc(m).pos.(wf{1}))==h.false)
            out=sprintf('%17s = FALSE',upper(f));
        elseif(data.head(h.lgc(m).pos.(wf{1}))==h.true)
            out=sprintf('%17s = TRUE',upper(f));  
        elseif(data.head(h.lgc(m).pos.(wf{1}))==h.undef.ntype)
            if(hide); return; end
            out=sprintf('%17s = UNDEFINED (%g)',upper(f),...
                data.head(h.lgc(m).pos.(wf{1})));
        else
            out=sprintf('%17s = INVALID (%g)',upper(f),...
                data.head(h.lgc(m).pos.(wf{1})));
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
                if(hide); return; end
                out=sprintf('%17s = UNDEFINED (%s)',upper(f),...
                    char(data.head(p(1):p(2))));
            else
                out=sprintf('%17s = %s',upper(f),...
                    char(data.head(p(1):p(2))));
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
                    if(hide); return; end
                    out=sprintf('%17s = UNDEFINED (%g)',upper(f),...
                        data.head(h.(h.ntype{m})(n).pos.(f)));
                else
                    out=sprintf('%17s = %-.10g',upper(f),...
                        data.head(h.(h.ntype{m})(n).pos.(f)));
                end
                return;
            end
        end
    end
else % multi
    for n=1:numel(h.real)
        % absolute time fields section
        if(isfield(h.real(n).pos,wf{1}))
            v=data.head(h.real(n).pos.(wf{1}));
            if(any(strcmp(wf{2},{'utc' 'utc6' 'tai' 'tai6'})))
                if(v==h.undef.ntype || isnan(v) || isinf(v))
                    if(hide); return; end
                    out=sprintf('%17s = UNDEFINED (%g)',upper(f),v);
                elseif(~good)
                    out=sprintf('%17s = NO REFTIME (%g)',upper(f),v);
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
                    out=sprintf(['%17s = %04d-%02d-%02d (%03d) '...
                        '%02d:%02d:%02d.%03d'],upper(f),utc(1),...
                        cal(2),cal(3),utc(2),utc(3),utc(4),...
                        fix(utc(5)/1000),mod(utc(5),1000));
                end
                return;
            end
        end
    end
end

% field not found
%warning('seizmo:listheader:fieldInvalid',...
%   'Filetype: %s, Version: %d\nInvalid field: %s',h.filetype,h.version,f);

end
