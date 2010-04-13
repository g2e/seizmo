function [data]=changeheader(data,varargin)
%CHANGEHEADER    Change SEIZMO data header values
%
%    Usage: data=changeheader(data,'field1',values1)
%           data=changeheader(data,'field1,'values1,...,'fieldN',valuesN)
%
%    Description: CHANGEHEADER(DATA,FIELD,VALUE) changes the header field 
%     FIELD to the value(s) in VALUE for each record in DATA and returns
%     an updated SEIZMO data structure.  FIELD must be a string 
%     corresponding to a valid header field.  A VALUE for normal header
%     fields may be a scalar (assigns same value to all), a column vector
%     of length equal to the number of records (assigns a separate value to
%     each record), or a char array (which is converted to a cellstr array
%     internally, and thus becomes a column vector).  YOU SHOULD NOT PASS A
%     ROW VECTOR UNLESS ALL YOUR RECORDS ARE A SINGLE FILETYPE!  If the
%     FIELD is a group field ('t','resp','user','kt','kuser') then VALUE
%     can be a scalar (assigns same value to all fields for all records), a
%     column vector (separate values for each record, but same value for
%     all fields in the group), a row vector (separate values for each
%     field in the group, but the same values across records), or an array
%     (separate values for all fields in all records).  Thus for group
%     field assignment, VALUE should be arranged so that columns delimit
%     values for each field in the group and rows delimit records 
%     (first row = first record, etc).  See the examples below to see how 
%     to replicate a set of values across several records.
%
%     CHANGEHEADER(DATA,FIELD1,VALUE1,...,FIELDN,VALUEN) allows changing
%     multiple fields in a single call.
%
%    Notes:
%     - Be careful with row vectors!
%     - Passing a nan, inf, -inf value will change a numeric field to its
%       undefined value.  Use 'nan', 'undef' or 'undefined' to do the same
%       for a character field.  This is useful for not having to remember
%       what the field's actual undefined value is.
%     - group fields:    t, kt, user, kuser, resp, dep, st, ev, nz, nzdttm,
%                         kname, {real group} utc, {real group} tai
%     - virtual fields:  nzmonth, nzcday, kzdttm, kzdate, kztime, z, ztai
%     - abs time fields: {real field} utc, {real field} tai
%
%    Header changes: Determined by input list.
%
%    Examples:
%     Some simple examples:
%      data=changeheader(data,'DELTA',getheader(data,'delta')*2);
%      data=changeheader(data,'STLA',lats,'STLO',lons)
%      data=changeheader(data,'KT0','sSKS');
%
%     Note that this will not split 'nan' to assign each letter to each
%     field in kuser because character arrays are first converted to cell
%     arrays.  So 'nan' is passed to each field in kuser and they all
%     become undefined (special behavior - see notes section):
%      data=changeheader(data,'kuser','nan')
%
%     Note that this case would bypass the cellstr conversion as the value
%     array is passed as a numeric array:
%      data=changeheader(data,'kuser',double('nan'))
%
%     To assign separate strings to fields in a character group field you
%     should pass a row vector cell array like this:
%      data=changeheader(data,'kuser',{'nan' 'blah' 'wow'})
%
%     But that replicates the same strings to all records in DATA, so to
%     assign separate values to all records and all fields (assuming there
%     are only 2 records):
%      data=changeheader(data,'kuser',{'nan' 'blah' 'wow'; 'another' ...
%           'set of' 'examples'});
%
%    See also:  LISTHEADER, GETHEADER, READHEADER, WRITEHEADER, GETLGC,
%               GETENUMID, GETENUMDESC, COMPAREHEADER

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
%        Apr. 23, 2009 - fix seizmocheck for octave, move usage up
%        Sep. 11, 2009 - fix splitting of values, changes behavior in
%                        multi-filetype cases and in group field cases,
%                        updated documentation to reflect this
%        Sep. 12, 2009 - added vgrp support
%        Sep. 15, 2009 - vf support, abs time support, doc update
%        Sep. 18, 2009 - 2nd pass at abs time support
%        Oct.  6, 2009 - dropped use of LOGICAL function
%        Oct. 16, 2009 - reftime code only used when necessary
%        Jan. 28, 2010 - eliminate extra struct checks
%        Jan. 29, 2010 - added VERSIONINFO cache support/hack
%        Feb. 16, 2010 - added informative errors for UTC/TAI fields
%        Apr. 13, 2010 - actually require fields are strings
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 13, 2010 at 12:00 GMT

% todo:

% throw error if unpaired fields
if (mod(nargin-1,2))
    error('seizmo:changeheader:badNargs','Unpaired Field/Value!');
end

% check data structure & get version info
[h,idx]=versioninfo(data);

% quick exit
if(nargin==1); return; end

% number of records
nrecs=numel(data);

% load SEIZMO info
global SEIZMO

% recursive section
% - break up into single filetype calls
nver=numel(h);
if(nver>1)
    % turn off struct checking
    oldseizmocheckstate=seizmocheck_state(false);
    oldversioninfocache=versioninfo_cache(true);
    
    try
        for i=1:nver
            % versioninfo cache hack
            SEIZMO.VERSIONINFO.H=h(i);
            SEIZMO.VERSIONINFO.IDX=ones(sum(idx==i),1);
            
            % need to parse varargin
            temp=cell(1,nargin-1);
            for j=1:2:nargin-2
                temp{j}=varargin{j};
            end
            for j=2:2:nargin-1
                sz=size(varargin{j});
                if(any(prod(sz)==[0 1]))
                    temp{j}=varargin{j};
                elseif(sz(1)==nrecs)
                    % dice up input only if a column vector
                    % - this prevents splitting up text unexpectedly
                    % - this is more stringent than previously
                    temp{j}=varargin{j}(idx==i,:);
                elseif(sz(1)==1)
                    temp{j}=varargin{j};
                else
                    error('seizmo:changeheader:invalidInputSize',...
                        'Value array for field %s incorrect size!',...
                        upper(varargin{j-1}));
                end
            end
            data(idx==i)=changeheader(data(idx==i),temp{:});
        end
        
        % toggle checking back
        seizmocheck_state(oldseizmocheckstate);
        versioninfo_cache(oldversioninfocache);
        
        % fix cache hack
        SEIZMO.VERSIONINFO.H=h;
        SEIZMO.VERSIONINFO.IDX=idx;
    catch
        % toggle checking back
        seizmocheck_state(oldseizmocheckstate);
        versioninfo_cache(oldversioninfocache);
        
        % fix cache hack
        SEIZMO.VERSIONINFO.H=h;
        SEIZMO.VERSIONINFO.IDX=idx;
        
        % rethrow error
        error(lasterror)
    end
    return;
end

% require all fields be strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:changeheader:badInput',...
        'FIELD(s) must be strings!');
end

% pull entire header
head=[data.head];

% preallocate ref
ref=[]; tai=[]; good=[];

% loop over field/value pairs
for i=1:2:(nargin-2)
    % force values into cell array
    if(isnumeric(varargin{i+1}) || islogical(varargin{i+1}))
        varargin{i+1}=num2cell(varargin{i+1});
    elseif(ischar(varargin{i+1}))
        varargin{i+1}=cellstr(varargin{i+1});
    end
    
    % force field name to be lowercase
    gf=strtrim(lower(varargin{i}));
    wf=getwords(gf);
    
    % check for group fields
    group=false; glen=1;
    if(isfield(h.vgrp,wf{1}))
        ngf=strtrim(strcat(h.vgrp.(wf{1}),{' '},joinwords(wf(2:end))));
        group=true; glen=numel(ngf);
    end
    
    % check & expand values
    if(isempty(varargin{i+1}))
        % skip changing if nothing to put there
        continue;
    elseif(isscalar(varargin{i+1}))
        % full expansion
        varargin{i+1}=varargin{i+1}(ones(nrecs,glen));
    elseif(isvector(varargin{i+1}))
        sz=size(varargin{i+1});
        
        % group or record expansion
        if(group)
            if(sz(1)==nrecs && nrecs~=1)
                % expand laterally (across group fields)
                varargin{i+1}=varargin{i+1}(:,ones(glen,1));
            elseif(sz(2)==glen && glen~=1)
                % expand vertically (across records)
                varargin{i+1}=varargin{i+1}(ones(nrecs,1),:);
            elseif(sz(2)==nrecs && nrecs~=1)
                % misoriented vector - transpose & expand laterally
                varargin{i+1}=varargin{i+1}(ones(glen,1),:).';
            elseif(sz(1)==glen && glen~=1)
                % misoriented vector - transpose & expand vertically
                varargin{i+1}=varargin{i+1}(:,ones(nrecs,1)).';
            else
                error('seizmo:changeheader:invalidInputSize',...
                    ['\nFiletype: %s, Version: %d\n'...
                    'Group value vector for field %s incorrect size!'],...
                    h.filetype,h.version,upper(varargin{i}));
            end
        elseif(prod(sz)==nrecs)
            % assure column vector
            varargin{i+1}=varargin{i+1}(:);
        else
            % not in a group, but is a vector of unusual length
            error('seizmo:changeheader:invalidInputSize',...
                ['\nFiletype: %s, Version: %d\n'...
                'Value vector for field %s incorrect size!'],...
                h.filetype,h.version,upper(varargin{i}));
        end
    % array
    else
        % no expansion, just check for correct orientation
        % - allow number of columns to exceed the number of
        %   fields in the group (to allow for variable group
        %   size across different filetypes)
        dim=size(varargin{i+1});
        if(dim(1)~=nrecs || dim(2)<glen)
            error('seizmo:changeheader:invalidInputSize',...
                ['\nFiletype: %s, Version: %d\n'...
                'Group value array for field %s not correct size'],...
                h.filetype,h.version,upper(varargin{i}))
        end
    end
    
    % group field loop
    for j=1:glen
        % modify field if group
        if(group); f=ngf{j};
        else f=gf;
        end
        
        % assign values to field
        [head,ref,tai,good]=a2h(head,h,f,varargin{i+1}(:,j),ref,tai,good);
    end
end

% put head back in
head=mat2cell(head,h.size,ones(1,nrecs));
[data.head]=deal(head{:});

end


function [head,ref,tai,good]=a2h(head,h,f,v,ref,tai,good)
%A2H    Subfunction that assigns values to header field

% virtual fields
if(isfield(h.vf,f))
    switch h.vf.(f).type
        case {'enum' 'lgc' 'int' 'real'}
            head=h.vf.(f).ch(h,head,cell2mat(v));
        case {'char' 'abs'}
            head=h.vf.(f).ch(h,head,v);
    end
    return;
end

% special treatment for enums
for m=1:numel(h.enum)
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
            if(numel(unique(v))==1)
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
                            h.filetype,h.version,upper(f));
                    end
                end
            else % not all are the same (slow b/c looped)
                for n=1:numel(v)
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
                                h.filetype,h.version,upper(f));
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
for m=1:numel(h.lgc)
    if(isfield(h.lgc(m).pos,f))
        if(iscellstr(v))
            % logic words (unknown => undefined here)
            lgctrue=strncmpi(v,'t',1);
            lgcfalse=strncmpi(v,'f',1);
            undef=strncmpi(v,'u',1) | strcmpi(v,'nan');
            head(h.lgc(m).pos.(f),lgctrue)=h.true;
            head(h.lgc(m).pos.(f),lgcfalse)=h.false;
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
for n=1:numel(h.stype)
    for m=1:numel(h.(h.stype{n}))
        if(isfield(h.(h.stype{n})(m).pos,f))
            p=h.(h.stype{n})(m).pos.(f);
            o=p(2)-p(1)+1;
            if(iscellstr(v))
                v(strcmpi(v,'undef'))={h.undef.stype};
                v(strcmpi(v,'undefined'))={h.undef.stype};
                v(strcmpi(v,'nan'))={h.undef.stype};
            else
                % assume numeric matrix (will be changed to cellstr)
                % (no check - breaks if mixed num/char cells)
                v=cell2mat(v);
                undeftmp1=isnan(v);
                undeftmp2=isinf(v);
                v(undeftmp1 | undeftmp2)=0;
                v=cellstr(char(v));
                v(sum(undeftmp1,2)~=0)={h.undef.stype};
                v(sum(undeftmp2,2)~=0)={h.undef.stype};
            end
            v=strnlen(v,o);
            head(p(1):p(2),:)=char(v).';
            return;
        end
    end
end

% remainder of numeric types
for n=1:numel(h.ntype)
    for m=1:numel(h.(h.ntype{n}))
        if(isfield(h.(h.ntype{n})(m).pos,f))
            v=cell2mat(v);
            v(isnan(v))=h.undef.ntype;
            v(isinf(v))=h.undef.ntype;
            head(h.(h.ntype{n})(m).pos.(f),:)=v;
            return;
        elseif(strcmp(h.ntype{n},'real'))
            % absolute time fields section
            wf=getwords(f);
            if(any(strcmpi(joinwords(wf(2:end)),...
                    {'utc' 'tai' '6utc' '6tai'})))
                if(isfield(h.real(m).pos,wf{1}))
                    % get reftimes
                    if(isempty(ref))
                        [ref,good]=vf_gh_z(h,head); good=good';
                        tai=ref; tai(good,:)=utc2tai(ref(good,:));
                    end
                    
                    % set all to undef
                    head(h.real(m).pos.(wf{1}),:)=h.undef.ntype;
                    
                    % who's (un)defined
                    nv=numel(v);
                    good5=false(nv,1);
                    good6=good5;
                    for i=1:nv
                        sv=size(v{i});
                        % - need errors to remind myself how this works
                        if(~isreal(v{i}))
                            error('seizmo:changeheader:badValue',...
                                ['Field: %s\nVALUE for a UTC/TAI ' ...
                                'field must be real!'],upper(f));
                        elseif(~isequal(v{i}(1:end-1),fix(v{i}(1:end-1))))
                            error('seizmo:changeheader:badValue',...
                                ['Field: %s\nVALUE for a UTC/TAI ' ...
                                'field must be whole numbers except ' ...
                                'for the seconds element!'],upper(f));
                        end
                        % leave as undef if undef, nan, or inf
                        if(~any(v{i}==h.undef.ntype) ...
                                && ~any(isnan(v{i}) | isinf(v{i})))
                            % require 1x5/6
                            if(isequal(sv,[1 5]))
                                good5(i)=true;
                            elseif(isequal(sv,[1 6]))
                                good6(i)=true;
                            else
                                error('seizmo:changeheader:badValue',...
                                    ['Field: %s\nVALUE for a UTC/TAI ' ...
                                    'field must be celled 1x5 or 1x6 ' ...
                                    'arrays!'],upper(f));
                            end
                        end
                    end
                    good5=good5 & good;
                    good6=good6 & good;
                    
                    % skip empty
                    if(any(good5))
                        % fix times
                        switch wf{2}
                            case {'utc' '6utc'}
                                head(h.real(m).pos.(wf{1}),good5)=...
                                    timediff(tai(good5,:),...
                                    utc2tai(cell2mat(v(good5,1))));
                            case {'tai' '6tai'}
                                head(h.real(m).pos.(wf{1}),good5)=...
                                    timediff(tai(good5,:),...
                                    cell2mat(v(good5,1)));
                        end
                    end
                    if(any(good6))
                        % fix times
                        switch wf{2}
                            case {'utc' '6utc'}
                                head(h.real(m).pos.(wf{1}),good6)=...
                                    timediff(tai(good6,:),...
                                    utc2tai(cell2mat(v(good6,1))));
                            case {'tai' '6tai'}
                                head(h.real(m).pos.(wf{1}),good6)=...
                                    timediff(tai(good6,:),...
                                    cell2mat(v(good6,1)));
                        end
                    end
                    return;
                end
            end
        end
    end
end
            
% field not found
warning('seizmo:changeheader:fieldInvalid',...
    'Filetype: %s, Version: %d\nInvalid field: %s !',...
    h.filetype,h.version,upper(f));

end
