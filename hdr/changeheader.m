function [data]=changeheader(data,varargin)
%CHANGEHEADER    Change SEIZMO data header values
%
%    Usage: data=changeheader(data,'field',values)
%           data=changeheader(data,'field abstime',values)
%           data=changeheader(data,'groupfield',values)
%           data=changeheader(data,'virtualfield',values)
%           data=changeheader(data,'field1',values1,...,'fieldN',valuesN)
%
%    Description:
%     DATA=CHANGEHEADER(DATA,FIELD,VALUES) changes the header field FIELD
%     to the value(s) in VALUES for each record in SEIZMO structure DATA
%     and returns the updated structure.  FIELD must be a string 
%     corresponding to a valid header field (you can use the LISTHEADER
%     command to get an idea of what header fields there are).  VALUES
%     may be a scalar (assigns the same value to all records), a column
%     vector of length equal to the number of records (assigns a separate
%     value to each record), or a char array (which is converted to a
%     cellstr array internally, and thus becomes a column vector).
%
%     DATA=CHANGEHEADER(DATA,'FIELD ABSTIME',VALUES) sets the field FIELD
%     to the absolute time VALUES using the standard ABSTIME.  ABSTIME must
%     be a string of the following choices:
%         UTC  - input time is in UTC as [yr doy hr min sec]
%         UTC6 - input time is in UTC as [yr mon cday hr min sec]
%         TAI  - input time is in TAI as [yr doy hr min sec]
%         TAI6 - input time is in TAI as [yr mon cday hr min sec]
%     Please note that while CHANGEHEADER may generally parse VALUES
%     correctly, passing VALUES as a cell array with each time in a
%     separate cell works best.  This is consistent with GETHEADER output.
%     Also you do not need to use UTC6 or TAI6 here (UTC or TAI can handle
%     both Nx5 & Nx6 input types) but it is available to be consistent with
%     the other header functions.
%
%     DATA=CHANGEHEADER(DATA,GROUPFIELD,VALUES) allows setting a predefined
%     group of header fields simultaneously.  The groups fields are:
%         T, KT, USER, KUSER, RESP, DEP, ST, EV, NZ, NZDTTM,
%         KNAME, CMP, DELAZ, REAL, INT, ENUM, LGC, CHAR
%     The members of the group field can be seen by using the VGRP command
%     on a group field.  VALUES can be a scalar (assigns same value to all
%     fields for all records), a column vector (separate values for each
%     record, but same value for all fields in the group), a row vector
%     (separate values for each field in the group, but the same values
%     across all records), or an array (separate values for all fields in
%     all records).  Thus for group field assignment, VALUE should be
%     arranged so that columns delimit values for each field in the group
%     and rows delimit records (first row = first record, etc).
%
%     DATA=CHANGEHEADER(DATA,VIRTUALFIELD,VALUES) alters a predefined
%     virtual field which is then alters 1 or more true header fields.
%     Available virtual fields:
%         NZMONTH  - calendar month of the reference time (from NZ* fields)
%         NZCDAY   - calendar day of the reference time (from NZ* fields)
%         KZDTTM   - formatted string of the reference date & time
%         KZDATE   - formatted string of the reference date
%         KZTIME   - formatted string of the reference time
%         Z        - reference time in UTC as [yr doy hr min sec]
%         Z6       - reference time in UTC as [yr mon cday hr min sec]
%         ZTAI     - reference time in TAI as [yr doy hr min sec]
%         ZTAI6    - reference time in TAI as [yr mon cday hr min sec]
%         GCP      - greater circle path azimuth (BAZ+180deg)
%         NCMP     - number of components (=1 for most cases)
%
%     DATA=CHANGEHEADER(DATA,FIELD1,VALUE1,...,FIELDN,VALUEN) changes
%     multiple fields in a single call.  FIELD may be any of the above
%     types.
%
%    Notes:
%     - Be careful with using row vectors as values!  Row vectors outside
%       of group field usage fails if there are multiple filetypes in use.
%       This isn't the usual case so you can get away with it most of the
%       time.  You have been warned!
%     - Passing nan, inf, -inf values will change a numeric field to its
%       undefined value.  'nan', 'undef' or 'undefined' will do the same
%       for a character field.  This is useful for not having to remember
%       what the field's actual undefined value is.
%     - Enumerated/Logical fields do not need modifiers as they are auto
%       detected in CHANGEHEADER, but Absolute Times do need 'utc' or 'tai'
%       after the field as a separate word (eg. 'b utc').
%     - Use SYNCHRONIZE or TIMESHIFT to update the reference time while
%       preserving the absolute timing of all timing fields.  CHANGEHEADER
%       will not do this!  You have been warned!
%
%    Header changes: Determined by input list.
%
%    Examples:
%     % Some simple examples:
%     data=changeheader(data,'user0',500);
%     data=changeheader(data,'iztype','iunkn');
%     data=changeheader(data,'t0',1650,'Kt0','sSKS');
%     data=changeheader(data,'StLA',lats,'STLo',lons)
%     data=changeheader(data,'DeLTA',getheader(data,'delta')*2);
%
%     % Clear picks:
%     data=changeheader(data,'t',nan,'kt','nan');
%
%     % Set begin time to right now:
%     data=changeheader(data,'b utc',datevec(now))
%
%     % Wrap absolute times in cells to assure proper parsing:
%     data=changeheader(data,'b utc',{[2012 12 21 0 0 0]});
%
%     % Change the month of the records by one:
%     data=changeheader(data,'nzmonth',getheader(data,'nzmonth')+1);
%
%    See also: LISTHEADER, GETHEADER, READHEADER, WRITEHEADER, QUERYHEADER,
%              COMPAREHEADER, VGRP

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
%        Jan. 28, 2012 - drop strnlen usage
%        Jan. 30, 2012 - doc update, 6utc/6tai changed to utc6/tai6,
%                        some support for abstimes without modifier, allow
%                        abstimes input to be uncelled
%        Mar.  1, 2012 - bugfix: forgot to define solo, drop bad cell2mat
%        Aug. 30, 2012 - big doc update for clarity
%        June 26, 2014 - bugfix: abstime correctly updates in case abstime
%                        then reftime then abstime modified (note this is
%                        probably not what you want to do as changing the
%                        reftime with changeheader doesn't update the other
%                        time fields, thus changing their abstime - so use
%                        synchronize or timeshift to properly update times)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 26, 2014 at 12:00 GMT

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
                        'VALUE array for field %s incorrect size!',...
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
        error(lasterror);
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
            elseif(isequal(sz,[1 5]) || isequal(sz,[1 6]))
                % probably uncelled AbsTime -- cell it properly!
                % and expand fully...
                varargin{i+1}={cell2mat(varargin{i+1})};
                varargin{i+1}=varargin{i+1}(ones(nrecs,glen));
            elseif(isequal(sz,[1 5*glen]))
                % probably uncelled AbsTime -- cell it properly!
                % and expand vertically...
                varargin{i+1}=mat2cell(cell2mat(varargin{i+1}),...
                    ones(sz(1),1),5*ones(1,glen));
                varargin{i+1}=varargin{i+1}(ones(nrecs,1),:);
            elseif(isequal(sz,[1 6*glen]))
                % probably uncelled AbsTime -- cell it properly!
                % and expand vertically...
                varargin{i+1}=mat2cell(cell2mat(varargin{i+1}),...
                    ones(sz(1),1),6*ones(1,glen));
                varargin{i+1}=varargin{i+1}(ones(nrecs,1),:);
            else
                error('seizmo:changeheader:invalidInputSize',...
                    ['\nFiletype: %s, Version: %d\n'...
                    'VALUE vector for field %s incorrect size!'],...
                    h.filetype,h.version,upper(varargin{i}));
            end
        elseif(prod(sz)==nrecs)
            % row vector to column vector
            varargin{i+1}=varargin{i+1}(:);
        elseif(isequal(sz,[1 5]) || isequal(sz,[1 6]))
            % probably uncelled AbsTime -- cell it properly!
            % and expand fully...
            varargin{i+1}={cell2mat(varargin{i+1})};
            varargin{i+1}=varargin{i+1}(ones(nrecs,glen));
        else
            % not in a group, but is a vector of unusual length
            error('seizmo:changeheader:invalidInputSize',...
                ['\nFiletype: %s, Version: %d\n'...
                'VALUE vector for field %s incorrect size!'],...
                h.filetype,h.version,upper(varargin{i}));
        end
    % array
    else
        % check for correct orientation
        sz=size(varargin{i+1});
        if(sz(1)~=nrecs)
            error('seizmo:changeheader:invalidInputSize',...
                ['\nFiletype: %s, Version: %d\n'...
                'VALUE array for field %s not correct size'],...
                h.filetype,h.version,upper(varargin{i}))
        elseif(sz(2)~=glen)
            if(any(sz(2)==[5 6]))
                % probably uncelled AbsTime -- cell it properly!
                % and expand horizontally...
                varargin{i+1}=mat2cell(cell2mat(varargin{i+1}),...
                    ones(sz(1),1));
                varargin{i+1}=varargin{i+1}(:,ones(1,glen));
            elseif(sz(2)==5*glen)
                % probably uncelled AbsTime -- cell it properly!
                varargin{i+1}=mat2cell(cell2mat(varargin{i+1}),...
                    ones(sz(1),1),5*ones(1,glen));
            elseif(sz(2)==6*glen)
                % probably uncelled AbsTime -- cell it properly!
                varargin{i+1}=mat2cell(cell2mat(varargin{i+1}),...
                    ones(sz(1),1),6*ones(1,glen));
            else
                error('seizmo:changeheader:invalidInputSize',...
                    ['\nFiletype: %s, Version: %d\n'...
                    'VALUE array for field %s not correct size'],...
                    h.filetype,h.version,upper(varargin{i}))
            end
        end
    end
    
    % group field loop
    for j=1:glen
        % modify field if group
        if(group); f=ngf{j};
        else f=gf;
        end
        
        % assign values to field
        head=a2h(head,h,f,varargin{i+1}(:,j));
    end
end

% put head back in
head=mat2cell(head,h.size,ones(1,nrecs));
[data.head]=deal(head{:});

end


function [head]=a2h(head,h,f,v)
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

% multiwords
wf=getwords(f); f=wf{1};
if(numel(wf)==1); solo=true; else solo=false; end

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
            v=char(v).';
            szv=size(v,1);
            if(o>szv); v(szv+1:o,:)=32; end
            head(p(1):p(2),:)=v(1:o,:);
            return;
        end
    end
end

% remainder of numeric types
for n=1:numel(h.ntype)
    for m=1:numel(h.(h.ntype{n}))
        if(isfield(h.(h.ntype{n})(m).pos,f))
            v=cell2mat(v); sv=size(v);
            % check for abstime
            if(any(sv(2)==[5 6]) && strcmp(h.ntype{n},'real'))
                % warn about no modifier
                if(solo)
                    warning('seizmo:changeheader:noModifier',...
                        ['Field: %s\n' ...
                        'No modifier given for AbsTime input! ' ...
                        'Assuming UTC!'],upper(f));
                    wf{2}='utc';
                end
                
                % get reftimes
                [tai,good]=vf_gh_ztai(h,head);
                good=good';
                
                % set all to undef
                head(h.real(m).pos.(f),:)=h.undef.ntype;
                
                % need errors to remind myself how this works
                if(~isreal(v))
                    error('seizmo:changeheader:badValue',...
                        ['Field: %s\nVALUE for a UTC/TAI ' ...
                        'field must be real!'],upper(f));
                elseif(~isequal(v(:,1:end-1),fix(v(:,1:end-1))))
                    error('seizmo:changeheader:badValue',...
                        ['Field: %s\nVALUE for a UTC/TAI ' ...
                        'field must be whole numbers except ' ...
                        'for the seconds element!'],upper(f));
                end
                
                % who's (un)defined
                good=good ...
                    & sum(v==h.undef.ntype | isnan(v) | isinf(v),2)==0;
                
                % skip empty
                if(any(good))
                    % fix times
                    switch wf{2}
                        case {'utc' 'utc6'}
                            head(h.real(m).pos.(f),good)=...
                                timediff(tai(good,:),utc2tai(v(good,:)));
                        case {'tai' 'tai6'}
                            head(h.real(m).pos.(f),good)=...
                                timediff(tai(good,:),v(good,:));
                    end
                end
                return;
            elseif(sv(2)==1) % normal case
                v(isnan(v))=h.undef.ntype;
                v(isinf(v))=h.undef.ntype;
                head(h.(h.ntype{n})(m).pos.(f),:)=v;
                return;
            else
                error('seizmo:changeheader:badInput',...
                    ['Filetype: %s, Version: %d\nField: %s\n!' ...
                    'VALUE array is incorrect size!'],...
                    h.filetype,h.version,upper(f));
            end
        end
    end
end
            
% field not found
warning('seizmo:changeheader:fieldInvalid',...
    'Filetype: %s, Version: %d\nInvalid field: %s !',...
    h.filetype,h.version,upper(f));

end
