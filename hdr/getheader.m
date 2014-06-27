function [varargout]=getheader(data,varargin)
%GETHEADER    Get SEIZMO data header values
%
%    Usage:    values=getheader(data,'field')
%              values=getheader(data,'field abstime')
%              values=getheader(data,'enumfield id')
%              values=getheader(data,'enumfield desc')
%              values=getheader(data,'logicfield lgc')
%              values=getheader(data,'groupfield')
%              values=getheader(data,'virtualfield')
%              [values1,...,valuesN]=getheader(data,'field1',...,'fieldN')
%              headers=getheader(data)
%
%    Description:
%     VALUES=GETHEADER(DATA,FIELD) returns the specified header field
%     values for all records stored in the SEIZMO data structure DATA.
%     FIELD must be a string corresponding to a valid header field and is
%     case insensitive.  Values are returned as a numeric array or
%     cellstring array (for character fields) oriented as a column vector
%     where each row corresponds to an individual record.
%
%     VALUES=GETHEADER(DATA,'FIELD ABSTIME') returns the field FIELD in
%     absolute time using the standard ABSTIME.  ABSTIME must be a string
%     as one of the following choices:
%         UTC  - return values in UTC as [yr doy hr min sec]
%         UTC6 - return values in UTC as [yr mon cday hr min sec]
%         TAI  - return values in TAI as [yr doy hr min sec]
%         TAI6 - return values in TAI as [yr mon cday hr min sec]
%     Please note that the output is a cell array with each cell
%     corresponding to a absolute time for each record/field.
%
%     VALUES=GETHEADER(DATA,'ENUMFIELD ID') retrieves the associated string
%     id of the enum field ENUMFIELD.  This is more useful than the default
%     output (the enumerated table index).  An example id output would be
%     'idisp' for the field 'idep' which is far easier to read than the
%     equivalent raw output of 6.  Output is a cell-string array.  The
%     string 'NaN' is returned for values that are outside the enum table.
%
%     VALUES=GETHEADER(DATA,'ENUMFIELD DESC') retrieves the associated
%     string description of the enum field ENUMFIELD.  This is useful for
%     descriptive output.  An example description output woud be
%     'Displacement (nm)' for the field 'idep'.  Output is a cell-string
%     array.  The string 'NaN' is returned for values that are outside the
%     enum table.
%
%     VALUES=GETHEADER(DATA,'LOGICFIELD LGC') retrieves logical strings
%     associated with the logical field LOGICFIELD.  VALUES is a cell
%     string array composed of 'true', 'false' and/or 'nan'.  This is
%     easier to read and helps to avoid confusion compared to the default
%     output of 0, 1 or nan.
%
%     VALUES=GETHEADER(DATA,'GROUPFIELD') returns the values of a
%     predefined group of header fields as a single array simultaneously.
%     Available group fields are:
%         T, KT, USER, KUSER, RESP, DEP, ST, EV, NZ, NZDTTM,
%         KNAME, CMP, DELAZ, REAL, INT, ENUM, LGC, CHAR
%     The members of the group field can be seen by using the VGRP command
%     on a group field.  VALUES is a numeric or cellstring array where each
%     row corresponds to a seperate record in DATA and each column
%     corresponds to a separate field.  For example, the group field 't'
%     returns a numeric array with 10 columns (for fields t0 to t9) and as
%     many rows as there are records in DATA.  Group field 'kuser' returns
%     a cellstring array with 3 columns (for kuser0 thru kuser2) and as
%     many rows as there are records in DATA.  Output modifiers described
%     above (abstime, id, desc, lgc) can be combined with group fields as
%     appropriate.
%
%     VALUES=GETHEADER(DATA,'VIRTUALFIELD') retrieves predefined virtual
%     field values defined by 1 or more true header fields.  Available
%     virtual fields:
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
%     [VALUES1,...,VALUESN]=GETHEADER(DATA,FIELD1,...,FIELDN) returns
%     multiple fields in a single call.  FIELD may be any of the above
%     types.
%
%     HEADERS=GETHEADER(DATA) returns the entire header for all records in
%     SEIZMO struct DATA as a single numeric array (raw values).  Rows in
%     the output array correspond to the header values of individual
%     records.  The order of the fields follows that of how they are stored
%     in the SEIZMO struct (see SEIZMO's function SEIZMODEF for details).
%
%    Notes:
%     - Undefined header field values (-12345 in SAC) are returned as nan
%       or 'NaN' for ease of use.
%
%    Examples:
%     % Extract the sample rates of records:
%     dt=getheader(data,'DeLtA');
%
%     % Get the station lat and lon for records:
%     [stla,stlo]=getheader(data,'stla','STLO');
%
%     % Retrieve all t series values as one array:
%     times=getheader(data,'t');
%
%     % Return the record datatype id:
%     iftype=getheader(data,'iftype id');
%
%     % Find records that are evenly spaced:
%     even=find(strcmp(getheader(data,'leven lgc'),'true'));
%
%     % Get the start time of records in UTC:
%     butc=cell2mat(getheader(data,'b utc'));
%
%    See also:  CHANGEHEADER, LISTHEADER, READHEADER, WRITEHEADER,
%               FINDPICKS, COMPAREHEADER, QUERYHEADER, VGRP

%     Version History:
%        Oct. 29, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 28, 2008 - new sachp support
%        Feb. 18, 2008 - rewrite - parses by version first
%        Feb. 23, 2008 - minor doc update
%        Feb. 28, 2008 - code cleanup
%        Feb. 29, 2008 - minor code cleaning
%        Mar.  4, 2008 - doc update
%        June 12, 2008 - doc update, full header dump fixes
%        Oct. 17, 2008 - added VINFO support, supports new struct layout
%        Nov. 16, 2008 - update for new name schema (now GETHEADER)
%        Apr. 23, 2009 - fix seizmocheck for octave, move usage up
%        Sep. 12, 2009 - added vgrp support
%        Sep. 15, 2009 - vf support, abs time support, doc update
%        Sep. 18, 2009 - 2nd pass at abs time support
%        Oct.  6, 2009 - dropped use of LOGICAL function
%        Oct. 16, 2009 - reftime code only used when necessary
%        Jan. 28, 2010 - eliminate extra struct checks
%        Jan. 29, 2010 - added VERSIONINFO cache support/hack
%        Apr. 13, 2010 - actually require fields are strings
%        Aug. 21, 2010 - doc update, NaN output masks undef (-12345)
%        Apr. 19, 2011 - bugfix for z6 & friends when 2+ header types
%        Jan. 28, 2012 - code refactoring in ph (for id/desc/lgc support)
%        Jan. 30, 2012 - doc update, 6utc/6tai changed to utc6/tai6,
%                        id/desc/lgc support
%        Mar. 15, 2012 - minor doc touch
%        Aug. 30, 2012 - big doc update for clarity
%        June 26, 2014 - bugfix: undefined/bad enum could cause abstime of
%                        following fields to be returned as undefined,
%                        bugfix: undefined abstime fields would cause
%                        abstime of following fields to be undefined
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 26, 2014 at 12:00 GMT

% todo:

% require at least one input
if(nargin<1)
    error('seizmo:getheader:notEnoughInputs',...
        'Not enough input arguments.')
end

% check data structure & grab header setup
[h,idx]=versioninfo(data);

% number of records
nrecs=numel(data);

% load SEIZMO info
global SEIZMO

% recursive section
% - break up into single filetype calls
nver=numel(h);
sfill=repmat({'NaN'},nrecs,1);
nfill=nan(nrecs,1);
if(nver>1)
    % turn off struct checking
    oldseizmocheckstate=seizmocheck_state(false);
    oldversioninfocache=versioninfo_cache(true);
    
    try
        nout=max([1 nargin-1]);
        for i=1:nver
            % versioninfo cache hack
            SEIZMO.VERSIONINFO.H=h(i);
            SEIZMO.VERSIONINFO.IDX=ones(sum(idx==i),1);
            
            varargoutn=cell(1,nout);
            [varargoutn{:}]=getheader(data(idx==i),varargin{:});
            
            % assign to varargout
            for j=1:nout
                % preallocate by type
                if(i==1)
                    if(iscell(varargoutn{j}))
                        varargout{j}=sfill; type=1;
                    else
                        varargout{j}=nfill; type=0;
                    end
                end
                % expand
                in=size(varargoutn{j},2);
                out=size(varargout{j},2);
                if(in>out)
                    if(type)
                        varargout{j}(:,out+1:in)=sfill(:,ones(1,in-out));
                    else
                        varargout{j}(:,out+1:in)=nfill(:,ones(1,in-out));
                    end
                end
                % assign
                varargout{j}(idx==i,1:in)=varargoutn{j};
            end
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

% push all headers into a single matrix
head=[data.head];

% push out entire header
if(nargin==1); varargout{1}=head.'; return; end

% require all fields be strings
if(~iscellstr(varargin))
    error('seizmo:getheader:badInput',...
        'FIELD(s) must be strings!');
end

% preallocate ref
ref=[]; ref6=[]; good=[];

% loop over fields
for i=1:nargin-1
    % force field to be lowercase
    gf=strtrim(lower(varargin{i}));
    wf=getwords(gf);
    
    % check for group fields
    group=false; glen=1;
    if(isfield(h.vgrp,wf{1}))
        ngf=strtrim(strcat(h.vgrp.(wf{1}),{' '},joinwords(wf(2:end))));
        group=true; glen=numel(ngf);
    end
    
    % group field loop
    for j=1:glen
        % modify field if group
        if(group); f=ngf{j};
        else f=gf;
        end
        
        % pull header values
        [val,type,ref,ref6,good]=ph(head,h,f,ref,ref6,good);
        
        % preallocate (did not know the type until now)
        if(j==1)
            if(type)
                varargout{i}=repmat({'NaN'},nrecs,glen);
            else
                varargout{i}=nan(nrecs,glen);
            end
        end
        
        % assign to varargout
        % - note this forces columnar output
        varargout{i}(:,j)=val;
    end
    
    % convert undefined values to nans
    if(type)
        varargout{i}(strcmpi(varargout{i},h.undef.stype))={'NaN'};
    else
        varargout{i}(varargout{i}==h.undef.ntype ...
            | isinf(varargout{i}))=nan;
    end
end

end

function [head,type,ref,ref6,good]=ph(head,h,f,ref,ref6,good)
%PH    Pull header values

% virtual fields
if(isfield(h.vf,f))
    switch h.vf.(f).type
        case {'enum' 'lgc' 'int' 'real'}
            head=h.vf.(f).gh(h,head);
            type=0;
        case {'char' 'abs'}
            head=h.vf.(f).gh(h,head);
            type=1;
    end
    return;
end

% output by type
for n=1:numel(h.stype)
    for m=1:numel(h.(h.stype{n}))
        if(isfield(h.(h.stype{n})(m).pos,f))
            p=h.(h.stype{n})(m).pos.(f);
            head=cellstr(char(head(p(1):p(2),:).'));
            type=1; return;
        end
    end
end
wf=getwords(f); solo=numel(wf)==1;
if(solo)
    for n=1:numel(h.ntype)
        for m=1:numel(h.(h.ntype{n}))
            if(isfield(h.(h.ntype{n})(m).pos,f))
                head=head(h.(h.ntype{n})(m).pos.(f),:).';
                type=0; return;
            end
        end
    end
else % multiword
    switch wf{2}
        case {'utc' 'tai' 'utc6' 'tai6'}
            for m=1:numel(h.real)
                if(isfield(h.real(m).pos,wf{1}))
                    % get reftimes
                    if(isempty(ref))
                        [ref,good]=vf_gh_z(h,head);
                        good=good';
                        ref6=ref(:,[1:2 2:5]);
                        ref6(good,1:3)=doy2cal(ref6(good,1:2));
                    end
                    
                    % deal with doy/cal split
                    if(wf{2}(end)=='6'); z=6; else z=5; end
                    
                    % get header values in a workable form
                    nrecs=size(head,2);
                    value=zeros(nrecs,z);
                    value(:,z)=head(h.real(m).pos.(wf{1}),:).';
                    
                    % default output to nan
                    head=nan(size(value,1),z);
                    
                    % who's (un)defined
                    vgood=good & value(:,z)~=h.undef.ntype ...
                        & ~isnan(value(:,z)) & ~isinf(value(:,z));
                    
                    % skip if empty
                    if(any(vgood))
                        switch wf{2}
                            case 'utc'
                                head(vgood,:)=fixtimes(ref(vgood,:)...
                                    +value(vgood,:),'utc');
                            case 'tai'
                                head(vgood,:)=fixtimes(utc2tai(...
                                    ref(vgood,:))+value(vgood,:));
                            case 'utc6'
                                head(vgood,:)=fixtimes(ref6(vgood,:)...
                                    +value(vgood,:),'utc');
                            case 'tai6'
                                head(vgood,:)=fixtimes(utc2tai(...
                                    ref6(vgood,:))+value(vgood,:));
                        end
                    end
                    head=mat2cell(head,ones(nrecs,1)); type=1;
                    return;
                end
            end
        case {'id' 'desc'}
            for m=1:numel(h.enum)
                if(isfield(h.enum(m).pos,wf{1}))
                    enum=head(h.enum(m).pos.(wf{1}),:).';
                    head=cell(size(enum));
                    head(:)={'NaN'};
                    idgood=fix(enum)==enum ...
                        & enum>=h.enum(1).minval & enum<=h.enum(1).maxval;
                    head(idgood)=h.enum(1).(wf{2})(enum(idgood)+1);
                    type=1; return;
                end
            end
        case 'lgc'
            for m=1:numel(h.lgc)
                if(isfield(h.lgc(m).pos,wf{1}))
                    lgc=head(h.lgc(m).pos.(wf{1}),:).';
                    head=cell(size(lgc));
                    head(:)={'NaN'};
                    head(lgc==h.true)={'true'};
                    head(lgc==h.false)={'false'};
                    type=1; return;
                end
            end
    end
end

% field not found
warning('seizmo:getheader:fieldInvalid',...
    'Filetype: %s, Version: %d\nInvalid field: %s',...
    h.filetype,h.version,f);
head=nan;
type=0;

end
