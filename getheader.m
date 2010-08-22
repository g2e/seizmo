function [varargout]=getheader(data,varargin)
%GETHEADER    Get SEIZMO data header values
%
%    Usage:  headers=getheader(data)
%            values=getheader(data,'field1')
%            [values1,...,valuesN]=getheader(data,'field1',...,'fieldN')
%
%    Description: GETHEADER(DATA) returns ALL header values for all records
%     in SEIZMO struct DATA as a single numeric array.  Rows in the output
%     array correspond to the header values of individual records.  The
%     order of the fields follows that of how they are stored in the SEIZMO
%     struct (see SEIZMO's function SEISMODEF for details).
%
%     GETHEADER(DATA,FIELD) returns the specified header field FIELD's
%     values for all records stored in the SEIZMO data structure DATA.
%     FIELD must be a string corresponding to a valid header field or a
%     valid group field (ie. t,kt,resp,user,kuser).  Values are returned in
%     numeric arrays or cellstring arrays oriented such that each column
%     corresponds to an individual header field and each row to an
%     individual record.  For example, the group field 't' would return a
%     numeric array with 10 columns (for t0 through t9) and as many rows as
%     there are records in DATA.  Group field 'kuser' would return a
%     cellstring array with 3 columns (for kuser0 thru kuser2) and as many
%     rows as there are records in DATA.  FIELD is case-insensitive.
%     
%     [VALUES1,...,VALUESN]=GETHEADER(DATA,FIELD1,...,FIELDN) returns the
%     values for each specified field or group field.
%
%    Notes:
%     - Enumerated fields return the value actually stored, an integer used
%       to look up the enum string id and description in a table.  To 
%       retrieve the associated string id or description see the functions 
%       GETENUMID or GETENUMDESC.
%     - Logical fields return the value actually stored, not a logical.  To
%       get a more useful value see GETLGC.
%     - group fields:    t, kt, user, kuser, resp, dep, st, ev, nz, nzdttm,
%                         kname, {real group} utc, {real group} tai
%     - virtual fields:  nzmonth, nzcday, kzdttm, kzdate, kztime, z, ztai
%     - abs time fields: {real field} utc, {real field} tai
%
%    Examples:
%     % Retrieve all t series values as one array:
%     times=getheader(data,'t')
%
%     % Extract the sample rates of records:
%     dt=getheader(data,'DeLtA')
%
%     % Get the station lat and lon for records:
%     [stla,stlo]=getheader(data,'stla','STLO')
%
%     % Enumerated fields return the table lookup index,
%     % which is the value stored in the header:
%     getheader(data,'iftype')
%
%    See also:  CHANGEHEADER, LISTHEADER, READHEADER, WRITEHEADER, GETLGC,
%               GETENUMID, GETENUMDESC, GETNCMP, GETARRIVAL, COMPAREHEADER

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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 21, 2010 at 12:00 GMT

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
                    if(iscellstr(varargoutn{j}))
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
        error(lasterror)
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
for n=1:numel(h.ntype)
    for m=1:numel(h.(h.ntype{n}))
        if(isfield(h.(h.ntype{n})(m).pos,f))
            head=head(h.(h.ntype{n})(m).pos.(f),:).';
            type=0; return;
        elseif(strcmpi(h.ntype{n},'real'))
            % absolute time fields section
            wf=getwords(f);
            if(numel(wf)>1 && isfield(h.real(m).pos,wf{1}))
                if(any(strcmpi(joinwords(wf(2:end)),{'utc' 'tai'})))
                    % get reftimes
                    if(isempty(ref))
                        [ref,good]=vf_gh_z(h,head); good=good';
                        ref6=ref(:,[1:2 2:5]);
                        ref6(good,1:3)=doy2cal(ref6(good,1:2));
                    end
                    
                    % get header values in a workable form
                    nrecs=size(head,2);
                    value=zeros(nrecs,5);
                    value(:,5)=head(h.(h.ntype{n})(m).pos.(wf{1}),:).';
                    
                    % default output to nan
                    head=nan(size(value,1),5);
                    
                    % who's (un)defined
                    good=good & value(:,5)~=h.undef.ntype ...
                        & ~isnan(value(:,5)) & ~isinf(value(:,5));
                    
                    % skip if empty
                    if(any(good))
                        switch wf{2}
                            case 'utc'
                                head(good,:)=fixtimes(ref(good,:)...
                                    +value(good,:),'utc');
                            case 'tai'
                                head(good,:)=fixtimes(utc2tai(...
                                    ref(good,:))+value(good,:));
                        end
                    end
                    head=mat2cell(head,ones(nrecs,1)); type=1;
                    return;
                elseif(any(strcmpi(joinwords(wf(2:end)),{'6utc' '6tai'})))
                    % get reftimes
                    if(isempty(ref6))
                        [ref,good]=vf_gh_z(h,head); good=good';
                        ref6=ref(:,[1:2 2:5]);
                        ref6(good,1:3)=doy2cal(ref6(good,1:2));
                    end
                    
                    % get header values in a workable form
                    nrecs=size(head,2);
                    value=zeros(nrecs,6);
                    value(:,6)=head(h.(h.ntype{n})(m).pos.(wf{1}),:).';
                    
                    % default output to nan
                    head=nan(size(value,1),6);
                    
                    % who's (un)defined
                    good=good & value(:,6)~=h.undef.ntype ...
                        & ~isnan(value(:,6)) & ~isinf(value(:,6));
                    
                    % skip if empty
                    if(any(good))
                        switch wf{2}
                            case '6utc'
                                head(good,:)=fixtimes(ref6(good,:)...
                                    +value(good,:),'utc');
                            case '6tai'
                                head(good,:)=fixtimes(utc2tai(...
                                    ref6(good,:))+value(good,:));
                        end
                    end
                    head=mat2cell(head,ones(nrecs,1)); type=1;
                    return;
                end
            end
        end
    end
end

% field not found
warning('seizmo:getheader:fieldInvalid',...
    'Filetype: %s, Version: %d\nInvalid field: %s',h.filetype,h.version,f);
head=nan;
type=0;

end
