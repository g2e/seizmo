function [varargout]=getheader(data,varargin)
%GETHEADER    Get SEIZMO data header values
%
%    Usage:  headers=getheader(data)
%            values=getheader(data,'field1')
%            [values1,...,valuesN]=getheader(data,'field1',...,'fieldN')
%
%    Description: GETHEADER(DATA) will return header values for all records
%     in DATA as a single numeric array.  Rows in the output array
%     correspond to the header values of individual records.  The order of
%     the fields follows that of how they are stored in the SEIZMO struct
%     (see SEIZMO's function SEISMODEF for details).  Character fields are 
%     returned as a series of their ascii number equivalents (0-255).
%
%     GETHEADER(DATA,FIELD) returns the specified header field FIELD's
%     values for each record stored in the SEIZMO data structure DATA.
%     FIELD must be a string corresponding to a valid header field or a
%     valid group field (ie. t,kt,resp,user,kuser).  Values are returned in
%     numeric arrays or cellstring arrays oriented such that each column
%     corresponds to an individual header field and each row to an
%     individual record.  So the group field 't' would return a numeric
%     array with 10 columns and as many rows as records in DATA while group
%     field 'kuser' would return a cellstring array with 3 columns and as
%     many rows as records in DATA.
%     
%     GETHEADER(DATA,FIELD1,...,FIELDN) returns one array of values per
%     field or group field.
%
%    Notes:
%     - Enumerated fields return the value actually stored, an integer used
%       to look up the enum string id and description in a table.  To 
%       retrieve the associated string id or description use the functions 
%       GETENUM or GETENUMDESC.
%     - Logical fields return the value actually stored, not a logical.  To
%       get a more useful value use GETLGC.
%     - group fields:    t, kt, user, kuser, resp, dep, st, ev, nz, nzdttm,
%                         kname, {real group} utc, {real group} tai
%     - virtual fields:  nzmonth, nzcday, kzdttm, kzdate, kztime, z, ztai
%     - abs time fields: {real field} utc, {real field} tai
%
%    Examples:
%     Put all t series values in one array:
%      times=getheader(data,'t')
%
%     Pull just the sample rates for records (fields are case insensitive):
%      dt=getheader(data,'DeLtA')
%
%     Get the station lat and lon for records:
%      [stla,stlo]=getheader(data,'stla','STLO')
%
%     Enumerated fields return the table lookup index which is
%     the value stored in the header:
%      getheader(data,'iftype')
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 16, 2009 at 17:50 GMT

% todo:

% require at least one input
if(nargin<1)
    error('seizmo:getheader:notEnoughInputs',...
        'Not enough input arguments.')
end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% number of records
nrecs=numel(data);

% recursive section
%   breaks up dataset with multiple calls so that the
%   rest of gh deals with only one header version
[h,idx]=versioninfo(data);
nver=numel(h);
sfill=repmat({'NaN'},nrecs,1);
nfill=nan(nrecs,1);
if(nver>1)
    nout=max([1 nargin-1]);
    for i=1:nver
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
    return;
end

% pull entire header
head=[data.head];

% push out entire header
if(nargin==1); varargout{1}=head.'; return; end

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
        varargout{i}(:,j)=val;
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
            if(isfield(h.real(m).pos,wf{1}))
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
                    
                    % default output to undef
                    head=ones(size(value,1),5)*h.undef.ntype;
                    
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
                    
                    % default output to undef
                    head=ones(size(value,1),6)*h.undef.ntype;
                    
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
