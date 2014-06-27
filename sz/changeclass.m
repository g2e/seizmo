function [data,origclass]=changeclass(data,newclass)
%CHANGECLASS    Change SEIZMO data storage (in memory)
%
%    Usage:    data=changeclass(data,newclass)
%              [data,origclass]=changeclass(data,newclass)
%
%    Description:
%     DATA=CHANGECLASS(DATA,NEWCLASS) changes the in-memory class type of
%     the independent and dependent data stored in the SEIZMO dataset DATA
%     to NEWCLASS.  NEWCLASS must be a string or cellstr of valid class(es)
%     ('double', 'single', etc).  This does not change the storage type of
%     data written to disk (requires a filetype/version change).
%
%     [DATA,ORIGCLASS]=CHANGECLASS(DATA,NEWCLASS) also returns the original
%     class type of the dependent data as a cellstr array.
%
%    Notes:
%     - Changing the storage type of written data requires a
%       filetype/version change (using CHANGEFILETYPE)
%
%    Header changes: NONE
%
%    Examples:
%     % Return double precision records and attempt to fix sample spacing:
%     data=fixdelta(changeclass(data,'double'))
%
%    See also: FIXDELTA

%     Version History:
%        Feb. 21, 2008 - initial version
%        Feb. 23, 2008 - uses GLGC now
%        Feb. 28, 2008 - uses SEISCHK now
%        Mar.  4, 2008 - doc update, fixed LEVEN bug, uses LGCCHK now
%        June  6, 2008 - minor code cleaning
%        Oct.  8, 2008 - changed name from DOUBLEIT to CLASSIT, doc update,
%                        dropped reclassing header (keep it double!), allow
%                        data to be changed to any Matlab supported class,
%                        drop LGCCHK
%        Nov. 22, 2008 - renamed from CLASSIT to CHANGECLASS
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Oct.  5, 2009 - doc update
%        Jan. 28, 2010 - seizmoverbose support
%        Aug. 16, 2010 - fix error usage
%        Apr.  3, 2012 - minor doc update
%        June 26, 2014 - also allow returning of original data class
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 26, 2014 at 14:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% verbosity
verbose=seizmoverbose;

% number of records
nrecs=numel(data);

% check class and convert to function handle
if(isa(newclass,'function_handle'))
    if(isscalar(newclass))
        newclass(1:nrecs,1)={newclass};
    elseif(numel(newclass)~=nrecs)
        error('seizmo:changeclass:badInput',...
            'NEWCLASS must be scalar or match number of records in DATA!');
    else
        newclass=mat2cell(newclass(:),ones(nrecs,1));
    end
elseif(ischar(newclass) || iscellstr(newclass))
    newclass=cellstr(newclass);
    for i=1:numel(newclass)
        newclass(i)={str2func(newclass{i})};
    end
else
    error('seizmo:changeclass:badInput',...
        'NEWCLASS must be a string, cellstr, or function handle!');
end
if(isscalar(newclass))
    newclass(1:nrecs,1)=newclass;
elseif(numel(newclass)~=nrecs)
    error('seizmo:changeclass:badInput',...
        'NEWCLASS must be scalar or match number of records in DATA!');
end

% detail message
if(verbose)
    disp('Changing Class of Dependent Data in Record(s)');
    print_time_left(0,nrecs);
end

% preallocate original class type output
origclass=cell(numel(newclass),1);

% reclass dependent data
for i=1:nrecs
    % get original class
    origclass{i}=class(data(i).dep);
    
    % reclass
    data(i).dep=newclass{i}(data(i).dep);
    
    % detail message
    if(verbose); print_time_left(i,nrecs); end
end

% detail message
if(verbose)
    disp('Changing Class of Independent Data in Record(s)');
    print_time_left(0,nrecs);
end

% reclass independent data
if(isfield(data,'ind'))
    for i=1:nrecs
        data(i).ind=newclass{i}(data(i).ind);
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
else
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end
end

end
