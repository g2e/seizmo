function [data]=sortbyfield(data,field,mode)
%SORTBYFIELD   Sort SEIZMO records by a header or SEIZMO struct field
%
%    Usage:    data=sortbyfield(data,field)
%              data=sortbyfield(data,field,mode)
%
%    Description:
%     DATA=SORTBYFIELD(DATA,FIELD) sorts records in SEIZMO struct DATA by
%     the header field FIELD.  Also will sort by any top-level struct field
%     in DATA such as 'name', 'version', 'byteorder', etc.  Data fields
%     'misc', 'head', 'dep', and 'ind' are not supported.  Group header
%     fields are also not allowed.
%
%     SORTBYFIELD(DATA,FIELD,MODE) sets the sorting order ('ascend' or
%     'descend' is allowed - 'ascend' is the default).
%
%    Notes:
%
%    Header changes: NONE
%
%    Examples: 
%     % Sort by descending degree distance:
%     data=sortbyfield(data,'gcarc','descend')
%
%    See also: SORT, GETHEADER

%     Version History:
%        Oct. 31, 2007 - initial version
%        Nov.  7, 2007 - doc update
%        Jan. 28, 2008 - code cleaning
%        Feb. 28, 2008 - better checking, doc update, any DATA field
%        Mar.  4, 2008 - minor doc update
%        Nov. 23, 2008 - updated for new name schema (now SORTBYHEADER),
%                        history fix
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        Oct. 13, 2009 - minor doc update, added .misc to bad fields
%        Feb.  3, 2010 - proper SEIZMO handling, versioninfo caching,
%                        seizmoverbose support
%        Feb. 11, 2011 - mass nargchk fix, dropped versioninfo caching
%        Apr.  2, 2012 - minor doc update, use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 15:05 GMT

% todo:

% check number of args
error(nargchk(2,3,nargin));

% check data structure
error(seizmocheck(data));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt sort
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % detail message
    if(verbose)
        disp('Sorting Record(s)');
        print_time_left(0,nrecs);
    end
    
    % set mode if none
    if(nargin==2 || isempty(mode)); mode='ascend'; end

    % get field values (or filenames/byte-orders/versions)
    bad={'head' 'dep' 'ind' 'misc'};
    if(~ischar(field) || any(strcmpi(field,bad)))
        error('seizmo:sortbyfield:badField','FIELD is bad!');
    elseif(isfield(data,field))
        if(isnumeric([data.(field)]))
            [values,indices]=sort([data.(field)]);
        else
            [values,indices]=sort({data.(field)});
        end
    else
        [values,indices]=sort(getheader(data,field));
    end

    % check indices size
    if(numel(indices)~=nrecs)
        error('seizmo:sortbyfield:tooManyIndices',...
            'Too many elements to sort by!')
    end

    % flip if descend mode
    if(strcmpi(mode,'descend')); indices=indices(end:-1:1); end

    % sort data
    data=data(indices);
    
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
