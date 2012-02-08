function [data]=writeparameters(data,varargin)
%WRITEPARAMETERS    Implements options passed to WRITE functions
%
%    Usage:    data=writeparameters(data,...,'name',name,...)
%              data=writeparameters(data,...,'prepend',string,...)
%              data=writeparameters(data,...,'append',string,...)
%              data=writeparameters(data,...,'delete',string,...)
%              data=writeparameters(data,...,'change',{orig replace},...)
%              data=writeparameters(data,...,'path',path,...)
%              data=writeparameters(data,...,'pathprepend',string,...)
%              data=writeparameters(data,...,'pathappend',string,...)
%              data=writeparameters(data,...,'pathdelete',string,...)
%              data=writeparameters(data,...,'pathchange',{orig repl},...)
%              data=writeparameters(data,...,'byteorder',endianness,...)
%
%    Description:
%     WRITEPARAMETERS is a wrapper function for several functions that
%     alter top level fields of a SEIZMO data structure such as the 'name',
%     'path' and 'byteorder' fields.  This is intended to provide the
%     functions WRITESEIZMO and WRITEHEADER with a large amount of options.
%     For options 'NAME', 'PREPEND', 'APPEND', 'DELETE', and 'CHANGE' see
%     CHANGENAME for usage.  For options 'PATH', 'PATHPREPEND',
%     'PATHAPPEND', 'PATHDELETE', and 'PATHCHANGE' see CHANGEPATH for
%     usage.  For option 'BYTEORDER' see CHANGEBYTEORDER for usage.
%
%    Notes:
%
%    Header changes: NONE
%
%    Examples:
%     % Read in a dataset, merge, and write out with new names:
%     writeseizmo(meld(readseizmo('*')),'append','.merged')
%
%    See also: CHANGEBYTEORDER, CHANGENAME, CHANGEPATH, WRITESEIZMO,
%              WRITEHEADER

%     Version History:
%        May  29, 2009 - initial version
%        Feb.  2, 2010 - proper SEIZMO handling
%        Apr. 25, 2010 - allow options like 'changepath', 'appendpath', etc
%        Feb. 11, 2011 - mass seizmocheck fix
%        Feb.  7, 2012 - merge to meld update, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  7, 2012 at 15:05 GMT

% check nargin
if(mod(nargin-1,2))
    error('seizmo:writeparameters:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
error(seizmocheck(data));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt name/path/endian changes
try
    % get options from command line
    pathidx=false(nargin-1,1);
    nameidx=false(nargin-1,1);
    for i=1:2:nargin-1
        if(~ischar(varargin{i}))
            error('seizmo:writeparameters:badInput',...
                'Options must be specified as a string!');
        end
        if(strcmpi(varargin{i},'byteorder'))
            data=changebyteorder(data,varargin{i+1});
        elseif(~isempty(strfind(varargin{i},'path')))
            pathidx([i i+1])=true;
        else
            nameidx([i i+1])=true;
        end
    end

    % make name and path changes
    if(any(nameidx))
        data=changename(data,varargin{nameidx});
    end
    if(any(pathidx))
        data=changepath(data,varargin{pathidx});
    end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
