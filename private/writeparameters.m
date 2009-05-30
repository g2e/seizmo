function [data]=writeparameters(data,varargin)
%WRITEPARAMETERS    Implements options passed to WRITE functions
%
%    Usage:    data=writeparameters(data,...,'name',name,...)
%              data=writeparameters(data,...,'prepend',string,...)
%              data=writeparameters(data,...,'append',string,...)
%              data=writeparameters(data,...,'delete',string,...)
%              data=writeparameters(data,...,'delete',{str1 ... strN},...)
%              data=writeparameters(data,...,'change',{orig replace},...)
%              data=writeparameters(data,...,'path',path,...)
%              data=writeparameters(data,...,'pathprepend',string,...)
%              data=writeparameters(data,...,'pathappend',string,...)
%              data=writeparameters(data,...,'pathdelete',string,...)
%              data=writeparameters(data,...,'pathdelete',{s1 ... sN},...)
%              data=writeparameters(data,...,'pathchange',{orig repl},...)
%              data=writeparameters(data,...,'byteorder',endianness,...)
%
%    Description: WRITEPARAMETERS is a wrapper function for several
%     funcitons that alter top level fields of a SEIZMO data structure such
%     as the 'name', 'path' and 'byteorder' fields.  This is intended to
%     provide the functions WRITESEIZMO and WRITEHEADER with a large amount
%     of options.  For options 'NAME', 'PREPEND', 'APPEND', 'DELETE', and
%     'CHANGE' see CHANGENAME for usage.  For options 'PATH',
%     'PATHPREPEND', 'PATHAPPEND', 'PATHDELETE', and 'PATHCHANGE' see
%     CHANGEPATH for usage.  For option 'BYTEORDER' see CHANGEBYTEORDER for
%     usage.
%
%    Notes:
%
%    Header changes: NONE
%
%    Examples:
%     Read in a dataset, merge, and write out with new names:
%      writeseizmo(merge(readseizmo('*')),'append','.merged')
%
%    See also: changebyteorder, changename, changepath

%     Version History:
%        May  29, 2009 - initial version
%
%     Testing History:
%        r72 - Linux Matlab (r2007b)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  29, 2009 at 21:55 GMT

% check nargin
if(mod(nargin-1,2))
    error('seizmo:writeparameters:badNumInputs',...
        'Bad number of arguments!');
end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

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
    elseif(strncmpi(varargin{i},'path',4))
        % remove path from string unless string is path
        if(~strcmpi(varargin{i},'path'))
            varargin{i}=varargin{i}(5:end);
        end
        pathidx([i i+1])=true;
    else
        nameidx([i i+1])=true;
    end
end

% make name and path changes
data=changename(data,varargin{nameidx});
data=changepath(data,varargin{pathidx});

end
