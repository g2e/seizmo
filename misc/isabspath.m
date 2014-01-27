function [lgc]=isabspath(path,iswindows)
%ISABSPATH    Determines if a path is an absolute path or not
%
%    Usage:    lgc=isabspath(path)
%              lgc=isabspath(path,iswindows)
%
%    Description:
%     LGC=ISABSPATH(PATH) checks if the path(s) in PATH are relative or
%     absolute and returns TRUE for those that are absolute paths.  PATHS
%     may be a string, char array or a cell string array.  LGC is a logical
%     array with one element per path in PATH.  This is useful to find
%     relative paths so you can convert them to absolute paths for
%     functions like EXIST.  The determination is done by discovering the
%     OS type of the current system using ISPC.
%
%     LGC=ISABSPATH(PATH,ISWINDOWS) allows setting the OS type for
%     determining if the paths are absolute or not when the paths are not
%     valid paths for the current machine.  For instance, set ISWINDOWS to
%     FALSE for Unix, Linux or MACOSX paths when you are using MicrosoftTM
%     WindowsTM.  ISWINDOWS must be TRUE or FALSE (scalar only).
%
%    Notes:
%     - The path is not required to exist or even to be valid!  This just
%       does a simple test on each path given the OS (e.g., is the first
%       character a '/' for unix).
%
%    Examples:
%     % Test a few relative paths:
%     isabspath('./somedir')
%     isabspath('../somedir')
%     isabspath('~/somedir')
%     isabspath('..\somewindir')
%
%     % Test a few absolute paths:
%     isabspath('/home')
%     isabspath('/usr/share/../bin')
%     isabspath('c:\Programs')
%
%     % And a few invalid ones:
%     isabspath('/\')                   % absolute path to the '\' dir?
%     isabspath('somedir\c:/somewhere') % win drive in a unix dir under pwd
%     isabspath('\\someserver\somedir') % maybe you can add this feature...
%
%    See also: ISPC, ISUNIX

%     Version History:
%        Jan. 27, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 11:15 GMT

% todo:

% check number of inputs
error(nargchk(1,2,nargin));

% check/fix path
if(ischar(path))
    path=cellstr(path);
elseif(~iscellstr(path))
    error('seizmo:isabspath:badInput',...
        'PATH must be a string, char array or a cellstr array!');
end

% check/default os
if(nargin<2 || isempty(iswindows)); iswindows=ispc; end
if(~islogical(iswindows) || ~isscalar(iswindows))
    error('seizmo:isabspath:badInput',...
        'ISWINDOWS must be TRUE or FALSE!');
end

% preallocate output as all relative paths
lgc=false(size(path));

% act by os
if(iswindows) % windows
    for i=1:numel(path)
        % require drive char to be a-z,A-Z
        if(isempty(path{i})); continue; end
        drive=double(upper(path{i}(1)));
        lgc(i)=drive>=65 && drive<=90 && strcmp(path{i}(2:3),':\');
    end
else % unix, linux, macosx
    for i=1:numel(path)
        if(isempty(path{i})); continue; end
        lgc(i)=strcmp(path{i}(1),'/');
    end
end

end
