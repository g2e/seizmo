function [list]=onefilelist(varargin)
%ONEFILELIST    Compiles multiple filelists into one
%
%    Description:  ONELIST(ARG1,...,ARGN) combines the input char/cellstr 
%     arrays into a single-column cellstr array and then pipes each cell 
%     through DIR to expand wildcard inputs.  This is mostly useful for 
%     taking file lists and wildcards in various array formats and 
%     combining them into one simple list of files.
%
%    Notes:
%     - does not support cell arrays of cellstr arrays (ie. no recursion!)
%     - directory entries are not expanded and will be omitted from the
%       returned list
%     - uses Matlab's DIR command to get the filelist and file info, which
%       has limited wildcard capability (only * is a valid wildcard and the
%       wildcards can only be applied to filenames, not the path to the
%       file)
%
%    Tested on: Matlab r2007b
%
%    Usage:    list=onefilelist(list1,...,listN)
%
%    Examples:
%     Compile a list of files from several directories
%      list=onefilelist('*','../*','~/*')
%
%    See also: strnlen

%     Version History:
%        Mar.  7, 2008 - initial version
%        Apr. 23, 2008 - now expands wildcards using DIR
%        Sep. 14, 2008 - doc update, input checks
%        Oct.  3, 2008 - fix paths being dropped
%        Oct.  5, 2008 - now handles directory entries in input filelist 
%                        without breaking (reads all files in directory)
%        Oct. 13, 2008 - skip directory entries in input filelist
%        Oct. 16, 2008 - fix files in current directory bug
%        Oct. 31, 2008 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 31, 2008 at 00:05 GMT

% todo:

% check nargin
if(nargin<1)
    error('SAClab:onefilelist:notEnoughInputs',...
        'Not enough input arguments.');
end

% check, organize, and compile char/cellstr arrays
for i=1:nargin
    % check that input is char or cellstr array
    if(~ischar(varargin{i}) && ~iscellstr(varargin{i}))
        error('SAClab:onefilelist:badInput',...
            'Inputs must be character and cellstr arrays only!');
    end
    
    % char to cellstr then array to column vector
    varargin{i}=cellstr(varargin{i});
    varargin{i}=varargin{i}(:).';
end

% concatinate arguments
varargin=[varargin{:}].';

% pump each cell through dir
list=[];
for i=1:length(varargin)
    % skip directories (so dir doesn't expand it)
    if(isdir(varargin{i})); continue; end
    
    % get this filelist
    files=dir(varargin{i});
    
    % add path if not empty list
    if(~isempty(files))
        pathstr=fileparts(varargin{i});
        % only add path if there is one
        if(~isempty(pathstr))
            fullfilenames=strcat(pathstr,filesep,{files.name}.');
            [files.name]=swap(fullfilenames{:});
        end
    end
    
    % build list
    list=[list; {files(~[files.isdir]).name}.']; %#ok<AGROW>
end

end
