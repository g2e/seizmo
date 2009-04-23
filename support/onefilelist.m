function [list]=onefilelist(varargin)
%ONEFILELIST    Compiles multiple filelists into one
%
%    Usage:    list=onefilelist(list1,...,listN)
%
%    Description:  ONELIST(ARG1,...,ARGN) combines the input char/cellstr 
%     arrays into a single-column cellstr array and then pipes each cell 
%     through XDIR to expand wildcard inputs.  This is mostly useful for 
%     taking file lists and wildcards in various array formats and 
%     combining them into one simple list of files.
%
%    Notes:
%     - XDIR is derived from RDIR (taken from the Matlab File Exchange)
%
%    Tested on: Matlab r2007b
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
%        Dec.  2, 2008 - use XDIR to expand wildcard ability
%        Apr. 23, 2009 - move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 22:15 GMT

% todo:

% check nargin
if(nargin<1)
    error('seizmo:onefilelist:notEnoughInputs',...
        'Not enough input arguments.');
end

% check, organize, and compile char/cellstr arrays
for i=1:nargin
    % check that input is char or cellstr array
    if(~ischar(varargin{i}) && ~iscellstr(varargin{i}))
        error('seizmo:onefilelist:badInput',...
            'Inputs must be character and cellstr arrays only!');
    end
    
    % char to cellstr then array to column vector
    varargin{i}=cellstr(varargin{i});
    varargin{i}=varargin{i}(:).';
end

% concatinate arguments
varargin=[varargin{:}].';

% pump each cell through xdir
list=[];
for i=1:length(varargin)
    % get this filelist
    files=xdir(varargin{i});
    
    % build list
    list=[list; files(~[files.isdir])]; %#ok<AGROW>
end

end
