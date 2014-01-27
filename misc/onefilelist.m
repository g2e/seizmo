function [list]=onefilelist(varargin)
%ONEFILELIST    Compiles multiple filelists into one
%
%    Usage:    list=onefilelist(list1,...,listN)
%              list=onefilelist
%
%    Description:
%     LIST=ONEFILELIST(ARG1,...,ARGN) combines input char/cellstr arrays
%     into a single-column cellstr array and then pipes each cell through
%     XDIR to expand wildcard inputs and to get detailed info for each
%     file.  This is mostly useful for taking file lists and wildcards in
%     various array formats and combining them into one structured list of
%     files.  See XDIR for the format of LIST.
%
%     LIST=ONEFILELIST presents a gui to select files.
%
%    Notes:
%     - The global SEIZMO.ONEFILELIST.FILTERSPEC allows access to the
%       filetype filtering in the gui menu.  See UIGETFILE for details.  An
%       all files spec is added to the bottom of the list automatically.
%
%    Examples:
%     % Compile a list of files from several directories
%     list=onefilelist('*','../*','~/*')
%
%    See also: STRNLEN, XDIR

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
%        Sep.  8, 2009 - minor doc fix
%        Oct. 16, 2009 - file select ui if no input
%        Jan. 27, 2010 - file select ui if empty first arg & only 1 arg
%        Feb. 14, 2010 - added SEIZMO global access to filterspec
%        Feb. 14, 2012 - LASTDIRECTORY saves last dir in gui, doc update
%        Jan. 27, 2014 - fix bad path bug in windows, reduce filesep calls
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 16:05 GMT

% todo:

% retrieve global settings
global SEIZMO

% directory separator
fs=filesep;

% grab filterspec
filterspec={}; lastdir=[pwd fs];
try
    if(~isempty(SEIZMO.ONEFILELIST.FILTERSPEC) ...
            && iscellstr(SEIZMO.ONEFILELIST.FILTERSPEC) ...
            && size(SEIZMO.ONEFILELIST.FILTERSPEC,2)==2)
        filterspec=SEIZMO.ONEFILELIST.FILTERSPEC;
    end
    if(~isempty(SEIZMO.ONEFILELIST.LASTDIRECTORY) ...
            && ischar(SEIZMO.ONEFILELIST.LASTDIRECTORY) ...
            && isdir(SEIZMO.ONEFILELIST.LASTDIRECTORY))
        lastdir=SEIZMO.ONEFILELIST.LASTDIRECTORY;
    end
catch
end

% check nargin
if(nargin<1 || (nargin==1 && isempty(varargin{1})))
    [files,path]=uigetfile(...
        [filterspec; {'*.*;*' 'All Files (*.* , *)'}],...
        'Select File(s)',lastdir,'MultiSelect','on');
    if(~isempty(path) && ischar(path))
        SEIZMO.ONEFILELIST.LASTDIRECTORY=[path fs];
    end
    if(isequal(0,files))
        error('seizmo:onefilelist:noFilesSelected','No files selected!');
    end
    varargin=strcat(path,fs,cellstr(files));
else
    % check, organize, and compile char/cellstr arrays
    for i=1:max(1,nargin)
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
end

% pump each cell through xdir
list=[];
for i=1:numel(varargin)
    % get this filelist
    files=xdir(varargin{i});
    
    % build list
    list=[list; files(~[files.isdir])]; %#ok<AGROW>
end

end
