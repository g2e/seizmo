function [varargout]=xdir(str)
%XDIR    Directory listing with recursion
%
%    Description: XDIR(PATH) lists the files matching PATH or within the
%     directory PATH.  The wildcard '*' may be used to specify a directory
%     or filename pattern (DIR only allows wildcards for the filename).
%     For example, XDIR */*.m lists all the M-files in the subdirectories
%     one level below the current directory.  The directory wildcard '**'
%     is used to match multiple directory levels for a file or directory
%     (including the current directory).  For example, XDIR **/*.m would
%     list all the M-files within the current directory as well as those
%     within all subdirectories below the current directory.  The listing
%     gives the directory logical first, then the size in bytes, followed
%     by the date and time of the last modification and finally the path to
%     the directory or file.
%
%     D=XDIR(PATH) returns the results in an M-by-1 structure with the
%     fields:
%       name    -- filename
%       date    -- modification date
%       bytes   -- number of bytes allocated to the file
%       isdir   -- 1 if name is a directory, 0 otherwise
%       datenum -- modification date as a MATLAB serial date number
%       path    -- path to the file or directory
%
%    Notes:
%     - DIR, and thus XDIR, do not handle properly filenames with a
%       wildcard ('*' or '**') or a file separator ('/' or '\') within
%       their name - they cannot be escaped.
%     - XDIR uses self-recursion to extend down directory structures, so
%       be aware of your recursion limit (typically 50) if accessing deep
%       trees
%
%    Tested on: Matlab r2007b
%
%    Usage:      xdir(path)
%              d=xdir(path)
%
%    Examples:
%     List all m-files in private directories below the current directory:
%      xdir **/private/*.m
%
%    See also: dir, ls, cd, delete, rmdir, mkdir

%     Version History:
%        Apr. 14, 2008 - initial version (Matlab FEX #19550)
%        Dec.  2, 2008 - bug fix for dir expansion of non-wildcard names,
%                        collapse multiple file separators, slight changes
%                        to the coding style, changed the listing to be a
%                        little more informative, added field 'path' to
%                        separate the path from the name, RDIR => XDIR
%
%     Written by Gus Brown ()
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  2, 2008 at 09:15 GMT

% todo:
% - specify a max recursion depth
% - handle '?' wildcard

% use the current directory if nothing is specified
if(nargin==0); str='*'; end

% collapse multiple file separators
if(strncmp(computer,'PC',2))
  newstr=regexprep(str,'\\{2,}','\');
  filesep='\';
else
  newstr=regexprep(str,'/{2,}','/');
  filesep='/';
end

% split the file path around the first wild card specifier
prepath='';       % the path before the wild card
wildpath='';      % the path wild card
postpath=newstr;  % the path after the wild card
I=find(newstr==filesep,1,'last');
if(~isempty(I))
  prepath=newstr(1:I);
  postpath=newstr(I+1:end);
  I=find(prepath=='*',1,'first');
  if(~isempty(I))
    postpath=[prepath(I:end) postpath];
    prepath=prepath(1:I-1);
    I=find(prepath==filesep,1,'last');
    if(~isempty(I))
      wildpath=prepath(I+1:end);
      prepath=prepath(1:I);
    end
    I=find(postpath==filesep,1,'first');
    if(~isempty(I))
      wildpath=[wildpath postpath(1:I-1)];
      postpath=postpath(I:end);
    end
  end
end

% if no directory wildcards then just get file list
if isempty(wildpath)
  % handle dir expanding directories
  if(isdir([prepath postpath]))
    D=dir([prepath postpath]);
    [D(:).path]=deal([prepath postpath]);
  else
    D=dir([prepath postpath]);
    [D(:).path]=deal(prepath);
  end
% a double wild directory means recurse down into sub directories
elseif strcmp(wildpath,'**')
  % first look for files in the current directory (remove extra filesep)
  D=xdir([prepath postpath(2:end)]);

  % then look for sub directories
  tmp=dir([prepath '*']);
  % process each directory
  for i=1:numel(tmp),
    if (tmp(i).isdir ...
        && ~strcmp(tmp(i).name,'.') && ~strcmp(tmp(i).name,'..') ),
      D=[D; xdir([prepath tmp(i).name filesep '**' postpath])]; %#ok<AGROW>
    end
  end
else
  % Process directory wild card looking for sub directories that match
  tmp=dir([prepath wildpath]);
  D=dir(''); [D(:).path]=deal([]);
  % process each directory found
  for i=1:numel(tmp),
    if((tmp(i).isdir ...
        && ~strcmp(tmp(i).name,'.') && ~strcmp(tmp(i).name,'..')))
      D=[D; xdir([prepath tmp(i).name postpath])]; %#ok<AGROW>
    end
  end
end

% display listing if no output variables are specified
if(nargout==0)
  for i=1:numel(D) 
      disp(sprintf('%d %16d %20s %-64s',...
          D(i).isdir,D(i).bytes,D(i).date,[D(i).path D(i).name]));
  end
else
  % send list out
  varargout{1}=D;
end

