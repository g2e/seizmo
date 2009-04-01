function [varargout]=xdir(str,depth)
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
%     D=XDIR(PATH,DEPTH) allows specifying the recursion limit for all '**'
%     elements in PATH.  If DEPTH is a single integer it is taken as the
%     maximum recursion depth with the minimum set to 0.  Specifying 0
%     (zero) will only look for matches in the current directory and will
%     not recurse.  Giving two integers specifies a recursion depth range.
%     For example, a range of [2 3] will only return matches where the '**'
%     elements have recursed 2 to 3 times.  It is not possible to specify
%     recursion limits for each '**' element separately.  The default is
%     [0 Inf].
%
%    Notes:
%     - XDIR uses self-recursion to extend down directory structures, so
%       be aware of Matlab's recursion limit (typically 50) if accessing
%       deep trees.
%
%    Tested on: Matlab r2007b
%
%    Usage:      xdir(path)
%                xdir(path,depth)
%              d=xdir(path)
%              d=xdir(path,depth)
%
%    Examples:
%     List all m-files in private directories below the current directory:
%      xdir **/private/*.m
%
%    See also: dir, ls, cd, delete, rmdir, mkdir

%     Version History:
%        Apr. 14, 2008 - initial version (RDIR - Matlab FEX #19550)
%        Dec.  2, 2008 - bug fix for dir expansion of non-wildcard names,
%                        collapse multiple file separators, slight changes
%                        to the coding style, changed the listing to be a
%                        little more informative, added field 'path' to
%                        separate the path from the name, RDIR => XDIR
%        Mar. 31, 2009 - rewrote path separation (fixes a few more bugs),
%                        now can specify recursion limits
%
%     Written by Gus Brown ()
%                Garrett Euler (ggeuler at seismo dot wustl dot edu)
%     Last Updated Mar. 31, 2009 at 17:45 GMT

% todo:
% - handle '?' wildcard
%   - something along the lines of:
%      d=dir('*')
%      for ...
%        regexp(filename,regexptranslate('wildcard',string),'match')
%      end
%   - but how to allow for '?' in a filename too?

% use the current directory if nothing is specified
if(nargin==0); str='*'; end

% check if string
if(~ischar(str))
    error('seizmo:xdir:badPath','PATH must be a string!');
end

% allow full range if depth not specified
if(nargin<2 || isempty(depth))
    depth=[0 inf];
elseif(~isreal(depth) || numel(depth)>2 || any(fix(depth)~=depth))
    error('seizmo:xdir:badDepth',...
        'DEPTH must be a 1 or 2 element array of real integers!')
elseif(isscalar(depth))
    depth=[0 depth];
end

% collapse multiple file separators
if(strncmp(computer,'PC',2))
  newstr=regexprep(str,'\\{2,}','\');
  filesep='\';
else
  newstr=regexprep(str,'/{2,}','/');
  filesep='/';
end

% assume easy case (no file separators, with or without a file wildcard)
prepath='';       % the path before the wild card
wildpath='';      % the path wild card
postpath=newstr;  % the path after the wild card

% find position of all file separators
F=find(newstr==filesep);

% if there are file separators then check for path wildcards
if(~isempty(F))
    % assume easy case (no path wildcards)
    prepath=newstr(1:F(end));
    postpath=newstr(F(end)+1:end);
    
    % find position of first wildcard
    W=find(newstr=='*',1,'first');
    
    % check for path wildcards
    if(~isempty(W))
        % split the file path around the first wild card specifier
        if(W<F(1))
            % no prepath, some wildpath, some postpath
            prepath='';
            wildpath=newstr(1:F(1)-1);
            postpath=newstr(F(1):end);
        elseif(W<F(end))
            % some prepath, some wildpath, some postpath
            L=F(find(W>F,1,'last'));
            H=F(find(W<F,1,'first'));
            prepath=newstr(1:L);
            wildpath=newstr(L+1:H-1);
            postpath=newstr(H:end);
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
  if(depth(1)<1)
    D=xdir([prepath postpath(2:end)],depth);
  else
    D=dir(''); [D(:).path]=deal([]);
  end

  % then look for sub directories (if recursion limit is not exceeded)
  if(depth(2))
    tmp=dir([prepath '*']);
    % process each directory
    for i=1:numel(tmp),
      if (tmp(i).isdir ...
          && ~strcmp(tmp(i).name,'.') && ~strcmp(tmp(i).name,'..') ),
        D=[D; xdir([prepath tmp(i).name filesep '**' postpath],depth-1)]; %#ok<AGROW>
      end
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
      D=[D; xdir([prepath tmp(i).name postpath],depth)]; %#ok<AGROW>
    end
  end
end

% to make all entries have a path (use './' or '.\' if there is none)
%[D(strcmp({D.path},'')).path]=deal(['.' filesep]);

% display listing if no output variables are specified
if(nargout==0)
  for i=1:numel(D) 
      disp(sprintf('%d %16d %20s %-s',...
          D(i).isdir,D(i).bytes,D(i).date,[D(i).path D(i).name]));
  end
else
  % send list out
  varargout{1}=D;
end

