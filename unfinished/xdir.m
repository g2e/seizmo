function [varargout]=xdir(str,depth)
%XDIR    Cross-Platform directory listing with recursion
%
%    Usage:      xdir()
%                xdir(path)
%                xdir(path,depth)
%              d=xdir(...)
%
%    Description: XDIR() lists the files/directories in the current
%     directory. The listing gives the directory logical first (true for
%     directories), then the file size in bytes, followed by the date and
%     time of the last modification and finally the name of the directory
%     or file.
%
%     XDIR(PATH) lists the files matching PATH or within the directory
%     PATH.  The wildcard '*' or '?' may be used to specify a directory or
%     filename pattern (Matlab's DIR only allows wildcards in the name).
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
%     XDIR(PATH,DEPTH) allows specifying the recursion limit for all '**'
%     elements in PATH.  If DEPTH is a single integer it is taken as the
%     maximum recursion depth with the minimum set to 0.  Specifying 0
%     (zero) will only look for matches in the current directory and will
%     not recurse.  Giving two integers specifies a recursion depth range.
%     For example, a range of [2 3] will only return matches where the '**'
%     elements have recursed 2 to 3 times.  It is not possible to specify
%     recursion limits for each '**' element separately.  The default is
%     [0 Inf].
%
%     D=XDIR(...) returns the results in an M-by-1 structure with the
%     fields:
%       name    -- filename or directory name
%       date    -- modification date
%       bytes   -- number of bytes allocated to the file or directory
%       isdir   -- 1 if is a directory, 0 otherwise
%       datenum -- modification date as a MATLAB serial date number
%       path    -- path to the file or directory
%
%    Notes:
%     - XDIR uses self-recursion to extend down directory structures, so
%       be aware of Matlab's recursion limit (typically 50) if accessing
%       deep trees.
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
%        Apr. 20, 2009 - add filesep to end of path for non-wildcard
%                        directory expansion
%        Apr. 23, 2009 - made to work with Octave too (fixes for assigning
%                        to struct and adding fields to empty struct)
%        June  4, 2009 - use GLOB to preempt or replace octave's dir
%        Aug.  4, 2009 - fix dir listing for Octave when no path given
%
%     Written by Gus Brown ()
%                Garrett Euler (ggeuler at seismo dot wustl dot edu)
%     Last Updated Aug.  4, 2009 at 09:45 GMT

% todo:
% - be mindful of octave/matlab differences in DIR
%   - wildcards before last pathsep
%       - matlab doesn't expand
%       - octave does
%       - xdir accounts for this!
%   - ? wildcard
%       - matlab doesn't treat as wild
%       - octave does
%       - xdir does NOT account for this!
%   - if wildcard entry is a directory
%       - octave expands the entry
%       - matlab lists the entry but does not expand it
%       - xdir accounts for this discrepancy!
%
% - handle '?' wildcard (only for matlab - octave already does this!)
%  * => .*     42 => 46 42
%  ? => .      63 => 46
%  . => \.     46 => 92 46
% \? => \?  92 63 => 92 63
%  \ => \\     92 => 92 92
%  ] => \]     93 => 92 93
%  [ => \[     91 => 92 91
%  $ => \$     36 => 92 36
%
%  *     ?     .     [     ]     $     \
% 42    63    46    91    93    36    92
%
% c:\/?\*  how to do this?
%
%   - something along the lines of:
%      d=dir('*')
%      for ...
%        regexp(filename,regexptranslate('wildcard',string),'match')
%      end
%   - but how to allow for '?' in a filename too?
%   - really we just need a glob for matlab too!

% take care of common tasks
persistent emptylist octave
if(isempty(octave) || ~islogical(octave))
  octave=strcmpi(getapplication,'OCTAVE');
end
if(~isstruct(emptylist))
  % make empty dir struct with path added in
  blah=[{'path'}; fieldnames(dir(''))];
  blah=[blah cell(size(blah))].';
  tmp([],1)=struct(blah{:});
  emptylist=tmp;
end

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
  W=find(newstr=='*' | newstr=='?',1,'first');
  
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

% if no path wildcards then just get listing
D=emptylist;
if isempty(wildpath)
  % handle dir expanding directories
  if(isdir([prepath postpath]))
    % no wildcards and is a directory
    % list directory and clean up ending filesep
    path=[prepath postpath];
    D=dir(path);
    if(strcmp(path(end),filesep)); path=path(1:end-1); end
    
    % workaround for new field to empty struct
    if(isempty(D))
      D=emptylist;
    else
      [D.path]=deal([path filesep]);
    end
  else
    % handle difference between octave/matlab
    if(octave)
      % list matches, list dir, get matches
      % - dir expands single dir case
      % - glob does not, so use glob
      tmp=glob([prepath postpath]);
      if(isempty(prepath))
        D=dir;
      else
        D=dir(prepath);
      end
      D=D(ismember(strcat(prepath,{D.name}),tmp));
      
      % workaround for new field to empty struct
      if(isempty(D))
        D=emptylist;
      else
        [D.path]=deal(prepath);
      end
    else % matlab
      % either a file or wildcarded files/dirs
      
      % get list of all in dir
      
      % translate postpath
      
      % find matches
      
      D=dir([prepath postpath]);
      % workaround for new field to empty struct
      if(isempty(D))
        D=emptylist;
      else
        [D.path]=deal(prepath);
      end
    end
  end
% a double wild directory means recurse down into sub directories
elseif strcmp(wildpath,'**')
  % first look for files in the current directory (remove extra filesep)
  if(depth(1)<1)
    D=xdir([prepath postpath(2:end)],depth);
  end

  % then look for sub directories (if recursion limit is not exceeded)
  if(depth(2))
    % handle difference between octave/matlab
    if(octave)
      tmp=glob([prepath '*']);
      
      % process each directory found
      for i=1:numel(tmp)
        if(isdir(tmp{i}) ...
            && ~strcmp(tmp{i}(end-1:end),[filesep '.']) ...
            && ~strcmp(tmp{i}(end-2:end),[filesep '..']))
          D=[D; xdir([tmp{i} filesep '**' postpath],depth-1)]; %#ok<AGROW>
        end
      end
    else % matlab
      % list directory
      tmp=dir([prepath '*']);
      
      % process each directory found
      for i=1:numel(tmp)
        if(tmp(i).isdir ...
            && ~strcmp(tmp(i).name,'.') && ~strcmp(tmp(i).name,'..'))
          D=[D; xdir([prepath tmp(i).name filesep '**' postpath],depth-1)]; %#ok<AGROW>
        end
      end
    end
  end
else
  % handle difference between octave/matlab
  if(octave)
    tmp=glob([prepath wildpath]);
    
    % process each directory found
    for i=1:numel(tmp)
      if(isdir(tmp{i}) ...
          && ~strcmp(tmp{i}(end-1:end),[filesep '.']) ...
          && ~strcmp(tmp{i}(end-2:end),[filesep '..']))
        D=[D; xdir([tmp{i} postpath],depth)]; %#ok<AGROW>
      end
    end
  else % matlab
    % Process directory wild card looking for sub directories that match
    tmp=dir([prepath wildpath]);
    
    % process each directory found
    for i=1:numel(tmp)
      if(tmp(i).isdir ...
          && ~strcmp(tmp(i).name,'.') && ~strcmp(tmp(i).name,'..'))
        D=[D; xdir([prepath tmp(i).name postpath],depth)]; %#ok<AGROW>
      end
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

end

function [application,version]=getapplication()
%GETAPPLICATION    Returns application running this script and its version
%
%    Usage:    [application,version]=getapplication()
%
%    Description: [APPLICATION,VERSION]=GETAPPLICATION() will determine and
%     return the name and version of the application running this script
%     (obviously only if the application can run this script in the first
%     place).  Both APPLICATION and VERSION are strings.
%
%    Notes:
%     - returns 'UNKNOWN' if it cannot figure out the application
%
%    Examples:
%     Matlab and Octave still behave quite differently for a number of
%     different functions so it is best in some cases to use different
%     function calls depending on which we are running:
%      [app,ver]=getapplication;
%      if(strcmp(app,'MATLAB'))
%        % do something via matlab routines
%      else
%        % do something via octave routines
%      end
%
%    See also: nativebyteorder, ver

%     Version History:
%        Nov. 13, 2008 - initial version
%        Mar.  3, 2009 - minor doc cleaning
%        Apr. 23, 2009 - move usage up
%        Sep.  8, 2009 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  8, 2009 at 19:55 GMT

% todo:

% checking for Matlab will throw an error in Octave
try
    % first check if we are in Matlab
    a=ver('matlab');
    
    % we are in Matlab
    application=a.Name;
    version=a.Version;
    return;
catch
    % check if we are in Octave
    if(exist('OCTAVE_VERSION','builtin')==5)
        application='OCTAVE';
        version=OCTAVE_VERSION;
        return;
    % ok I have no clue what is running
    else
        application='UNKNOWN';
        version='UNKNOWN';
        return;
    end
end

end

function [string]=translator(string)
% - handle '?' wildcard (only for matlab - octave already does this!)
%
%  * => .*     42 => 46 42
%  ? => .      63 => 46
%
%  . => \.     46 => 92 46
% \? => \?  92 63 => 92 63
% \* => \*  92 42 => 92 42
%  \ => \\     92 => 92 92
%  ] => \]     93 => 92 93
%  [ => \[     91 => 92 91
%  $ => \$     36 => 92 36
%
%  *     ?     .     [     ]     $     \
% 42    63    46    91    93    36    92
%
% c:\???  how to do this?
% - avoid so it should never happen!!

% position of chars to be escaped
period=string==46;
brkbgn=string==91;
brkend=string==93;
dollar=string==36;

% find *,?,\
stars=string==42;
quest=string==63;
slash=string==92;

% remove \* from *,  \? from ?,  \* & \? from \
sp=regexp(string,'\\\*');
qp=regexp(string,'\\\?');
stars(find(stars)-1==sp)=false;
quest(find(quest)-1==qp)=false;
slash(find(slash)==sp | find(slash)==qp)=false;

% push string into cell array with each char in a separate cell
string=num2cell(string);

end
