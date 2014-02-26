function [txt,file]=readtxt(file,filterspec)
%READTXT    Reads in a text file as a single string
%
%    Usage:    txt=readtxt(file)
%              txt=readtxt([],filterspec)
%              [txt,file]=readtxt(...)
%
%    Description:
%     TXT=READTXT(FILE) reads in the ascii file given by FILE (the location
%     of the file on the filesystem) as a single string TXT.  TXT is a row
%     vector of characters.  End of line characters are not removed.
%     Calling READTXT without FILE or with FILE set to '' will present a
%     graphical file selection menu.
%
%     TXT=READTXT([],FILTERSPEC) sets the file filter specifier for gui-
%     based file selection.  The default FILTERSPEC is:
%      {'*.txt;*.TXT' 'TXT Files (*.txt,*.TXT)';
%       '*.*' 'All Files (*.*)'}
%     This is mainly so other functions do not need to include their own ui
%     to select an ascii file (that will get passed to here anyways).
%
%     [TXT,FILE]=READTXT(...) returns the file (with absolute path) as the
%     character string FILE.  This is useful when no filename was given as
%     an input.
%
%    Notes:
%
%    Examples:
%     % The purpose of READTXT is to simplify the reading of text files:
%     readtxt('somefile.txt')
%
%    See also: GETWORDS

%     Version History:
%        Dec. 30, 2009 - initial version
%        Jan. 26, 2010 - add graphical selection
%        Feb.  5, 2010 - improved file checks
%        July 30, 2010 - nargchk fix
%        Aug. 10, 2010 - filterspec option added
%        Nov.  1, 2011 - doc update
%        Jan. 26, 2014 - abs path exist fix
%        Feb.  9, 2014 - file output
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2014 at 13:00 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% directory separator
fs=filesep;

% default/check filterspec
if(nargin<2 || isempty(filterspec))
    filterspec={'*.txt;*.TXT' 'TXT Files (*.txt,*.TXT)';
        '*.*' 'All Files (*.*)'};
end
if(~iscellstr(filterspec))
    error('seizmo:readtxt:badFilterSpec',...
        'FILTERSPEC must be a cellstr array for UIGETFILE!');
end

% graphical selection
if(nargin<1 || isempty(file))
    [file,path]=uigetfile(filterspec,'Select File');
    if(isequal(0,file))
        error('seizmo:readtxt:noFileSelected','No input file selected!');
    end
    file=[path fs file];
else
    % check file
    if(~isstring(file))
        error('seizmo:readtxt:fileNotString',...
            'FILE must be a string!');
    end
    if(~isabspath(file)); file=[pwd fs file]; end
    if(~exist(file,'file'))
        error('seizmo:readtxt:fileDoesNotExist',...
            'File: %s\nDoes Not Exist!',file);
    elseif(exist(file,'dir'))
        error('seizmo:readtxt:dirConflict',...
            'File: %s\nIs A Directory!',file);
    end
end

% open file for reading
fid=fopen(file);

% check if file is openable
if(fid<0)
    error('seizmo:readtxt:cannotOpenFile',...
        'File: %s\nNot Openable!',file);
end

% read in file and close
txt=fread(fid,'*char');
fclose(fid);

% row vector
txt=txt';

end
