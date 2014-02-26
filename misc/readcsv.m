function [s]=readcsv(file,delimiter,flag)
%READCSV    Read in .csv formatted file as a structure
%
%    Usage:    s=readcsv(file)
%              s=readcsv(file,delimiter)
%              s=readcsv(string,delimiter,true)
%
%    Description:
%     S=READCSV(FILE) reads in a comma-separated values (CSV) text file
%     FILE as a scalar struct S.  FILE is expected to have 1 header line
%     of comma-separated field names.  All values are treated as text and
%     are not processed so the output struct fields are cellstr arrays of
%     size NLx1 where NL is the number of lines (minus the header line).
%     The function SSIDX is useful for accessing specific entries of scalar
%     structs (see the Examples section below).  Calling READCSV without
%     the FILE input or with FILE set to '' will present a graphical file
%     selection menu.
%
%     STRUCT=READCSV(FILE,DELIMITER) uses the single character DELIMITER to
%     parse the .csv file.  The default delimiter is ',' (comma).  Line
%     termination characters (linefeed & carriage return) are not allowed.
%
%     S=READCSV(STRING,DELIMITER,TRUE) will take the first input to be a
%     character string of the contents of an .csv file (rather than the
%     filename) if the third input is set to the logical TRUE.  The string
%     should be the same as if the file was read with READTXT (a single row
%     char vector with linefeeds included).  This is also useful when
%     combined with URLREAD to directly grab info from the web rather than
%     having to save it to your filesystem first.
%
%    Notes:
%     - Does not handle text entries that contain the delimiter or a line
%       terminator!
%     - White space before/after an entry is not preserved!
%
%    Examples:
%     % Get the first 3 lines of a csv file:
%     ssidx(readcsv,1:3)
%
%    See also: WRITECSV, CSVREAD, CSVWRITE

%     Version History:
%        Sep. 16, 2009 - initial version
%        Jan. 26, 2010 - add graphical selection
%        Feb.  5, 2010 - improved file checks
%        Jan. 12, 2011 - nargchk fix, fix for stricter getwords
%        Jan. 28, 2011 - handle empty entries (,,), allow alternate
%                        field delimiter, require nonempty field name
%        Feb. 28, 2012 - output is scalar struct now
%        Jan. 26, 2014 - abs path exist fix
%        Feb.  8, 2014 - use readtxt, doc update
%        Feb.  9, 2014 - contents as a string input allowed now
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2014 at 17:25 GMT

% todo:
% - text delimiter
% - handle field delimiter in text
%   - find delimiters between text delimiters and remove them from
%     delimiters list
%       - too advanced for getwords :(
% - linefeeds are not allowed to be delimiters

% check nargin
error(nargchk(0,3,nargin));

% default flag
if(nargin<3 || isempty(flag)); flag=false; end

% skip if string input
if(~flag)
    % read in csv file
    if(nargin<1); file=[]; end
    txt=readtxt(file,{'*.csv;*.CSV' 'CSV Files (*.csv,*.CSV)';
    '*.*' 'All Files (*.*)'});
else
    % just copy file to txt
    if(nargin<1 || isempty(file) || ~isstring(file))
        error('seizmo:readcsv:emptyStr',...
            'STRING must be non-empty!');
    else
        txt=file;
    end
end

% default/check delimiter
if(nargin<2 || isempty(delimiter)); delimiter=','; end
if(~isstring(delimiter) || ~isscalar(delimiter))
    error('seizmo:readcsv:badDelimiter',...
        'DELIMITER must be a single character!');
end

% replace last lf/cr with nothing
while(isspace(txt(end)))
    txt(end)=[];
end

% get header
lf=txt==10;
cr=txt==13;
hde=find(lf | cr,1,'first');
hd=txt(1:hde-1);

% breakup header into individual fields
% - treat ',,' as an empty entry
f=strtrim(getwords(hd,delimiter,false));

% require all fields have a name
for i=1:numel(f)
    if(numel(f{i})==0)
        error('seizmo:readcsv:badFieldName',...
            'CSV File: %s\nField %d is not named!',file,i);
    end
end

% replace line terminators with commas
if(any(lf))
    % replace linefeeds with commas
    txt(lf)=delimiter;
    if(any(cr))
        % replace carriage returns with nothing
        txt(cr)=[];
    end
elseif(any(cr))
    % replace carriage returns with commas
    txt(cr)=delimiter;
end

% breakup str into individual values
% - treat ',,' as an empty entry
v=strtrim(getwords(txt,delimiter,false))';

% get number of fields, values, events
nf=numel(f);
nv=numel(v)-nf;
nev=nv/nf;

% just to make sure (not definitive but it makes sure the last part works)
if(nev~=fix(nev))
    error('seizmo:readcsv:malformedEventCSV',...
        ['CSV file: %s\n'...
        'Some lines have an incorrect number of entries!'],file);
end

% create scalar struct
for i=1:nf
    s.(f{i})=v(nf+i:nf:end);
end

end
