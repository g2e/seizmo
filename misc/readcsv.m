function [lines]=readcsv(file,delimiter)
%READCSV    Read in .csv formatted file as a structure
%
%    Usage:    struct=readcsv(file)
%              struct=readcsv(file,delimiter)
%
%    Description:
%     STRUCT=READCSV(FILE) reads in a comma-separated values (CSV) text
%     file FILE as a struct.  FILE is expected to have 1 header line of
%     comma-separated field names.  All values returned in STRUCT are text.
%     Each line in FILE is placed as a separated index in STRUCT, such that
%     all values from line 4 (counting the header line) can be accessed in
%     STRUCT(3).  Calling READCSV without FILE or with FILE set to '' will
%     present a graphical file selection menu.
%
%     STRUCT=READCSV(FILE,DELIMITER) uses the single character DELIMITER to
%     parse the .csv file.  The default delimiter is ',' (comma).  Line
%     termination characters (linefeed & carriage return) are not allowed.
%
%    Notes:
%     - does not handle text entries with the field delimiter or
%       line terminators!
%     - white space before/after an entry is not preserved!
%
%    Examples:
%     Read a SOD (Standing Order for Data) generated event csv file:
%      events=readcsv('events.csv')
%
%    See also: WRITECSV, CSVREAD, CSVWRITE

%     Version History:
%        Sep. 16, 2009 - initial version
%        Jan. 26, 2010 - add graphical selection
%        Feb.  5, 2010 - improved file checks
%        Jan. 12, 2011 - nargchk fix, fix for stricter getwords
%        Jan. 28, 2011 - handle empty entries (,,), allow alternate
%                        field delimiter, require nonempty field name
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2011 at 17:25 GMT

% todo:
% - text delimiter
% - handle field delimiter in text
%   - find delimiters between text delimiters and remove them
%       - too advanced for getwords :(
% - linefeeds are not allowed to be delimiters

% check nargin
error(nargchk(0,2,nargin));

% graphical selection
if(nargin<1 || isempty(file))
    [file,path]=uigetfile(...
        {'*.csv;*.CSV' 'CSV Files (*.csv,*.CSV)';
        '*.*' 'All Files (*.*)'},...
        'Select CSV File');
    if(isequal(0,file))
        error('seizmo:readcsv:noFileSelected',...
            'No input file selected!');
    end
    file=strcat(path,filesep,file);
else
    % check file
    if(~ischar(file))
        error('seizmo:readcsv:fileNotString',...
            'FILE must be a string!');
    end
    if(~exist(file,'file'))
        error('seizmo:readcsv:fileDoesNotExist',...
            'CSV File: %s\nDoes Not Exist!',file);
    elseif(exist(file,'dir'))
        error('seizmo:readcsv:dirConflict',...
            'CSV File: %s\nIs A Directory!',file);
    end
end

% default/check delimiter
if(nargin<2 || isempty(delimiter)); delimiter=','; end
if(~isstring(delimiter) || ~isscalar(delimiter))
    error('seizmo:readcsv:badDelimiter',...
        'DELIMITER must be a single character!');
end

% open file for reading
fid=fopen(file);

% check if file is openable
if(fid<0)
    error('seizmo:readcsv:cannotOpenFile',...
        'CSV File: %s\nNot Openable!',file);
end

% read in file and close
str=fread(fid,'*char').';
fclose(fid);

% replace last lf/cr with nothing
while(isspace(str(end)))
    str(end)=[];
end

% get header
lf=str==10;
cr=str==13;
hde=find(lf | cr,1,'first');
hd=str(1:hde-1);

% breakup header into individual fields
% - treat ',,' as an empty entry
f=strtrim(getwords(hd,delimiter,false));

% require all fields have a name
for i=1:numel(f)
    if(~numel(f{i}))
        error('seizmo:readcsv:badFieldName',...
            'CSV File: %s\nField %d is not named!',file,i);
    end
end

% replace line terminators with commas
if(any(lf))
    % replace linefeeds with commas
    str(lf)=delimiter;
    if(any(cr))
        % replace carriage returns with nothing
        str(cr)=[];
    end
elseif(any(cr))
    % replace carriage returns with commas
    str(cr)=delimiter;
end

% breakup str into individual values
% - treat ',,' as an empty entry
v=strtrim(getwords(str,delimiter,false));

% get number of fields, values, events
nf=numel(f);
nv=numel(v)-nf;
nev=nv/nf;

% just to make sure (not definitive but it makes sure the last part works)
if(nev~=round(nev))
    error('seizmo:readcsv:malformedEventCSV',...
        ['CSV file: %s\n'...
        'Some lines have an incorrect number of entries!'],file);
end

% create struct
lines(1:nev,1)=struct();
for i=1:nf
    [lines.(f{i})]=deal(v{nf+i:nf:end});
end

end
