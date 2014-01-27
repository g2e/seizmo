function []=writecsv(file,struct,delimiter,overwrite)
%WRITECSV    Write out .csv formatted file from a structure
%
%    Usage:    writecsv(file,struct)
%              writecsv(file,struct,delimiter)
%              writecsv(file,struct,delimiter,overwrite)
%
%    Description:
%     WRITECSV(FILE,STRUCT) writes a comma-separated values (CSV) text file
%     FILE using the scalar struct array STRUCT.  STRUCT is expected to be
%     a single-level structure with all fields containing equally sized
%     cellstr arrays (the values).  This corresponds to the output from
%     READCSV.  If FILE is an empty string a graphical file creation menu
%     is presented.
%
%     WRITECSV(FILE,STRUCT,DELIMITER) uses the character DELIMITER to
%     delimit the values in FILE.  Useful when entries contain commas.
%
%     WRITECSV(FILE,STRUCT,DELIMITER,OVERWRITE) quietly overwrites a pre-
%     existing CSV file without confirmation when OVERWRITE is set to TRUE.
%     By default OVERWRITE is FALSE.  OVERWRITE is ignored in the graphical
%     file creation menu.
%
%    Notes:
%     - Beware of text entries with commas or line terminators.  They
%       will not be read back in correctly!
%
%    Examples:
%     % Read a SOD (Standing Order for Data) generated event csv file,
%     % change the locations a bit, and write out the updated version:
%     events=readcsv('events.csv')
%     [events.latitude,events.longitude]=...
%         deal(events.longitude,events.latitude);
%     writecsv('events_flip.csv',events)
%
%     % Make your own .csv:
%     a=struct('yo',{'woah' 'dude!'},'another',{'awe' 'some'});
%     writecsv('easy.csv',a)
%
%     % Now use the pipe character '|' as the delimiter:
%     writecsv('easypipe.csv',a,'|')
%
%    See also: READCSV, CSVREAD, CSVWRITE

%     Version History:
%        Sep. 16, 2009 - initial version
%        Sep. 22, 2009 - fixed dir check bug, skip confirmation option
%        Jan. 26, 2010 - add graphical selection
%        Feb.  5, 2010 - add check for overwrite flag
%        Feb. 11, 2011 - mass nargchk fix, use fprintf
%        Feb. 28, 2012 - input is scalar struct now, delimiter input
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% directory separator
fs=filesep;

% defaults
if(nargin<3 || isempty(delimiter)); delimiter=', '; end
if(nargin<4 || isempty(overwrite)); overwrite=false; end
if(~ischar(delimiter) || size(delimiter,1)~=1 || ndims(delimiter)>2)
    error('seizmo:writecsv:badInput',...
        'DELIMITER must be a string!');
elseif(~isscalar(overwrite) || ~islogical(overwrite))
    error('seizmo:writecsv:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check structure
if(~isstruct(struct) || ~isscalar(struct))
    error('seizmo:writecsv:badInput',...
        'STRUCT must be a scalar struct!');
end

% graphical selection
if(isempty(file))
    [file,path]=uiputfile(...
        {'*.csv;*.CSV' 'CSV Files (*.csv,*.CSV)';
        '*.*' 'All Files (*.*)'},...
        'Save CSV File as');
    if(isequal(0,file))
        error('seizmo:writecsv:noFileSelected','No output file selected!');
    end
    file=[path fs file];
else
    % check file
    if(~isstring(file))
        error('seizmo:writecsv:fileNotString',...
            'FILE must be a string!');
    end
    if(~isabspath(file)); file=[pwd fs file]; end
    if(exist(file,'file'))
        if(exist(file,'dir'))
            error('seizmo:writecsv:dirConflict',...
                'CSV File: %s\nIs A Directory!',file);
        end
        if(~overwrite)
            fprintf('CSV File: %s\nFile Exists!\n',file);
            reply=input('Overwrite? Y/N [N]: ','s');
            if(isempty(reply) || ~strncmpi(reply,'y',1))
                disp('Not overwriting!');
                return;
            end
            disp('Overwriting!');
        end
    end
end

% get the field names
f=fieldnames(struct);
nf=numel(f);

% build the cellstr array
nlines=numel(struct.(f{1}))+1;
tmp(1,1:2*nf*nlines)={delimiter};
tmp(1:2:2*nf)=f;
for i=1:nf
    if(numel(struct.(f{i}))~=nlines-1)
        error('seizmo:writecsv:badInput',...
            'Values in STRUCT are not the same size!');
    end
    tmp(2*nf+(2*i-1):2*nf:end)=struct.(f{i});
end

% add line terminators
tmp(2*nf:2*nf:end)={sprintf('\n')};

% check is cellstr
if(~iscellstr(tmp))
    error('seizmo:writecsv:badInput',...
        'Not all values in STRUCT are be cellstr!');
end

% build the char vector
tmp=cell2mat(tmp);

% open file for writing
fid=fopen(file,'w');

% check if file is openable
if(fid<0)
    error('seizmo:writecsv:cannotOpenFile',...
        'CSV File: %s\nNot Openable!',file);
end

% write to file
cnt=fwrite(fid,tmp,'char');

% check count
if(cnt~=numel(tmp))
    error('seizmo:writecsv:writeFailed',...
        'CSV File: %s\nCould not write (some) text to file!',file);
end

% close file
fclose(fid);

end
