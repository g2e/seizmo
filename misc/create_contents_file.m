function []=create_contents_file(mdir,desc,file,o)
%CREATE_CONTENTS_FILE    Create Contents.m file for a directory of m-files
%
%    Usage:    create_contents_file(mdir,desc,file,overwrite)
%
%    Description:
%     CREATE_CONTENTS_FILE(MDIR,DESC,FILE,OVERWRITE) allows for quick
%     creation of a Contents.m file for a directory containing some mfiles.
%     All inputs are optional (graphical prompts are presented to the user
%     for the directory and output file selection).  MDIR is the directory
%     of mfiles to list.  DESC is a string describing the directory of
%     mfiles.  FILE is the output Contents.m file.  OVERWRITE is a logical
%     flag indicating if over writing an existing file can be done
%     silently.  The Contents.m file is created using the H1 help lines of
%     the mfiles (using function H1).
%
%    Notes:
%
%    Examples:
%     % Giving no inputs will prompt the user for both the directory and
%     % the output file.  The description is generic, making this the
%     % quickest ways to create a Contents.m file:
%     create_contents_file;
%
%    See also: HELP, H1, LOOKFOR, CLEAN_CONTENTS_FILE

%     Version History:
%        Jan.  3, 2011 - initial version
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 23:00 GMT

% todo:

% check nargin
error(nargchk(0,4,nargin));

% directory separator
fs=filesep;

% check directory if given
if(nargin>0 && ~isempty(mdir))
    % check directory
    if(~isstring(mdir))
        error('seizmo:create_contents_file:dirNotString',...
            'MDIR must be a string!');
    end
    if(~isabspath(mdir)); mdir=[pwd fs mdir]; end
    if(~exist(mdir,'dir'))
        error('seizmo:create_contents_file:dirDoesNotExist',...
            'Directory: %s\nDoes Not Exist!',mdir);
    end
end

% check description if given
if(nargin>1 && ~isempty(desc) && ~isstring(desc))
    error('seizmo:create_contents_file:badInput',...
        'DESC must be a string!');
end

% default overwrite to false
if(nargin<4 || isempty(o)); o=false; end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:create_contents_file:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check file if given
if(nargin>2 && ~isempty(file))
    % check file
    if(~isstring(file))
        error('seizmo:create_contents_file:fileNotString',...
            'FILE must be a string!');
    end
    if(exist(file,'file'))
        if(exist(file,'dir'))
            error('seizmo:create_contents_file:dirConflict',...
                'Contents File: %s\nIs A Directory!',file);
        end
        if(~o)
            fprintf('Contents File: %s\nFile Exists!\n',file);
            reply=input('Overwrite? Y/N [N]: ','s');
            if(isempty(reply) || ~strncmpi(reply,'y',1))
                disp('Not overwriting!');
                return;
            end
            disp('Overwriting!');
        end
    end
end

% graphical directory selection
if(nargin<1 || isempty(mdir))
    mdir=uigetdir([],'Select a directory');
    if(isequal(0,mdir))
        error('seizmo:create_contents_file:noDirSelected',...
            'No input directory selected!');
    end
end

% default description
if(nargin<2 || isempty(desc)); desc=['Contents of ' mdir ':']; end

% graphical file selections
if(nargin<3 || isempty(file))
    [file,path]=uiputfile(...
        {'*.m' 'M Files (*.m)';
        '*.*' 'All Files (*.*)'},...
        'Save M File as',[mdir fs 'Contents.m']);
    if(isequal(0,file))
        error('seizmo:create_contents_file:noFileSelected',...
            'No output file selected!');
    end
    file=strcat(path,fs,file);
end

% get the mfiles
list=onefilelist([mdir fs '*.m']);

% preallocate output array
nfiles=numel(list);
tmp=cell(nfiles+1,1);
tmp{1}=['% ' desc sprintf('\n')];

% setup individual function lines
for a=1:nfiles
    % retrieve h1 line
    tmp{a+1}=h1([list(a).path fs list(a).name]);
    
    % add comment char & line terminators
    tmp{a+1}=['%' tmp{a+1} sprintf('\n')];
end

% open file for writing
fid=fopen(file,'w');

% check if file is openable
if(fid<0)
    error('seizmo:create_contents_file:cannotOpenFile',...
        'Contents File: %s\nNot Openable!',file);
end

% write to file
for a=1:nfiles+1
    cnt=fwrite(fid,tmp{a},'char');
    
    % check count
    if(cnt~=numel(tmp{a}))
        error('seizmo:create_contents_file:writeFailed',...
            'Contents File: %s\nCould not write (some) text!',file);
    end
end

% close file
fclose(fid);

end
