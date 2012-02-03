function [data]=readheader(varargin)
%READHEADER    Read header info from SEIZMO datafiles
%
%    Usage:    data=readheader(filelist1,...,filelistN)
%
%    Description:
%     READHEADER(FILELIST1,...,FILELISTN) reads the headers of SEIZMO
%     compatible datafiles into a SEIZMO structure.  Accepts char arrays of
%     filenames (one filename per row) and/or cell arrays of filenames (one
%     filename per cell).  Wildcards are valid.
%
%     Fields of output SEIZMO structure:
%      path      - directory of file
%      name      - file name
%      filetype  - type of file
%      version   - version of filetype
%      byteorder - byte-order of file (ieee-le or ieee-be)
%      head      - header data
%      hasdata   - logical indicating if data is read in (false here)
%      ind       - independent component data (for uneven)
%      dep       - dependent component data
%      misc      - place for miscellaneous record info
%
%    Notes:
%
%    Header changes: NONE
%
%    Examples:
%     % Some basic examples:
%     data=readheader('KATH.R');
%     data=readheader('SQRL.R','AAK.R');
%
%     % Maybe you have several big filelists to read in:
%     data=readheader(USARRAYdatafiles,FLEDdatafiles,MOMAdatafiles);
%
%     % Read in headers of all SEIZMO readible files in current directory:
%     data=readheader('*')
%
%     % Read some datafile's headers, modify them, and write out changes:
%     writeheader(changeheader(readheader('*.SAC'),'kuser0','QCed'))
%
%    See also: READDATA, READDATAWINDOW, WRITEHEADER, CHANGEHEADER
%              GETHEADER, LISTHEADER, READSEIZMO, WRITESEIZMO, BSEIZMO

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 18, 2008 - fixed argument parse bug
%        Feb. 29, 2008 - doc update, minor code cleaning
%        Mar.  2, 2008 - code cleaning
%        Mar.  4, 2008 - minor code cleaning
%        Apr. 23, 2008 - minor doc update (wildcards now valid)
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update, error on no files
%        Oct. 15, 2008 - added hasdata field
%        Oct. 17, 2008 - added dir and filetype fields
%        Oct. 27, 2008 - update for struct changes
%        Nov. 15, 2008 - new name schema (now READHEADER)
%        Mar.  3, 2009 - update for GETFILEVERSION
%        Apr. 23, 2009 - move usage up
%        Sep. 11, 2009 - added misc field
%        Sep. 29, 2009 - added dep & ind (to allow dataset concatenation)
%        Oct.  5, 2009 - reordered struct fields, minor doc update
%        Oct. 16, 2009 - dropped fullfile usage, took filesep out of loop
%                        (for speed), reduced seizmodef calls
%        Jan. 30, 2010 - seizmoverbose support, update for XDIR fixes
%        Feb. 14, 2010 - use ONEFILELIST filterspec global option
%        Jan. 30, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 30, 2012 at 16:05 GMT

% todo:

% set filterspec appropriately
global SEIZMO
SEIZMO.ONEFILELIST.FILTERSPEC=...
    {'*.sac;*.SAC;sac.*;SAC.*' 'SAC Files (*.sac,*.SAC,sac.*,SAC.*)';
     '*.sz;*.SZ;sz.*;SZ.*' 'SEIZMO Files (*.sz,*.SZ,sz.*,SZ.*)'};

% compile file lists
varargin=onefilelist(varargin{:});
nfiles=numel(varargin);

% error if no files
if(nfiles<1)
    error('seizmo:readheader:noFilesToRead','No files to read!');
end

% pre-allocating SEIZMO structure
data(nfiles,1)=struct('path',[],'name',[],'filetype',[],...
    'version',[],'byteorder',[],'head',[],'hasdata',[],...
    'ind',[],'dep',[],'misc',[]);

% "preallocate" definition info
filetypes(1:1,1)={''};
versions=zeros(1,1);
def=cell(1,1); c=1;

% verbosity
verbose=seizmoverbose;

% detail message
if(verbose)
    disp('Reading in Header Portion of Record(s)');
    print_time_left(0,nfiles);
end

% loop for each file
destroy=false(nfiles,1);
for i=1:nfiles
    % get filetype, version and byte-order
    name=[varargin(i).path varargin(i).name];
    [filetype,version,endian]=getfileversion(name);
    
    % validity check
    if(isempty(filetype))
        % invalid - jump to next file
        destroy(i)=true;
        % detail message
        if(verbose); print_time_left(i,nfiles,true); end
        continue;
    end
    
    % retrieve header setup
    idx=strcmpi(filetype,filetypes) & version==versions;
    if(any(idx))
        h=def{idx};
    else
        % get definition
        filetypes{c}=filetype;
        versions(c)=version;
        def{c}=seizmodef(filetype,version);
        h=def{c};
        c=c+1;
    end
    
    % handle types separately
    if(strcmpi(filetype,'SAC Binary')...
            || strcmpi(filetype,'SEIZMO Binary'))
        % open file for reading
        fid=fopen(name,'r',endian);
        
        % fid check
        if(fid<0)
            % non-existent file or directory
            warning('seizmo:readheader:badFID',...
                'File not openable, %s !',name);
            destroy(i)=true;
            % detail message
            if(verbose); print_time_left(i,nfiles,true); end
            continue;
        end
        
        % get file size
        fseek(fid,0,'eof');
        bytes=ftell(fid);
        
        % minimum size check
        if(bytes<h.data.startbyte)
            % smaller than header size
            fclose(fid);
            warning('seizmo:readheader:fileTooShort',...
                'File too short, %s !',name);
            destroy(i)=true;
            % detail message
            if(verbose); print_time_left(i,nfiles,true); end
            continue;
        end
        
        % save directory and filename
        data(i).path=varargin(i).path;
        data(i).name=varargin(i).name;
        
        % save filetype, version, endian and hasdata
        data(i).filetype=filetype;
        data(i).version=version;
        data(i).byteorder=endian;
        data(i).hasdata=false;
        
        % preallocate header
        data(i).head=nan(h.size,1,h.store);
        
        % reading in header
        n=h.types;
        for m=1:numel(n)
            for k=1:numel(h.(n{m}))
                fseek(fid,h.(n{m})(k).startbyte,'bof');
                data(i).head(h.(n{m})(k).minpos:h.(n{m})(k).maxpos)=...
                    fread(fid,h.(n{m})(k).size,h.(n{m})(k).store);
            end
        end
        
        % closing file
        fclose(fid);
    end
    
    % detail message
    if(verbose); print_time_left(i,nfiles); end
end

% remove unread entries
data(destroy)=[];

end
