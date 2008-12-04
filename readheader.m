function [data]=readheader(varargin)
%READHEADER    Read header info from SEIZMO datafiles
%
%    Description: READHEADER(FILELIST1,...,FILELISTN) reads the headers of
%     SEIZMO compatible datafiles into a SEIZMO structure.  Accepts char 
%     arrays of filenames (one filename per row) and/or cell arrays of 
%     filenames (one filename per cell).  Wildcards are valid.
%
%     Fields of output SEIZMO structure:
%      path - directory of file
%      name - file name
%      filetype - type of datafile
%      version - version of filetype
%      byteorder - byte-order of file (ieee-le or ieee-be)
%      hasdata - logical indicating if data is read in (false here)
%      head - contains header data
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: NONE
%
%    Usage:    data=readheader(filelist1,...,filelistN)
%
%    Examples:
%     Some basic examples:
%      data=readheader('KATH.R');
%      data=readheader('SQRL.R','AAK.R');
%
%     Maybe you have several big filelists to read in:
%      data=readheader(USARRAYdatafiles,FLEDdatafiles,MOMAdatafiles);
%
%     Read in headers of all SEIZMO readible files in current directory:
%      data=readheader('*')
%
%     Read in some datafile's headers, modify them, and write out changes:
%      writeheader(changeheader(readheader('*.SAC'),'kuser0','QCed'))
%
%    See also: readdata, readdatawindow, writeheader, changeheader
%              getheader, listheader, readseizmo, writeseizmo, bseizmo

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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 27, 2008 at 03:40 GMT

% todo:

% compile file lists
varargin=onefilelist(varargin{:});
nfiles=numel(varargin);

% error if no files
if(nfiles<1)
    error('seizmo:readheader:noFilesToRead','No files to read!');
end

% pre-allocating SEIZMO structure
data(nfiles,1)=struct('path',[],'name',[],'filetype',[],...
    'version',[],'byteorder',[],'hasdata',[],'head',[]);

% loop for each file
destroy=false(nfiles,1);
for i=1:nfiles
    % get filetype, version and byte-order
    name=fullfile(varargin(i).path,varargin(i).name);
    [filetype,version,endian]=getversion(name);
    
    % validity check
    if(isempty(filetype))
        % invalid - jump to next file
        destroy(i)=true;
        continue;
    end
    
    % retrieve header setup
    h=seizmodef(filetype,version);
    
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
            continue;
        end
        
        % save directory and filename
        if(isempty(varargin(i).path)); varargin(i).path=['.' filesep]; end
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
        for m=1:length(n)
            for k=1:length(h.(n{m}))
                fseek(fid,h.(n{m})(k).startbyte,'bof');
                data(i).head(h.(n{m})(k).minpos:h.(n{m})(k).maxpos)=...
                    fread(fid,h.(n{m})(k).size,h.(n{m})(k).store);
            end
        end
        
        % closing file
        fclose(fid);
    end
end

% remove unread entries
data(destroy)=[];

end
