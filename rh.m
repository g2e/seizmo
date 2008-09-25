function [data]=rh(varargin)
%RH    Read SAClab headers from datafiles
%
%    Description: RH(FILELIST1,...,FILELISTN) reads the headers of SAClab
%     compatible datafiles into a SAClab structure.  Accepts character 
%     arrays of filenames (one filename per row) and/or cell arrays of 
%     filenames (one filename per cell).  Wildcards are valid.
%
%     Fields of output SAClab structure:
%      head - contains header data
%      name - filename (may include path)
%      endian - byte-order of file (ieee-le or ieee-be)
%      version - version of datafile
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: arguments must be char or cellstr arrays
%
%    Header changes: NONE
%
%    Usage:    data=rh(filelist1,...,filelistN)
%
%    Examples:
%     Some basic examples:
%      data=rh('KATH.R');
%      data=rh('SQRL.R','AAK.R');
%
%     Maybe you have several big filelists to read in:
%      data=rh(USARRAYdata,FLEDdata,MOMAdata);
%
%     Read in headers of all SAClab readible files in current directory:
%      data=rh('*')
%
%     Read in some datafile's headers, modify them, and write out changes:
%      wh(ch(rh('*.SAC'),'kuser0','QCed'))
%
%    See also: rdata, rpdw, rseis, wseis, wh, bseis, seisdef, gv

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 18, 2008 - fixed argument parse bug
%        Feb. 29, 2008 - doc update, minor code cleaning
%        Mar.  2, 2008 - code cleaning
%        Mar.  4, 2008 - minor code cleaning
%        Apr. 23, 2008 - minor doc update (wildcards now valid)
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update, error on no files
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2008 at 23:20 GMT

% compile file lists
varargin=onefilelist(varargin{:});
nfiles=length(varargin);

% error if no files
if(nfiles<1)
    error('SAClab:rh:noFilesToRead','No files to read!');
end

% pre-allocating SAClab structure
data(nfiles,1)=struct('name',[],'head',[],'endian',[],'version',[]);

% loop for each file
destroy=false(nfiles,1);
for i=1:nfiles
    % get version and byte-order
    [version,endian]=gv(varargin{i});
    
    % validity check
    if(version==0)
        % invalid - jump to next file
        destroy(i)=true;
        continue;
    end
    
    % open file for reading
    fid=fopen(varargin{i},'r',endian);
    
    % fid check
    if(fid<0)
        % non-existent file or directory
        warning('SAClab:rh:badFID',...
            'File not openable, %s !',varargin{i});
        destroy(i)=true;
        continue;
    end
    
    % retrieve header setup
    h=seisdef(version);
    
    % get file size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    
    % minimum size check
    if(bytes<h.data.startbyte)
        % smaller than header size
        fclose(fid);
        warning('SAClab:rh:fileTooShort',...
            'File too short, %s !',varargin{i});
        destroy(i)=true;
        continue;
    end
    
    % save filename, byte-order, version
    data(i).name=varargin{i};
    data(i).endian=endian;
    data(i).version=version;
    
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

% remove unread entries
data(destroy)=[];

end
