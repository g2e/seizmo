function [data]=rh(varargin)
%RH    Read binary seismic datafile header
%
%    Description: Reads in binary seismic datafile headers into a seislab
%     structure.  Accepts character arrays of filenames (one filename per 
%     row) and/or cell arrays of filenames (one filename per cell).
%
%     Structure fields of output:
%      head - contains header data
%      name - filename (may include path)
%      endian - byte-order of file (ieee-le or ieee-be)
%      version - version of file
%
%    Usage:    seislab_struct=rh(['seisfile1'; 'seisfile2'; ...],...
%                                 {'seisfile3' 'seisfile4' ...},...
%                                 'seisfile5','seisfile6',...)
%
%    Examples:
%     data=rh('KATH.R');
%     data=rh('SQRL.R','AAK.R');
%     data=rh(cellarray1,chrarray1,chrarray2,cellarray2);
%
%     read in headers of files in current directory
%      files=dir();
%      data=rh(files.name);
%
%    See also: rdata, rpdw, rseis, wseis, bseis, seishi, gv

% compile file lists
varargin=onelist(varargin{:});
nfiles=length(varargin);

% empty data if no files
if(nfiles<1); nfiles=[]; end

% pre-allocating seislab structure
data(nfiles,1)=struct('name',[],'head',[],'endian',[],'version',[]);

% loop for each file
destroy=false(nfiles,1);
for i=1:nfiles
    % get version and byte-order
    [version,endian]=gv(varargin{i});
    
    % validity check
    if(version==0)
        % invalid - jump to next file
        destroy(i)=1;
        continue;
    end
    
    % open file for reading
    fid=fopen(varargin{i},'r',endian);
    
    % fid check
    if(fid<0)
        % non-existent file or directory
        warning('seislab:rh:badFID',...
            'File not openable, %s',varargin{i});
        destroy(i)=1;
        continue;
    end
    
    % retrieve header setup
    h=seishi(version);
    
    % get file size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    
    % minimum size check
    if(bytes<h.data.startbyte)
        % smaller than header size
        fclose(fid);
        warning('seislab:rh:fileTooShort',...
            'File too short, %s',varargin{i});
        destroy(i)=1;
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
