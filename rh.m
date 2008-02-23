function [data]=rh(varargin)
%RH    Read SAC binary file headers
%
%    Description: Reads in SAC (seismic analysis code) binary file headers
%     into a SAClab structure.  Accepts character arrays of filenames (one
%     filename per row) and/or cell arrays of filenames (one filename per
%     cell).
%
%     Structure fields for all SAC files:
%      data.head - header
%      data.name - filename (may include path)
%      data.endian - byte-order of file
%      data.version - version of file
%
%    Usage:    saclab_struct=rh(['sacfile1'; 'sacfile2'; ...],...
%                                 {'sacfile3' 'sacfile4' ...},...
%                                 'sacfile5','sacfile6',...)
%
%    Examples:
%     data=rh('KATH.R');
%     data=rh('SQRL.R','AAK.R');
%     data=rh(cellarray1,chrarray1,chrarray2,cellarray2);
%
%     read in headers of SAC files from current directory
%      files=dir();
%      data=rh(files.name);
%
%    See also: rdata, rpdw, rsac, wsac, bsac, sachi, gh, lh, ch, wh, gv

% compile file lists
varargin=onelist(varargin{:});
nfiles=length(varargin);

% empty data if no files
if(nfiles<1); data=[]; return; end

% pre-allocating SAClab structure
data(nfiles,1)=struct('name',[],'head',[],'endian',[],'version',[]);

% initializing count of invalid SAC files
j=0;

% loop for each file
for i=1:nfiles
    % get version and byte-order
    [version,endian]=gv(varargin{i});
    
    % validity check
    if(version==0)
        % invalid - jump to next file
        j=j+1;
        continue;
    end
    
    % open file for reading
    fid=fopen(varargin{i},'r',endian);
    
    % invalid fid check (non-existent file or directory)
    if(fid<0)
        warning('SAClab:badFID','File not openable, %s',varargin{i});
        
        % jump to next file
        j=j+1;
        continue;
    end
    
    % retrieve header setup
    h=sachi(version);
    
    % file size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    
    % minimum size check
    if(bytes<h.data.startbyte)
        % smaller than header size
        fclose(fid);
        warning('SAClab:fileTooShort','File too short, %s',varargin{i});
        
        % jump to next file
        j=j+1;
        continue;
    end
    
    % save filename, byte-order, version
    data(i-j).name=varargin{i};
    data(i-j).endian=endian;
    data(i-j).version=version;
    
    % preallocate header
    data(i-j).head=zeros(h.size,1,h.store);
    
    % reading in header
    n=h.types;
    for m=1:length(n)
        for k=1:length(h.(n{m}))
            fseek(fid,h.(n{m})(k).startbyte,'bof');
            data(i-j).head(h.(n{m})(k).minpos:h.(n{m})(k).maxpos)=...
                fread(fid,h.(n{m})(k).size,h.(n{m})(k).store);
        end
    end
    
    % closing file
    fclose(fid);
end

% remove unused entries
data(nfiles-j+1:nfiles)=[];

end
