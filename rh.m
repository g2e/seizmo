function [data]=rh(varargin)
%RH    Read SAClab headers from binary datafiles
%
%    Description: Reads headers from binary datafiles into a SAClab
%     structure.  Accepts character arrays of filenames (one filename per 
%     row) and/or cell arrays of filenames (one filename per cell).
%     Wildcards are valid.
%
%     Fields of output SAClab structure:
%      head - contains header data
%      name - filename (may include path)
%      endian - byte-order of file (ieee-le or ieee-be)
%      version - version of file
%
%    Usage:    SAClab_struct=rh(['seisfile1'; 'seisfile2'; ...],...
%                                 {'seisfile3' 'seisfile4' ...},...
%                                 'seisfile5','seisfile6',...)
%
%    Examples:
%     Some basic examples:
%      data=rh('KATH.R');
%      data=rh('SQRL.R','AAK.R');
%
%     Maybe you have several big filelists to read in:
%      data=rh(cellarray1,chrarray1,chrarray2,cellarray2);
%
%     Read in headers of all SAClab readible files in current directory:
%      data=rh('*')
%
%     Read in some datafile's headers, modify them, and write out changes:
%      wh(ch(rh('*.SAC'),'kuser0','QCed'))
%
%    See also: rdata, rpdw, rseis, wseis, wh, bseis, seisdef, gv

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - Documentation Update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2008 at 05:50 GMT

% compile file lists
varargin=onelist(varargin{:});
nfiles=length(varargin);

% empty data if no files
if(nfiles<1); nfiles=[]; end

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
            'File not openable, %s',varargin{i});
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
            'File too short, %s',varargin{i});
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
