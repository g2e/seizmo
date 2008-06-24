function [data,failed]=rdata(data,trim)
%RDATA    Read SAClab data from datafiles on disk
%
%    Description: [OUTDATA,FAILED]=RDATA(DATA,TRIM) reads the data from
%     files on disk utilizing the header info in DATA, and returns the 
%     combined dataset as OUTDATA.  Optional parameter TRIM determines how 
%     RDATA handles data that had errors.  By default TRIM is set to true 
%     which deletes from OUTDATA any records that had errors while reading.  
%     Setting TRIM to FALSE will preserve records in OUTDATA that had 
%     errors.  Optional output FAILED returns a logical matrix equal in 
%     size to DATA with entries set to TRUE for those records which had 
%     reading errors.
%
%     SAClab data structure setup:
%
%     Fields for all files:
%      head - contains header data
%      name - filename (may include path)
%      endian - byte-order of file (ieee-le or ieee-be)
%      version - version of datafile
%
%     Fields for timeseries files:
%      x(:,1) - amplitudes
%      t(:,1) - times (if uneven spacing)
%
%     Fields for spectral amp/phase files:
%      x(:,1) - spectral amplitudes
%      x(:,2) - spectral phase
%
%     Fields for spectral real/imag files:
%      x(:,1) - spectral real
%      x(:,2) - spectral imaginary
%
%     Fields for general xy files:
%      x(:,1) - dependent component
%      t(:,1) - independent component (if uneven spacing)
%
%     Fields for xyz grid files:
%      x(:,1) - matrix data (nodes evenly spaced; advances l2r,b2t)
%     
%    Notes:
%     - Multi-component files will replicate the number of columns by the
%       number of components.  So a three component spectral file will have
%       six columns total in field x.  Components share the same timing.
%     - Currently LEVEN is ignored for spectral and xyz files.  Later
%       versions may require LEVEN to be set to true for these files.
%
%    System requirements: Matlab 7
%
%    Data requirements: DATA has header, endian, name and version fields.
%
%    Header changes: NONE
%
%    Usage: data=rdata(data)
%           data=rdata(data,trim)
%           [data,failed]=rdata(...)
%
%    Examples:
%     Read in datafiles (headers only) from the current directory, subset
%     it to include only time series files, and then read in the associated
%     time series data:
%      data=rh('*')
%      data=data(strcmpi(genumdesc(data,'iftype'),'Time Series File'))
%      data=rdata(data)
%
%    See also: rh, rpdw, rseis, wseis, bseis, seisdef, gv, seissize

%     Version History:
%        Jan. 28, 2008 - initial version
%        Feb. 18, 2008 - works with new GH
%        Feb. 28, 2008 - works with new GLGC, GENUMDESC
%        Feb. 29, 2008 - works with SEISSIZE and new definitions
%        Mar.  3, 2008 - dataless support, workaround for SAC bug
%        Mar.  4, 2008 - documentation update
%        June 12, 2008 - documentation update
%        June 23, 2008 - documentation update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 23, 2008 at 17:40 GMT

% todo:
% .t and .x => .ind and .dep
% leven checks
% stricter check on SAC bug
% type field
% trim option => 'trim',true|false

% check number of inputs
error(nargchk(1,2,nargin))

% check data structure
error(seischk(data,'name','endian'))

% default trim
if(nargin==1 || isempty(trim) || ...
        ~islogical(trim) || ~isscalar(trim))
    trim=true;
end

% number of records
nrecs=length(data);

% estimated filesize from header
est_bytes=seissize(data);

% header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
iftype=genumdesc(data,'iftype');
warning('off','SAClab:gh:fieldInvalid')
[npts,ncmp]=gh(data,'npts','ncmp');
warning('on','SAClab:gh:fieldInvalid')

% clean up and check ncmp
ncmp(isnan(ncmp))=1;
if(any(ncmp<1 | fix(ncmp)~=ncmp))
    error('SAClab:rdata:badNumCmp',...
        'field ncmp must be a positive integer')
end

% headers setup
vers=unique([data.version]);
nver=length(vers);
h(nver)=seisdef(vers(nver));
for i=1:nver-1
    h(i)=seisdef(vers(i));
end

% read loop
failed=false(nrecs,1);
for i=1:nrecs
    % logical index of header info
    v=(data(i).version==vers);
    
    % open file for reading
    fid=fopen(data(i).name,'r',data(i).endian);
    
    % fid check
    if(fid<0)
        % non-existent file or directory
        warning('SAClab:rdata:badFID',...
            'File not openable, %s',data(i).name);
        failed(i)=true;
        continue;
    end
    
    % file size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    
    % byte size check
    if(bytes>est_bytes(i))
        % size big enough but inconsistent - read anyways (SAC bugfix)
        warning('SAClab:rdata:badFileSize',...
            ['Filesize does not match header info (gt)\n'...
            'File: %s'],data(i).name);
    elseif(bytes<est_bytes(i))
        % size too small - skip
        fclose(fid);
        warning('SAClab:rdata:badFileSize',...
            ['Filesize does not match header info (lt)\n'...
            'File: %s'],data(i).name);
        failed(i)=true;
        continue;
    end
    
    % preallocate data record with NaNs, deallocate timing
    data(i).x=nan(npts(i),ncmp(i),h(v).data.store); data(i).t=[];
    
    % skip if npts==0 (dataless)
    if(npts(i)<1); fclose(fid); continue; end
    
    % act by file type (any new filetypes will have to be added here)
    fseek(fid,h(v).data.startbyte,'bof');
    if(strcmpi(iftype(i),'Time Series File'))
        % time series file - amplitude and time
        for k=1:ncmp(i)
            data(i).x(:,k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
        
        % timing of amp data if uneven
        if(strcmpi(leven(i),'false'))
            data(i).t(:,1)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    elseif(strcmpi(iftype(i),'Spectral File-Real/Imag'))
        % preallocate data record with NaNs
        data(i).x=nan(npts(i),2*ncmp(i),h(v).data.store);
        
        % spectral file - real and imaginary
        for k=1:ncmp(i)
            data(i).x(:,2*k-1)=fread(fid,npts(i),['*' h(v).data.store]);
            data(i).x(:,2*k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    elseif(strcmpi(iftype(i),'Spectral File-Ampl/Phase'))
        % preallocate data record with NaNs
        data(i).x=nan(npts(i),2*ncmp(i),h(v).data.store);
        
        % spectral file - amplitude and phase
        for k=1:ncmp(i)
            data(i).x(:,2*k-1)=fread(fid,npts(i),['*' h(v).data.store]);
            data(i).x(:,2*k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    elseif(strcmpi(iftype(i),'General X vs Y file'))
        % general x vs y data (x is 'dependent')
        for k=1:ncmp(i)
            data(i).x(:,k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
        
        % independent data (if uneven)
        if(strcmpi(leven(i),'false'))
            data(i).t(:,1)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    elseif(strcmpi(iftype(i),'General XYZ (3-D) file'))
        % general xyz (3D) grid - nodes are evenly spaced
        for k=1:ncmp(i)
            data(i).x(:,k)=fread(fid,npts(i),['*' h(v).data.store]);
        end
    else
        % unknown filetype
        fclose(fid)
        warning('SAClab:rdata:iftypeBad',...
            'File: %s\nBad filetype: %s',data(i).name,iftype(i));
        failed(i)=true;
        continue;
    end
    
    % closing file
    fclose(fid);
end

% remove unread entries
if(trim); data(failed)=[]; end

end
