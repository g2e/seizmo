function [data,failed]=rdata(data,varargin)
%RDATA    Read SAClab data from datafiles
%
%    Description: [OUTDATA,FAILED]=RDATA(DATA) reads in data from SAClab 
%     compatible datafiles utilizing the header info in DATA, and returns 
%     the combined dataset as OUTDATA.  Datafiles that are not SAClab
%     compatible or have errors will be removed from the returned dataset.
%     Optional output FAILED returns a logical matrix equal in size to DATA
%     with entries set to TRUE for those records which had reading errors.
%
%     [OUTDATA,FAILED]=RDATA(DATA,'TRIM',LOGICAL) sets optional parameter
%     TRIM to determine how RDATA handles data that had errors. By default 
%     TRIM is set to TRUE, which deletes any records from OUTDATA that had 
%     errors while reading.  Setting TRIM to FALSE will preserve records in
%     OUTDATA that had errors.
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
%      dep(:,1) - amplitudes
%      ind(:,1) - times (if uneven spacing)
%
%     Fields for spectral amp/phase files:
%      dep(:,1) - spectral amplitudes
%      dep(:,2) - spectral phase
%
%     Fields for spectral real/imag files:
%      dep(:,1) - spectral real
%      dep(:,2) - spectral imaginary
%
%     Fields for general xy files:
%      dep(:,1) - dependent component
%      ind(:,1) - independent component (if uneven spacing)
%
%     Fields for xyz grid files:
%      dep(:,1) - matrix data (nodes evenly spaced; advances l2r,b2t)
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
%    Header changes: NONE
%
%    Usage: data=rdata(data)
%           data=rdata(data,'trim',true|false)
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
%        Mar.  4, 2008 - doc update
%        June 12, 2008 - doc update
%        June 23, 2008 - doc update
%        Sep. 15, 2008 - minor doc update, negative NPTS check, enforce
%                        LEVEN=TRUE for spectral and xyz files, .dep and
%                        .ind rather than .x and .t, trim option made to
%                        match RPDW and CUTIM
%        Sep. 26, 2008 - doc update
%        Sep. 27, 2008 - updated for GET_N_CHECK and VINFO
%        Oct.  7, 2008 - combined read code for similar filetypes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  7, 2008 at 03:05 GMT

% todo:

% check number of inputs
error(nargchk(1,3,nargin))

% check data structure
error(seischk(data,'name','endian'))

% default trim
trim=true;

% legacy trim option
if(nargin==2 && ~isempty(varargin{1}))
    if((~islogical(varargin{1}) && ~isnumeric(varargin{1})) ...
            || ~isscalar(varargin{1}))
        error('SAClab:rdata:badInput',...
            'TRIM option not able to be evaluated!');
    else
        trim=varargin{1};
    end
end

% trim option
if(nargin==3)
    if(strcmpi(varargin{1},'trim'))
        if(~iscalar(varargin{2}) || ...
                (~islogical(varargin{2}) && ~isnumeric(varargin{2})))
            error('SAClab:rdata:badInput',...
                'TRIM option unable to be evaluated!');
        else
            trim=varargin{2};
        end
    else
        error('SAClab:rdata:badInput','Bad option!');
    end
end

% number of records
nrecs=numel(data);

% header info
ncmp=gncmp(data);
npts=gh(data,'npts');
iftype=genumdesc(data,'iftype');
leven=glgc(data,'leven');
error(lgcchk('leven',leven(npts>1)))

% estimated filesize from header
est_bytes=seissize(data);

% headers setup
[h,vi]=vinfo(data);

% read loop
failed=false(nrecs,1);
for i=1:nrecs
    % open file for reading
    fid=fopen(data(i).name,'r',data(i).endian);
    
    % fid check
    if(fid<0)
        % non-existent file or directory
        warning('SAClab:rdata:badFID',...
            'File not openable: %s !',data(i).name);
        failed(i)=true;
        continue;
    end
    
    % file size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    
    % byte size check
    if(bytes>est_bytes(i))
        % size big enough but inconsistent - read anyways (SAC bugfix)
        % SAC BUG: converting a spectral file to a time series file does
        % not deallocate the second component, thus the written file has
        % twice as much data.
        warning('SAClab:rdata:badFileSize',...
            ['Filesize of file %s does not match header info!\n'...
            '%d (estimated) > %d (on disk) --> Reading Anyways!'...
            'This is usually caused by a SAC bug and can be ignored.'],...
            data(i).name,est_bytes(i),bytes);
    elseif(bytes<est_bytes(i))
        % size too small - skip
        fclose(fid);
        warning('SAClab:rdata:badFileSize',...
            ['Filesize of file %s does not match header info!\n'...
            '%d (estimated) < %d (on disk) --> Skipping!'],...
            data(i).name,est_bytes(i),bytes);
        failed(i)=true;
        continue;
    end
    
    % preallocate data record with NaNs, deallocate timing
    data(i).dep=nan(npts(i),ncmp(i),h(vi(i)).data.store); 
    data(i).ind=[];
    
    % skip if npts==0 (dataless)
    if(npts(i)==0); fclose(fid); continue; end
    
    % skip if npts<0 (bad)
    if(npts(i)<0)
        fclose(fid);
        warning('SAClab:rdata:nptsBad',...
            'NPTS for file %s can not be set negative!',data(i).name);
        failed(i)=true;
        continue;
    end
    
    % act by file type (any new filetypes will have to be added here)
    fseek(fid,h(vi(i)).data.startbyte,'bof');
    if(strcmpi(iftype(i),{'Time Series File' 'General X vs Y file'}))
        % time series file - amplitude and time
        for k=1:ncmp(i)
            data(i).dep(:,k)=fread(fid,npts(i),['*' h(vi(i)).data.store]);
        end
        
        % timing of amp data if uneven
        if(strcmpi(leven(i),'false'))
            data(i).ind(:,1)=fread(fid,npts(i),['*' h(vi(i)).data.store]);
        end
    elseif(strcmpi(iftype(i),{'Spectral File-Real/Imag' 'Spectral File-Ampl/Phase'}))
        % preallocate data record with NaNs
        data(i).dep=nan(npts(i),2*ncmp(i),h(vi(i)).data.store);
        
        % spectral file - real and imaginary
        for k=1:ncmp(i)
            data(i).dep(:,2*k-1)=fread(fid,npts(i),['*' h(vi(i)).data.store]);
            data(i).dep(:,2*k)=fread(fid,npts(i),['*' h(vi(i)).data.store]);
        end
        
        % check leven
        if(strcmpi(leven(i),'false'))
            fclose(fid);
            warning('SAClab:rh:badLeven',...
                'LEVEN for Spectral file %s must be TRUE!',data(i).name);
            failed(i)=true;
            continue;
        end
    elseif(strcmpi(iftype(i),'General XYZ (3-D) file'))
        % general xyz (3D) grid - nodes are evenly spaced
        for k=1:ncmp(i)
            data(i).dep(:,k)=fread(fid,npts(i),['*' h(vi(i)).data.store]);
        end
        
        % check leven
        if(strcmpi(leven(i),'false'))
            fclose(fid);
            warning('SAClab:rh:badLeven',...
                'LEVEN for XYZ file %s must be TRUE!',data(i).name);
            failed(i)=true;
            continue;
        end
    else
        % unknown filetype
        fclose(fid)
        warning('SAClab:rdata:iftypeBad',...
            'File: %s\nBad filetype: %s !',data(i).name,iftype(i));
        failed(i)=true;
        continue;
    end
    
    % closing file
    fclose(fid);
end

% remove unread entries
if(trim); data(failed)=[]; end

end
