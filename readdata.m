function [data,failed]=readdata(data,varargin)
%READDATA    Read SEIZMO data from datafiles
%
%    Usage: data=readdata(data)
%           data=readdata(data,'trim',true|false)
%           [data,failed]=readdata(...)
%
%    Description: [OUTDATA,FAILED]=READDATA(DATA) reads in data from SEIZMO 
%     compatible datafiles utilizing the header info in DATA, and returns 
%     the combined dataset as OUTDATA.  Datafiles that are not SEIZMO
%     compatible or have errors will be removed from the returned dataset.
%     Optional output FAILED returns a logical matrix equal in size to DATA
%     with entries set to TRUE for those records which had reading errors.
%
%     [OUTDATA,FAILED]=READDATA(DATA,'TRIM',LOGICAL) sets option parameter
%     TRIM to determine how to handle data that had errors.  By default 
%     TRIM is set to TRUE, which deletes any records from OUTDATA that had 
%     errors while reading.  Setting TRIM to FALSE will preserve records in
%     OUTDATA that had errors.
%
%     SEIZMO data structure setup:
%
%     Fields for all files:
%      path - directory of file
%      name - file name
%      filetype - type of file
%      version - version of filetype
%      byteorder - byte-order of file (ieee-le or ieee-be)
%      head - header data
%      hasdata - logical indicating if data is read in
%      ind - independent component data (for uneven)
%      dep - dependent component data
%      misc - place for miscellaneous record info
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
%       six columns of data.  All dependent components ('dep') share the 
%       same independent component ('ind').
%
%    Header changes: see CHECKHEADER
%
%    Examples:
%     Read in datafiles (headers only) from the current directory, subset
%     it to include only time series files, and then read in the associated
%     time series data:
%      data=readheader('*')
%      data=data(strcmpi(getenumid(data,'iftype'),'itime'))
%      data=readdata(data)
%
%    See also: READHEADER, READDATAWINDOW, READSEIZMO, WRITESEIZMO, GETLGC,
%              WRITEHEADER, BSEIZMO, SEIZMODEF, GETFILEVERSION, SEIZMOSIZE,
%              GETHEADER, LISTHEADER, CHANGEHEADER, GETENUMID, GETENUMDESC

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
%        Oct.  8, 2008 - fixed a couple bugs from last change
%        Oct. 15, 2008 - update to new SEISCHK, hasdata field
%        Oct. 17, 2008 - added CHKHDR, support new struct layout
%        Oct. 27, 2008 - minor doc update, updated for struct changes,
%                        make use of state changes for SEISCHK & CHKHDR
%        Nov. 12, 2008 - fix allocation of ind cmp for uneven records
%        Nov. 17, 2008 - update for new name schema (now READDATA)
%        Mar.  3, 2009 - update for GETFILEVERSION
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 12, 2009 - error msg for empty data
%        Sep. 11, 2009 - added misc field to doc
%        Sep. 29, 2009 - minor doc update
%        Oct.  5, 2009 - minor doc update
%        Oct. 16, 2009 - drop fullfile usage
%        Jan. 30, 2010 - better SAC bug workaround, versioninfo caching,
%                        proper SEIZMO handling, seizmoverbose support
%        Feb.  3, 2010 - update for XDIR fixes, fix for npts==0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  3, 2010 at 23:10 GMT

% todo:

% check number of inputs
msg=nargchk(1,3,nargin);
if(~isempty(msg)); error(msg); end

% default trim
trim=true;

% legacy trim option
if(nargin==2 && ~isempty(varargin{1}))
    if((~islogical(varargin{1}) && ~isnumeric(varargin{1})) ...
            || ~isscalar(varargin{1}))
        error('seizmo:readdata:badInput',...
            'TRIM option must be logical!');
    else
        trim=varargin{1};
    end
end

% trim option
if(nargin==3)
    if(strcmpi(varargin{1},'trim'))
        if(~iscalar(varargin{2}) || ...
                (~islogical(varargin{2}) && ~isnumeric(varargin{2})))
            error('seizmo:readdata:badInput',...
                'TRIM option must be logical!');
        else
            trim=varargin{2};
        end
    else
        error('seizmo:readdata:badInput','Unknown option!');
    end
end

% balk at empty data structure (no override)
if(isempty(data))
    error('seizmo:readdata:noRecords','No records to read data from!');
end

% get SEIZMO settings
global SEIZMO

% headers setup (also checks struct)
[h,vi]=versioninfo(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
oldversioninfocache=versioninfo_cache(true);

% attempt reading in data
try
    % estimated filesize from header
    [est_bytes,hbytes,dbytes]=seizmosize(data);

    % header info
    [npts,ncmp]=getheader(data,'npts','ncmp');
    iftype=getenumid(data,'iftype');
    leven=getlgc(data,'leven');
    
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % detail message
    if(verbose)
        disp('Reading in Data Portion of Record(s)');
        print_time_left(0,nrecs);
    end

    % read loop
    failed=false(nrecs,1);
    for i=1:nrecs
        % construct fullname
        name=[data(i).path data(i).name];

        % open file for reading
        fid=fopen(name,'r',data(i).byteorder);

        % fid check
        if(fid<0)
            % non-existent file or directory
            warning('seizmo:readdata:badFID',...
                'Record: %d, File not openable: %s !',i,name);
            failed(i)=true;
            % detail message
            if(verbose); print_time_left(i,nrecs,true); end
            continue;
        end

        % file size
        fseek(fid,0,'eof');
        bytes=ftell(fid);

        % byte size check
        if(bytes>est_bytes(i))
            if(bytes==(hbytes(i)+2*dbytes(i)))
                % SAC BUG: converting a spectral record to a time series
                %      record and then writing it out, writes out a 2nd
                %      component that contains gibberish.
            else
                % size big enough but inconsistent - skip
                warning('seizmo:readdata:badFileSize',...
                    ['Record: %d\nFile: %s\n'...
                    'Filesize does not match header info!\n'...
                    '%d (estimated) > %d (on disk) --> Skipping!'],...
                    i,name,est_bytes(i),bytes);
                failed(i)=true;
                % detail message
                if(verbose); print_time_left(i,nrecs,true); end
                continue;
            end
        elseif(bytes<est_bytes(i))
            % size too small - skip
            fclose(fid);
            warning('seizmo:readdata:badFileSize',...
                ['Record: %d\nFile: %s\n'...
                'Filesize does not match header info!\n'...
                '%d (estimated) < %d (on disk) --> Skipping!'],...
                i,name,est_bytes(i),bytes);
            failed(i)=true;
            % detail message
            if(verbose); print_time_left(i,nrecs,true); end
            continue;
        end

        % preallocate data record with NaNs
        % deallocate timing if even sampling
        data(i).dep=nan(npts(i),ncmp(i),h(vi(i)).data.store);
        if(isfield(data,'ind'))
            if(strcmpi(leven(i),'true')); data(i).ind=[];
            else data(i).ind=nan(npts(i),1,h(vi(i)).data.store);
            end
        end

        % skip if npts==0
        if(npts(i)==0)
            % all data read in
            data(i).hasdata=true;
            % closing file
            fclose(fid);
            % detail message
            if(verbose); print_time_left(i,nrecs); end
            continue;
        end

        % act by file type (any new filetypes will have to be added here)
        fseek(fid,h(vi(i)).data.startbyte,'bof');
        if(any(strcmpi(iftype(i),{'itime' 'ixy'})))
            % time series file - amplitude and time
            for k=1:ncmp(i)
                data(i).dep(:,k)=fread(fid,npts(i),...
                    ['*' h(vi(i)).data.store]);
            end

            % timing of amp data if uneven
            if(strcmpi(leven(i),'false'))
                data(i).ind(:,1)=fread(fid,npts(i),...
                    ['*' h(vi(i)).data.store]);
            end
        elseif(any(strcmpi(iftype(i),{'irlim' 'iamph'})))
            % preallocate data record with NaNs
            data(i).dep=nan(npts(i),2*ncmp(i),h(vi(i)).data.store);

            % spectral file - real and imaginary
            for k=1:ncmp(i)
                data(i).dep(:,2*k-1)=fread(fid,npts(i),...
                    ['*' h(vi(i)).data.store]);
                data(i).dep(:,2*k)=fread(fid,npts(i),...
                    ['*' h(vi(i)).data.store]);
            end
        elseif(any(strcmpi(iftype(i),'ixyz')))
            % general xyz (3D) grid - nodes are evenly spaced
            for k=1:ncmp(i)
                data(i).dep(:,k)=fread(fid,npts(i),...
                    ['*' h(vi(i)).data.store]);
            end
        end

        % all data read in
        data(i).hasdata=true;

        % closing file
        fclose(fid);
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end

    % remove unread entries
    if(trim)
        data(failed)=[];
        SEIZMO.VERSIONINFO.IDX(failed)=[];
    end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
catch
    % toggle struct checking back
    seizmocheck_state(oldseizmocheckstate);
    versioninfo_cache(oldversioninfocache);
    
    % rethrow error
    error(lasterror)
end

end
