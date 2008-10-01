function []=wseis(data)
%WSEIS    Write SAClab data to datafiles
%
%    Description: Writes SAClab data to datafiles.  For files that do not 
%     already exist the 'name' field is used as the output filename.  To 
%     write files to a specific directory, either add the path to the 
%     'name' field or browse there using the CD command before using WSEIS.  
%     Output file byte-order can be set using the 'endian' field.
%
%    Notes:
%     - unnamed data will be given the name SAClab.N.sac where N is the
%       index number of the data in the SAClab structure
%     - data without an endianness given will be set to that of the current
%       architecture
%
%    System requirements: Matlab 7
%
%    Header changes: NONE (maybe NPTS,NCMP,NVHDR)
%
%    Usage:    wseis(data)
%
%    Examples:
%     Read in some files, clean them up and write over:
%      wseis(taper(rtrend(rseis('*'))))
%
%     To write out all records in big endian:
%      wseis(cendian(data,'ieee-be'))
%
%     Alter the first records path/filename and write it:
%      data(1).name='../../myfile';
%      wseis(data(1))
%
%    See also:  rseis, bseis, seisdef, gv, rpdw, rdata, rh, wh

%     Version History:
%        Oct. 29, 2007 - initial version, supports struct data
%        Nov.  7, 2007 - doc update
%        Jan. 27, 2008 - SACHP support
%        Feb. 11, 2008 - SACHI support
%        Feb. 29, 2008 - minor doc fix
%        Mar.  4, 2008 - renamed from WSAC to WSEIS, better error messages
%        June 12, 2008 - added examples, doc update
%        Sep. 26, 2008 - .dep and .ind rather than .x and .t, history fix,
%                        doc update, NCMP/NPTS/NVHDR checks, error msg
%                        fixes, dataless support
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 26, 2008 at 19:05 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'name','endian'))

% header info
[ncmp,npts,iftype,leven]=get_n_check(data);

% estimated filesize from header
est_bytes=seissize(data);

% headers setup
[h,vi]=vinfo(data);

% platform's endianness
endian=nativeendian;

% loop over records
for i=1:length(data)
    % check for empty filename
    if(isempty(data(i).name))
        warning('SAClab:wseis:namelessData',...
            ['Record %d has no associated filename!\n'...
            '==> Set name as SAClab.%d.sac !'],i,i);
        data(i).name=['SAClab.' num2str(i) '.sac'];
    end
    
    % check for empty endianness
    if(isempty(data(i).name))
        warning('SAClab:wseis:endianlessData',...
            ['Record %d has no associated byteorder!\n'...
            '==> Set byteorder as %s to match platform!'],i,endian);
        data(i).endian=endian;
    end
    
    % check that data matches npts, ncmp
    [nrows,ncols]=size(data(i).dep);
    if(nrows~=npts(i))
        warning('SAClab:wseis:nptsInconsistent',...
            ['NPTS does not match data for record %d !\n'...
             'Changing NPTS to match data! (%d ==> %d)'],...
             i,npts(i),nrows);
        npts(i)=nrows;
        data(i)=ch(data(i),'npts',nrows);
        est_bytes(i)=seissize(data(i));
    end
    if(ncols~=ncmp(i))
        warning('SAClab:wseis:nptsInconsistent',...
            ['NCMP does not match data for record %d !\n'...
             'Changing NCMP to match data! (%d ==> %d)'],...
             i,ncmp(i),ncols);
        ncmp(i)=ncols;
        % check that header can handle multiple components
        if(ncmp(i)>1)
            if(~h(vi(i)).mulcmp.valid)
                % change version
                warning('SAClab:wseis:versNotMulCmp',...
                    ['SAClab version %d cannot handle multiple '...
                    'components!\nChanging record %d to version %d !'],...
                    data(i).version,i,h(vi(i)).mulcmp.altver);
                data(i)=ch(data(i),'nvhdr',h(vi(i)).mulcmp.altver);
                data(i).version=h(vi(i)).mulcmp.altver;
                [h,vi]=vinfo(data);
            end
        end
        data(i)=ch(data(i),'ncmp',ncols);
        est_bytes(i)=seissize(data(i));
    end
    nipts=numel(data(i).ind);
    if(strcmpi(leven(i),'false'))
        if(nipts~=npts(i))
            warning('SAClab:wseis:nptsInconsistent',...
                ['NPTS does not match data for record %d !\n'...
                 'Truncating larger to size of smaller!'],i);
            if(nipts>npts(i))
                data(i).ind=data(i).ind(1:npts(i));
            else
                npts(i)=nipts;
                data(i)=ch(data(i),'npts',nipts);
                data(i).dep=data(i).dep(1:nipts,:);
                est_bytes(i)=seissize(data(i));
            end
        end
    else
        if(nipts~=0)
            warning('SAClab:wseis:unnecessaryData',...
                ['Ignoring unnecessary independent component data\n'...
                 'for evenly-spaced record %d !'],i);
            data(i).ind=[];
        end
    end
    
    % open file for writing
    fid=fopen(data(i).name,'w',data(i).endian);
    
    % check fid
    if(fid<0)
        % unopenable file for writing (permissions/directory?)
        error('SAClab:wseis:badFID',...
            ['File not openable for writing!\n'...
            '(Permissions problem / conflict with directory?)\n'...
            'File: %s\n'],data(i).name);
    end
    
    % fill header with dummy bytes (so we can seek around)
    fseek(fid,0,'bof');
    count=fwrite(fid,zeros(h(vi(i)).data.startbyte,1),'char');
    
    % verify write
    if(count<h(vi(i)).data.startbyte)
        % write failed
        fclose(fid);
        error('SAClab:wseis:writeFailed',...
            'Writing failed!\nFile: %s !',data(i).name);
    end
    
    % write header
    n=h(vi(i)).types;
    for m=1:length(n)
        for k=1:length(h(vi(i)).(n{m}))
            fseek(fid,h(vi(i)).(n{m})(k).startbyte,'bof');
            fwrite(fid,data(i).head(h(vi(i)).(n{m})(k).minpos:h(vi(i)).(n{m})(k).maxpos),h(vi(i)).(n{m})(k).store);
        end
    end
    
    % skip if npts==0 (dataless)
    if(npts(i)==0); fclose(fid); continue; end
    
    % act by file type (new filetypes will have to be added here)
    fseek(fid,h(vi(i)).data.startbyte,'bof');
    if(strcmpi(iftype(i),'Time Series File'))
        % time series file - amplitude and time
        for k=1:ncmp(i)
            fwrite(fid,data(i).dep(:,k),h(vi(i)).data.store);
        end
        
        % timing of amp data if uneven
        if(strcmpi(leven(i),'false'))
            fwrite(fid,data(i).ind(:),h(vi(i)).data.store);
        end
    elseif(strcmpi(iftype(i),'Spectral File-Real/Imag'))
        % spectral file - real and imaginary
        for k=1:ncmp(i)
            fwrite(fid,data(i).dep(:,2*k-1),h(vi(i)).data.store);
            fwrite(fid,data(i).dep(:,2*k),h(vi(i)).data.store);
        end
    elseif(strcmpi(iftype(i),'Spectral File-Ampl/Phase'))
        % spectral file - amplitude and phase
        for k=1:ncmp(i)
            fwrite(fid,data(i).dep(:,2*k-1),h(vi(i)).data.store);
            fwrite(fid,data(i).dep(:,2*k),h(vi(i)).data.store);
        end
    elseif(strcmpi(iftype(i),'General X vs Y file'))
        % general x vs y data (x is 'dependent')
        for k=1:ncmp(i)
            fwrite(fid,data(i).dep(:,k),h(vi(i)).data.store);
        end
        
        % independent data (if uneven)
        if(strcmpi(leven(i),'false'))
            fwrite(fid,data(i).ind(:),h(vi(i)).data.store);
        end
    elseif(strcmpi(iftype(i),'General XYZ (3-D) file'))
        % general xyz (3D) grid - nodes are evenly spaced
        for k=1:ncmp(i)
            fwrite(fid,data(i).dep(:,k),h(vi(i)).data.store);
        end
    else
        % unknown filetype
        fclose(fid)
        error('SAClab:wseis:iftypeBad',...
            'File: %s\nBad filetype: %s !',data(i).name,iftype(i));
    end
    
    % verify written file's size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    if(bytes~=est_bytes(i))
        % write failed/incomplete
        fclose(fid);
        error('SAClab:wseis:badFileSize',...
            ['Output file''s written size does not match expected size!\n'...
            'File: %s'],data(i).name);
    end
    
    % close file
    fclose(fid);
end

end
