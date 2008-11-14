function []=wseis(data)
%WSEIS    Write SAClab data to datafiles
%
%    Description: WSEIS(DATA) writes SAClab records in DATA to datafiles.
%     The 'name' field is used as the output filename, the 'location' field
%     gives the path and the 'endian' field specifies the byte-order.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: Any inconsistent field found with CHKHDR
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
%     Alter the first record's filename and write it:
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
%        Oct.  7, 2008 - combined write code for similar filetypes
%        Oct.  8, 2008 - fixed a couple bugs caused by last change
%        Oct. 15, 2008 - fixed spectral data support, drop support for
%                        blank fields, support hasdata field, better
%                        separation of even vs uneven checks
%        Oct. 16, 2008 - moved checks to CHKHDR, support new struct layout
%        Oct. 27, 2008 - update for struct change, use state changes for
%                        SEISCHK & CHKHDR
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 27, 2008 at 02:50 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seischk(data,'dep'))

% turn off struct checking
oldseischkstate=get_seischk_state;
set_seischk_state(true);

% headers setup
[h,vi]=vinfo(data);

% check headers
data=chkhdr(data);

% turn off header checking
oldchkhdrstate=get_chkhdr_state;
set_chkhdr_state(true);

% estimated filesize from header
est_bytes=seissize(data);

% header info
ncmp=gncmp(data);
npts=gh(data,'npts');
iftype=genumdesc(data,'iftype');
leven=glgc(data,'leven');

% loop over records
for i=1:length(data)
    % skip writing if dataless
    if(~data(i).hasdata)
        warning('SAClab:wseis:dataless',...
            'Use WH to write only headers! Record %d Skipped.',i);
        continue; 
    end
    
    % construct fullname
    name=fullfile(data(i).location,data(i).name);

    % open file for writing
    fid=fopen(name,'w',data(i).endian);
    
    % check fid
    if(fid<0)
        % unopenable file for writing (permissions/directory?)
        error('SAClab:wseis:badFID',...
            ['File not openable for writing!\n'...
            '(Permissions problem / conflict with directory?)\n'...
            'Record: %d, File: %s\n'],i,name);
    end
    
    % fill header with dummy bytes (so we can seek around)
    fseek(fid,0,'bof');
    count=fwrite(fid,zeros(h(vi(i)).data.startbyte,1),'char');
    
    % verify write
    if(count<h(vi(i)).data.startbyte)
        % write failed
        fclose(fid);
        error('SAClab:wseis:writeFailed',...
            'Writing failed!\nRecord: %d, File: %s',i,name);
    end
    
    % write header
    n=h(vi(i)).types;
    for m=1:length(n)
        for k=1:length(h(vi(i)).(n{m}))
            fseek(fid,h(vi(i)).(n{m})(k).startbyte,'bof');
            fwrite(fid,data(i).head(h(vi(i)).(n{m})(k).minpos:...
                h(vi(i)).(n{m})(k).maxpos),h(vi(i)).(n{m})(k).store);
        end
    end
    
    % skip data writing if npts==0
    if(npts(i)==0); fclose(fid); continue; end
    
    % act by file type (new filetypes will have to be added here)
    fseek(fid,h(vi(i)).data.startbyte,'bof');
    if(any(strcmpi(iftype(i),{'Time Series File' 'General X vs Y file'})))
        % dependent component(s) of data
        for k=1:ncmp(i)
            fwrite(fid,data(i).dep(:,k),h(vi(i)).data.store);
        end
        
        % independent component of data if uneven
        if(strcmpi(leven(i),'false'))
            fwrite(fid,data(i).ind(:),h(vi(i)).data.store);
        end
    elseif(any(strcmpi(iftype(i),...
            {'Spectral File-Real/Imag' 'Spectral File-Ampl/Phase'})))
        % spectral file
        for k=1:ncmp(i)
            fwrite(fid,data(i).dep(:,2*k-1),h(vi(i)).data.store);
            fwrite(fid,data(i).dep(:,2*k),h(vi(i)).data.store);
        end
    elseif(any(strcmpi(iftype(i),'General XYZ (3-D) file')))
        % general xyz (3D) grid - nodes are evenly spaced
        for k=1:ncmp(i)
            fwrite(fid,data(i).dep(:,k),h(vi(i)).data.store);
        end
    end
    
    % verify written file's size
    fseek(fid,0,'eof');
    bytes=ftell(fid);
    if(bytes~=est_bytes(i))
        % write failed/incomplete
        fclose(fid);
        error('SAClab:wseis:badFileSize',...
            ['Output file disksize does not match expected size!\n'...
            'Record: %d, File: %s'],i,name);
    end
    
    % close file
    fclose(fid);
end

% toggle checking back
set_seischk_state(oldseischkstate);
set_chkhdr_state(oldchkhdrstate);

end
