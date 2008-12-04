function []=writeseizmo(data)
%WRITESEIZMO    Write SEIZMO records to datafiles
%
%    Description: WRITESEIZMO(DATA) writes SEIZMO records in DATA to
%     datafiles.  Uses the fields in the structure to determine how and
%     where to write the records.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: see CHECKHEADER
%
%    Usage:    writeseizmo(data)
%
%    Examples:
%     Read in some files, clean them up and write over:
%      writeseizmo(taper(removetrend(readseizmo('*'))))
%
%     To write out all records in big endian:
%      writeseizmo(changebyteorder(data,'ieee-be'))
%
%     Alter the location of where the files are written:
%      [data.path]=deal('some/new/directory');
%      writeseizmo(data)
%
%    See also:  readseizmo, bseizmo, seizmodef, getversion, readdata,
%               readdatawindow, readheader, writeheader

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
%        Nov. 17, 2008 - update for new name schema (now WRITESEIZMO)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 17, 2008 at 18:10 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% headers setup
[h,vi]=versioninfo(data);

% check headers
data=checkheader(data);

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% estimated filesize from header
est_bytes=seizmosize(data);

% header info
ncmp=getncmp(data);
npts=getheader(data,'npts');
iftype=getenumdesc(data,'iftype');
leven=getlgc(data,'leven');

% loop over records
for i=1:length(data)
    % skip writing if dataless
    if(~data(i).hasdata)
        warning('seizmo:writeseizmo:dataless',...
            'Use WRITEHEADER to write only headers! Record %d Skipped.',i);
        continue; 
    end
    
    % construct fullname
    name=fullfile(data(i).path,data(i).name);

    % open file for writing
    fid=fopen(name,'w',data(i).byteorder);
    
    % check fid
    if(fid<0)
        % unopenable file for writing (permissions/directory?)
        error('seizmo:writeseizmo:badFID',...
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
        error('seizmo:writeseizmo:writeFailed',...
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
        error('seizmo:writeseizmo:badFileSize',...
            ['Output file disksize does not match expected size!\n'...
            'Record: %d, File: %s'],i,name);
    end
    
    % close file
    fclose(fid);
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
