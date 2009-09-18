function []=writeseizmo(data,varargin)
%WRITESEIZMO    Write SEIZMO records to datafiles
%
%    Usage:    writeseizmo(data)
%              writeseizmo(data,...,'name',name,...)
%              writeseizmo(data,...,'prepend',string,...)
%              writeseizmo(data,...,'append',string,...)
%              writeseizmo(data,...,'delete',string,...)
%              writeseizmo(data,...,'delete',{str1 ... strN},...)
%              writeseizmo(data,...,'change',{original replacement},...)
%              writeseizmo(data,...,'path',path,...)
%              writeseizmo(data,...,'pathprepend',string,...)
%              writeseizmo(data,...,'pathappend',string,...)
%              writeseizmo(data,...,'pathdelete',string,...)
%              writeseizmo(data,...,'pathdelete',{str1 ... strN},...)
%              writeseizmo(data,...,'pathchange',{original replacemnt},...)
%              writeseizmo(data,...,'byteorder',endianness,...)
%
%    Description: WRITESEIZMO(DATA) writes SEIZMO records in DATA to
%     datafiles.  Uses the fields in the structure to determine how and
%     where to write the records.
%
%     For options 'NAME', 'PREPEND', 'APPEND', 'DELETE', and 'CHANGE' see
%     CHANGENAME for usage.  For options 'PATH', 'PATHPREPEND',
%     'PATHAPPEND', 'PATHDELETE', and 'PATHCHANGE' see CHANGEPATH for
%     usage.  For option 'BYTEORDER' see CHANGEBYTEORDER for usage.
%
%    Notes:
%     - Requires field LOVROK to be set to TRUE for overwriting.
%
%    Header changes: see CHECKHEADER
%
%    Examples:
%     Read in some files, clean them up and write over:
%      writeseizmo(taper(removetrend(readseizmo('*'))))
%
%     To write out all records in big endian:
%      writeseizmo(data,'byteorder','ieee-be')
%
%     Alter the location of where the files are written:
%      writeseizmo(data,'path','some/new/directory')
%
%    See also:  readseizmo, bseizmo, seizmodef, getfileversion, readdata,
%               readdatawindow, readheader, writeheader, changename,
%               changepath

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
%        Dec. 13, 2008 - added mkdir call to make sure path exists
%        Apr.  7, 2009 - LOVROK support, better messages/checks
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        May  29, 2009 - allow options via WRITEPARAMETERS
%        Sep.  5, 2009 - improved the examples
%        Sep. 18, 2009 - added empty data shortcut
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 18, 2009 at 14:35 GMT

% todo:

% check nargin
if(mod(nargin-1,2))
    error('seizmo:writeseizmo:badNumInputs',...
        'Bad number of arguments!');
end

% shortcut for empty data
if(isempty(data))
    disp('Nothing to write!')
    return;
end

% handle options
data=writeparameters(data,varargin{:});

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

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
[leven,lovrok]=getlgc(data,'leven','lovrok');

% loop over records
for i=1:length(data)
    % skip writing if dataless
    if(~data(i).hasdata)
        warning('seizmo:writeseizmo:dataless',...
            ['Record: %d, File: %s\n' ...
            'Use WRITEHEADER to write only headers!'],i,name);
        continue; 
    end
    
    % construct fullname
    name=fullfile(data(i).path,data(i).name);
    
    % make sure directory exists
    [ok,msg,msgid]=mkdir(data(i).path);
    if(~ok)
        warning(msgid,msg);
        error('seizmo:writeseizmo:pathBad',...
            ['Record: %d, File: %s\n' ...
            'Cannot write record to path!'],i,name);
    end

    % check if existing file
    if(exist(name,'file'))
        % check lovrok
        if(strcmpi(lovrok(i),'false'))
            warning('seizmo:writeseizmo:lovrokBlock',...
                ['Record: %d, File: %s\n' ...
                'LOVROK set to FALSE!\n' ...
                'Cannot Overwrite ==> Skipping!'],i,name);
            continue;
        end
        
        % check if directory
        if(exist(name,'dir'))
            error('seizmo:writeseizmo:badFID',...
                ['Record: %d, File: %s\n' ...
                'File not openable for writing!\n'...
                '(Directory Conflict!)'],i,name);
        end
    end
    
    % open file for writing
    fid=fopen(name,'w',data(i).byteorder);
    
    % check fid
    if(fid<0)
        % unopenable file for writing (permissions/directory?)
        error('seizmo:writeseizmo:badFID',...
            ['Record: %d, File: %s\n' ...
            'File not openable for writing!\n'...
            '(Permissions problem?)'],i,name);
    end
    
    % fill header with dummy bytes (so we can seek around)
    fseek(fid,0,'bof');
    count=fwrite(fid,zeros(h(vi(i)).data.startbyte,1),'char');
    
    % verify write
    if(count<h(vi(i)).data.startbyte)
        % write failed
        fclose(fid);
        error('seizmo:writeseizmo:writeFailed',...
            'Record: %d, File: %s\nWriting failed!',i,name);
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
            ['Record: %d, File: %s\n'...
            'Output file disksize does not match expected size!\n'],...
            i,name);
    end
    
    % close file
    fclose(fid);
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end
