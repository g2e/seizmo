function [data]=changebyteorder(data,endianness)
%CHANGEBYTEORDER    Change the byteorder of SEIZMO data records
%
%    Usage:    data=changebyteorder(data,endianness)
%
%    Description: CHANGEBYTEORDER(DATA,ENDIANNESS) changes the byte-order
%     that the records in the SEIZMO struct DATA will be written as to
%     ENDIANNESS.  ENDIANNESS must be the string 'ieee-le' or 'ieee-be' or
%     it may be a char/cellstr array of those strings to define each
%     record's endianness individually.  ENDIANNESS may also be specified
%     as 'little', 'big', 'intel', 'pc', or 'sun'.
%
%    Notes:
%     - empty string ('') will preserve the byteorder of records
%
%    Examples:
%     Change records in current directory to the platform's byte-ordering:
%      writeseizmo(changebyteorder(readseizmo('*'),nativebyteorder))
%
%    See also: NATIVEBYTEORDER, WRITESEIZMO, READSEIZMO, BSEIZMO,
%              READDATA, READDATAWINDOW, READHEADER, WRITEHEADER

%     Version History:
%        Sep. 25, 2008 - initial version
%        Nov. 16, 2008 - rename from CENDIAN to CHANGEBYTEORDER
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        May  28, 2009 - fixed to update byteorder field, minor doc update
%        May  29, 2009 - allow empty endianness (no change)
%        Jan. 28, 2010 - seizmoverbose support
%        Apr. 25, 2010 - added new strings
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 25, 2010 at 14:05 GMT

% todo:

% check number of inputs
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% verbosity
verbose=seizmoverbose;

% number of records
nrecs=numel(data);

% detail message
if(verbose)
    disp('Changing Byte Order of Record(s)');
    print_time_left(0,nrecs);
end

% fast exit
if(isempty(endianness))
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end
    return;
end

% check and fix type
if(~iscellstr(endianness))
    if(~ischar(endianness))
        error('seizmo:changebyteorder:badInput',...
            'ENDIANNESS must be a char or cellstr array!');
    else
        endianness=cellstr(endianness);
    end
end

% expand scalar
if(isscalar(endianness))
    endianness=endianness(ones(numel(data),1));
elseif(numel(data)~=numel(endianness))
    error('seizmo:changebyteorder:badInput',...
        'ENDIANNESS must be scalar or match the size of DATA!');
end

% find those to preserve
preserve=strcmp(endianness,'');
[endianness{preserve}]=deal(data(preserve).byteorder);

% check endianness
endianness=lower(endianness);
little=ismember(endianness,{'ieee-le' 'intel' 'little' 'pc'});
endianness(little)={'ieee-le'};
big=ismember(endianness,{'ieee-be' 'sun' 'big'});
endianness(big)={'ieee-be'};
if(any(~(little | big)))
    error('seizmo:changebyteorder:badEndian',...
        'ENDIANNESS must be ''ieee-le'' or ''ieee-be''!');
end

% change endianness
[data.byteorder]=deal(endianness{:});

% detail message
if(verbose); print_time_left(nrecs,nrecs); end

end
