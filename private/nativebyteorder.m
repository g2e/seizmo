function [endianness]=nativebyteorder()
%NATIVEBYTEORDER    Returns native endianness of present platform
%
%    Description: NATIVEBYTEORDER returns the endianness of the current
%     computer using the builtin Matlab function COMPUTER.  Output is a
%     string of either 'ieee-le' or 'ieee-be' (works with FREAD/FWRITE
%     statements).
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    endianness=nativebyteorder
%
%    Examples:
%     Change SEIZMO files in current directory to current platform's
%     byte-ordering:
%      writeseizmo(changebyteorder(readseizmo('*'),nativebyteorder))
%
%    See also: changebyteorder, writeseizmo, readseizmo, bseizmo

%     Version History:
%        Sep. 25, 2008 - initial version
%        Nov. 15, 2008 - update for new name schema (now NATIVEBYTEORDER)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 15, 2008 at 21:35 GMT

% todo:

% get this platform's native byte-order
[platform,maxint,endianness]=computer;
clear platform maxint

% put in a SEIZMO consumable format (for fread/fwrite)
if(strcmpi(endianness,'L')); endianness='ieee-le';
else endianness='ieee-be'; end

end