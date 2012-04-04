function [endianness]=nativebyteorder()
%NATIVEBYTEORDER    Returns native endianness of present platform
%
%    Usage:    endianness=nativebyteorder
%
%    Description:
%     ENDIANNESS=NATIVEBYTEORDER returns the endianness of the current
%     computer using the built-in Matlab function COMPUTER.  Output is a
%     string of either 'ieee-le' or 'ieee-be' (works with FREAD/FWRITE
%     statements).
%
%    Notes:
%
%    Examples:
%     % Change SEIZMO files in current directory to current platform's
%     % byte-ordering:
%     writeseizmo(readseizmo('*'),'byteorder',nativebyteorder)
%
%    See also: CHANGEBYTEORDER, WRITESEIZMO, READSEIZMO, BSEIZMO

%     Version History:
%        Sep. 25, 2008 - initial version
%        Nov. 15, 2008 - update for new name schema (now NATIVEBYTEORDER)
%        Apr. 23, 2009 - move usage up
%        Sep.  8, 2009 - update example
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 05:45 GMT

% todo:

% get this platform's native byte-order
[platform,maxint,endianness]=computer;
clear platform maxint

% put in a SEIZMO consumable format (for fread/fwrite)
if(strcmpi(endianness,'L')); endianness='ieee-le';
else endianness='ieee-be'; end

end
