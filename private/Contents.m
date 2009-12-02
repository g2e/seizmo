% Seismology Toolbox - seizmo_internal
% Version 0.6.0-r105 Ararat 2-Dec-2009
%
% Low-level internal functions
%CHECKPARAMETERS       - Parses options passed to CHECKHEADER
%CUTPARAMETERS         - Parses inputs defining the data window(s)
%GET_CHECKHEADER_STATE - Returns TRUE if CHECKHEADER is on, FALSE if not
%GETFILEVERSION        - Get filetype, version and byte-order of SEIZMO datafile
%GET_SEIZMOCHECK_STATE - Returns TRUE if SEIZMOCHECK is on, FALSE if not
%ISSEIZMO              - True for SEIZMO data structures
%ISVALIDSEIZMO         - TRUE for valid filetype/version combinations
%PLOTCONFIGFIX         - Fixes the SEIZMO plot configuration struct
%SEIZMOCHECK           - Validate SEIZMO data structure
%SEIZMODEF             - Returns specified SEIZMO definition structure
%SEIZMOSIZE            - Returns header-estimated disksize of SEIZMO records in bytes
%SEIZMOVERBOSE         - Turn verbose SEIZMO output on (TRUE) or off (FALSE)
%SET_CHECKHEADER_STATE - Turn CHECKHEADER on (TRUE) or off (FALSE)
%SET_SEIZMOCHECK_STATE - Turn SEIZMOCHECK on (TRUE) or off (FALSE)
%VALIDSEIZMO           - Returns valid SEIZMO datafile filetypes or versions
%VERSIONINFO           - Returns version info for SEIZMO data records
%WRITEPARAMETERS       - Implements options passed to WRITE functions

