function [varargout]=checkoperr(varargin)
%CHECKOPERR    Controls behavior of SEIZMO CHECKHEADER function
%
%    Usage:    checkoperr()
%              checkoperr('defaults')
%              checkoperr('option1','error|warn|ignore|fix|warnfix',...,
%                         'optionN','error|warn|ignore|fix|warnfix')
%              s=checkoperr(...)
%
%    Description: CHECKOPERR provides access (through display, output, or
%     modification) to the action options of CHECKHEADER.  Modifying these
%     options alters the behavior of CHECKHEADER for the entire session or
%     until another CHECKOPERR call changes the setting.  All options may
%     be set to one of five values: 'error', 'warn', 'ignore', 'fix', and
%     'warnfix'.  The first three values ('error' 'warn' and 'ignore') are
%     always available (for better or worse) while the last two ('fix' and
%     'warnfix') may not be implemented (a warning message indicates why).
%     See individual options to find out what is checked/fixed.
%
%     CHECKOPERR() displays the current CHECKHEADER settings.
%
%     CHECKOPERR('defaults') clears any previous settings to assure that
%     CHECKHEADER will behave in the default manner.
%
%     CHECKOPERR('all','error|warn|ignore|fix|warnfix') changes all
%     action options to the specified value.  Useful when chained with more
%     specific settings to focus on specialized header problems.
%
%     CHECKOPERR('option','error|warn|ignore|fix|warnfix') alters a
%     specific option to the specified state.  See the table in the Notes
%     section for valid options.
%
%     CHECKOPERR('OPTION1','ERROR|WARN|IGNORE|FIX|WARNFIX',...,
%                'OPTIONN','ERROR|WARN|IGNORE|FIX|WARNFIX') changes
%                multiple options in a single command.
%
%     S=CHECKOPERR(...) outputs the current option settings in a struct.
%
%    Notes:
%     - Valid Options:
%       INVALID_NVHDR
%           CHECKS:     VALID VERSION FIELD WITHIN HEADER
%           FIX:        CHANGE NVHDR TO STRUCT .VERSION FIELD
%           DEFAULT:    ERROR
%       INVALID_IFTYPE
%           CHECKS:     VALID DATATYPE FIELD WITHIN HEADER
%           FIX:        CHANGE IFTYPE TO GENERAL XY DATATYPE
%           DEFAULT:    ERROR
%       MULTIPLE_IFTYPE
%           CHECKS:     DATASET HAS MULTIPLE DATATYPES
%           FIX:        NONE
%           DEFAULT:    WARN
%       INVALID_UNEVEN
%           CHECKS:     UNEVENLY SAMPLED SPECTRAL OR XYZ DATATYPE
%           FIX:        NONE
%           DEFAULT:    ERROR
%       INVALID_IZTYPE
%           CHECKS:     VALID REFERENCE TIME TYPE FIELD
%           FIX:        CHANGE IZTYPE TO UNKNOWN
%           DEFAULT:    WARN
%       MULTIPLE_IZTYPE
%           CHECKS:     DATASET HAS MULTIPLE REFERENCE TIME TYPES
%           FIX:        NONE
%           DEFAULT:    WARN
%       NONZERO_IZTYPE
%           CHECKS:     ASSOCIATED FIELD IS NONZERO (TOLERANCE OF 0.001s)
%           FIX:        CHANGE IZTYPE TO UNKNOWN
%           DEFAULT:    WARN
%       INVALID_IDEP
%           CHECKS:     VALID DEPENDENT COMPONENT TYPE
%           FIX:        CHANGE IDEP TO UNKNOWN
%           DEFAULT:    WARN
%       MULTIPLE_IDEP
%           CHECKS:     DATASET HAS MULTIPLE DEPENDENT COMPONENT TYPES
%           FIX:        NONE
%           DEFAULT:    WARN
%       UNSET_LEVEN
%           CHECKS:     EVENLY SPACED SAMPLING LOGICAL, LEVEN, IS DEFINED
%           FIX:        SET TO TRUE
%           DEFAULT:    ERROR
%       MULTIPLE_LEVEN
%           CHECKS:     DATASET HAS SOME EVENLY AND SOME UNEVENLY SAMPLED
%           FIX:        NONE
%           DEFAULT:    WARN
%       MULTIPLE_DELTA
%           CHECKS:     DATASET HAS MULTIPLE SAMPLERATES
%           FIX:        NONE
%           DEFAULT:    WARN
%       NEGATIVE_DELTA
%           CHECKS:     NEGATIVE SAMPLE SPACING
%           FIX:        NONE
%           DEFAULT:    ERROR
%       ZERO_DELTA
%           CHECKS:     SAMPLE SPACING IS ZERO
%           FIX:        NONE
%           DEFAULT:    ERROR
%       NEGATIVE_NPTS
%           CHECKS:     NEGATIVE NUMBER OF POINTS
%           FIX:        SET NPTS TO ZERO
%           DEFAULT:    WARNFIX
%       ZERO_NPTS
%           CHECKS:     ZERO NUMBER OF POINTS
%           FIX:        NONE
%           DEFAULT:    WARN
%       MULCMP_DEP
%           CHECKS:     NCMP SET >1
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       INVALID_MULCMP_DEP
%           CHECKS:     FILETYPE SUPPORTS MULTIPLE DEPENDENT COMPONENTS
%           FIX:        CHANGE TO MULTI-COMPONENT CAPABLE FILETYPE/VERSION
%           DEFAULT:    WARNFIX
%       NEGATIVE_NCMP
%           CHECKS:     NCMP SET NEGATIVE
%           FIX:        SET NCMP TO ZERO
%           DEFAULT:    WARNFIX
%       ZERO_NCMP
%           CHECKS:     NCMP SET TO ZERO
%           FIX:        NONE
%           DEFAULT:    WARN
%       BAD_SPECTRAL_B
%           CHECKS:     B FIELD FOR SPECTRAL DATA IS ZERO
%           FIX:        SET B TO ZERO
%           DEFAULT:    ERROR
%       BAD_SPECTRAL_DELTA
%           CHECKS:     DELTA FIELD FOR SPECTRAL DATA IS CONSISTENT
%           FIX:        SET DELTA BASED ON E/NPTS
%           DEFAULT:    ERROR
%       BAD_SPECTRAL_E
%           CHECKS:     E FIELD FOR SPECTRAL DATA IS NYQUIST
%           FIX:        SET E BASED ON SDELTA, OR 0.5 IF SDELTA UNDEFINED
%           DEFAULT:    ERROR
%       BAD_SPECTRAL_NPTS
%           CHECKS:     NPTS FOR SPECTRAL DATA IS A POWER OF 2
%           FIX:        NONE
%           DEFAULT:    ERROR
%       BAD_SPECTRAL_NSPTS
%           CHECKS:     NSPTS FOR SPECTRAL DATA IS DEFINED, & <NPTS
%           FIX:        SET NSPTS TO NPTS
%           DEFAULT:    ERROR
%       BAD_SPECTRAL_SB
%           CHECKS:     SB FOR SPECTRAL DATA IS DEFINED
%           FIX:        SET SB TO ZERO
%           DEFAULT:    ERROR
%       BAD_SPECTRAL_SDELTA
%           CHECKS:     SDELTA IS DEFINED
%           FIX:        SET SDELTA BASED ON E, OR 1 IF E UNDEFINED
%           DEFAULT:    ERROR
%       NONINTEGER_REFTIME
%           CHECKS:     REFTIME FIELDS (NZ*) ARE IN INTEGERS
%           FIX:        SET NZ* TO ZERO (NZJDAY IS SET TO 1)
%           DEFAULT:    WARN
%       UNSET_REFTIME
%           CHECKS:     REFTIME FIELDS (NZ*) ARE DEFINED
%           FIX:        SET NZ* TO ZERO (NZJDAY IS SET TO 1)
%           DEFAULT:    WARN
%       OUTOFRANGE_REFTIME
%           CHECKS:     REFTIME FIELDS (NZ*) ARE IN DEFINED RANGES
%           FIX:        SET NZ* TO ZERO (NZJDAY IS SET TO 1)
%           DEFAULT:    WARN
%       INCONSISTENT_E
%           CHECKS:     E IS CONSISTENT WITH B & DELTA FOR EVENLY SAMPLED
%           FIX:        CHANGE E TO MATCH B & DELTA
%           DEFAULT:    FIX
%       INACCURATE_TIMING
%           CHECKS:     TIME PRECISION DEGRADATION
%           FIX:        NONE
%           DEFAULT:    WARN
%       OUTOFRANGE_LAT
%           CHECKS:     STLA/EVLA ARE IN RANGE -90 TO 90
%           FIX:        SET STLA/EVLA TO UNDEFINED
%           DEFAULT:    WARNFIX
%       OUTOFRANGE_LON
%           CHECKS:     STLO/EVLO ARE IN RANGE -180 TO 180
%           FIX:        WRAP STLO/EVLO TO -180 TO 180
%           DEFAULT:    WARNFIX
%       OLD_DELAZ
%           CHECKS:     GCARC/AZ/BAZ/DIST (IF LCALDA TRUE)
%           FIX:        UPDATE GCARC/AZ/BAZ/DIST
%           DEFAULT:    FIX
%       KM_DEPTH
%           CHECKS:     EVDP >0m & <1000m
%           FIX:        MULTIPLY EVDP BY 1000
%           DEFAULT:    WARN
%       UNSET_ELEV
%           CHECKS:     EVEL/STEL ARE DEFINED
%           FIX:        SET EVEL/STEL TO ZERO
%           DEFAULT:    IGNORE
%       UNSET_DEPTH
%           CHECKS:     EVDP/STDP ARE DEFINED
%           FIX:        SET EVDP/STDP TO ZERO
%           DEFAULT:    IGNORE
%       INCONSISTENT_DEP_NPTS
%           CHECKS:     NPTS IN HEADER MATCHES NPTS IN DEP
%           FIX:        UPDATE NPTS IN HEADER
%           DEFAULT:    WARNFIX
%       INCONSISTENT_DEP_NCMP
%           CHECKS:     NCMP IN HEADER MATCHES NCMP IN DEP
%           FIX:        UPDATE NCMP IN HEADER
%           DEFAULT:    WARNFIX
%       BAD_SPECTRAL_DEP
%           CHECKS:     DEPENDENT DATA IS 2^N BY 2*NCMP
%           FIX:        NONE
%           DEFAULT:    ERROR
%       MISSING_IND
%           CHECKS:     IND DATA CORRESPONDING TO DEP DATA EXISTS
%           FIX:        SET LEVEN TO TRUE
%           DEFAULT:    WARNFIX
%       EVEN_IND
%           CHECKS:     EVENLY SPACED DATA W/ NON-EMPTY IND DATA
%           FIX:        SET .IND TO EMPTY
%           DEFAULT:    FIX
%       MULCMP_IND
%           CHECKS:     MULTIPLE INDEPENDENT DATA IN IND
%           FIX:        NONE
%           DEFAULT:    ERROR
%       INCONSISTENT_IND_NPTS
%           CHECKS:     NUMBER OF POINTS IN IND AND DEP MATCH
%           FIX:        TRUNCATE TO SHORTER
%           DEFAULT:    ERROR
%       CMPLX_IND
%           CHECKS:     IND IS COMPLEX
%           FIX:        SET IND REAL
%           DEFAULT:    ERROR
%       NONMONOTONIC_IND
%           CHECKS:     TIME POINTS IN IND ALWAYS INCREASING/DECREASING
%           FIX:        PARALLEL SORT IND & DEP
%           DEFAULT:    ERROR
%       REPEAT_IND
%           CHECKS:     REPEAT TIME POINTS IN IND
%           FIX:        DROP ALL BUT FIRST
%           DEFAULT:    ERROR
%       INCONSISTENT_IND_B
%           CHECKS:     B MATCHES IND DATA
%           FIX:        UPDATE B
%           DEFAULT:    FIX
%       INCONSISTENT_IND_E
%           CHECKS:     E MATCHES IND DATA
%           FIX:        UPDATE E
%           DEFAULT:    FIX
%       INCONSISTENT_IND_DELTA
%           CHECKS:     DELTA MATCHES IND DATA
%           FIX:        UPDATE DELTA
%           DEFAULT:    FIX
%       CMPLX_DEP
%           CHECKS:     DEP IS COMPLEX
%           FIX:        SET DEP REAL
%           DEFAULT:    ERROR
%       OLD_DEP_STATS
%           CHECKS:     DEPMAX/DEPMEN/DEPMIN
%           FIX:        UPDATE DEP* STATS
%           DEFAULT:    FIX
%       CMPLX_HEAD
%           CHECKS:     HEAD IS COMPLEX
%           FIX:        SET HEAD REAL
%           DEFAULT:    ERROR
%
%    Examples:
%     Make CHECKHEADER only require that the version field is valid:
%      checkoperr('all','ignore','invalid_nvhdr','error')
%
%    See also: CHECKHEADER, CHECKPARAMETERS, BINOPERR

%     Version History:
%        June  4, 2009 - initial version
%        Sep. 28, 2009 - major revision for upcoming checkheader rewrite
%        Sep. 29, 2009 - uses CHECKPARAMETERS now
%        Oct.  5, 2009 - doc update
%        Oct. 16, 2009 - added complex data checks
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 16, 2009 at 19:40 GMT

% todo:

% check nargout
if(nargout>1)
    error('seizmo:checkoperr:tooManyOutputs',...
        'Too many output arguments!');
end

% pass arguments on to checkparameters
option=checkparameters(true,varargin{:});
fields=fieldnames(option).';

% no inputs/outputs = display options
if(~nargin && ~nargout)
    for i=fields; disp(sprintf('%25s = %s',i{:},option.(i{:}))); end
end

% output
if(nargout); varargout{1}=option; end

end
