function [option]=checkparameters(setglobal,varargin)
%CHECKPARAMETERS    Parses options passed to CHECKHEADER
%
%    Usage:    options=checkparameters(setglobal,'def|default|defaults')
%              options=checkparameters(setglobal,...
%                         'option1','error|warn|ignore|fix|warnfix',...,
%                         'optionN','error|warn|ignore|fix|warnfix')
%
%    Description: OPTIONS=CHECKPARAMETERS(SETGLOBAL,'DEF|DEFAULT|DEFAULTS')
%     returns the default options for CHECKHEADER.  If SETGLOBAL is TRUE,
%     then the default settings are saved to the global
%     SEIZMO.CHECKPARAMETERS.  This overwrites any previous modifications
%     from options passed with the next usage form (see below).  This will
%     cause CHECKHEADER to behave in the default manner when no options are
%     passed to it.  If SETGLOBAL is FALSE, then the default settings are
%     returned without modifying the global options.
%
%     OPTIONS=CHECKPARAMETERS(SETGLOBAL,...
%                         'OPTION1','ERROR|WARN|IGNORE|FIX|WARNFIX',...,
%                         'OPTIONN','ERROR|WARN|IGNORE|FIX|WARNFIX')
%     returns a modified set of options for CHECKHEADER, with the options
%     indicated by OPTION1 to OPTIONN set to the state specified.  Setting
%     SETGLOBAL to TRUE will also modify the behavior of CHECKHEADER for
%     all subsequent calls until the options are changed back or the
%     previous usage form resets the behavior back to the defaults.
%     Setting SETGLOBAL to FALSE will not change the saved options.  See
%     CHECKHEADER or the source code below for available options.
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
%     Make all subsequent calls to CHECKHEADER only check the version
%     stored in the header:
%      checkparameters(true,'all','ignore','invalid_nvhdr','error');
%
%    See also: CHECKHEADER, CHECKOPERR

%     Version History:
%        Sep. 29, 2009 - initial version taken from CHECKOPERR
%        Oct.  5, 2009 - added option table to notes, new option EVEN_IND
%        Oct. 16, 2009 - added complex data checks
%        Oct. 21, 2009 - dropped rmfield usage (slow)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 21, 2009 at 07:25 GMT

% todo:

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end

% check setglobal
if(~islogical(setglobal))
    error('seizmo:checkparameters:badInput','SETGLOBAL must be logical!');
end

% defaults
%version
option.INVALID_NVHDR='ERROR';
%enums
option.INVALID_IFTYPE='ERROR';
option.INVALID_UNEVEN='ERROR';
option.MULTIPLE_IFTYPE='WARN';
option.INVALID_IZTYPE='WARN';
option.NONZERO_IZTYPE='WARN';
option.MULTIPLE_IZTYPE='WARN';
option.INVALID_IDEP='WARN';
option.MULTIPLE_IDEP='WARN';
%delta
option.UNSET_LEVEN='ERROR';
option.MULTIPLE_LEVEN='WARN'; %CHECK W/ VSDATA TOO
option.MULTIPLE_DELTA='WARN';
option.NEGATIVE_DELTA='ERROR';
option.ZERO_DELTA='ERROR'; % CHECK W/ VSDATA TOO
%npts
option.NEGATIVE_NPTS='WARNFIX';
option.ZERO_NPTS='WARN'; % CHECK W/ VSDATA TOO
%ncmp
option.MULCMP_DEP='IGNORE'; % CHECK W/ VSDATA TOO
option.INVALID_MULCMP_DEP='WARNFIX'; % CHECK W/ VSDATA TOO
option.NEGATIVE_NCMP='WARNFIX';
option.ZERO_NCMP='WARN'; % CHECK W/ VSDATA TOO
%spectral
option.BAD_SPECTRAL_B='ERROR';
option.BAD_SPECTRAL_SB='ERROR';
option.BAD_SPECTRAL_SDELTA='ERROR';
option.BAD_SPECTRAL_NPTS='ERROR';
option.BAD_SPECTRAL_NSPTS='ERROR';
option.BAD_SPECTRAL_DELTA='ERROR';
option.BAD_SPECTRAL_E='ERROR';
%timing
option.NONINTEGER_REFTIME='WARN';
option.UNSET_REFTIME='WARN';
option.OUTOFRANGE_REFTIME='WARN';
option.INCONSISTENT_E='FIX';
option.INACCURATE_TIMING='WARN';
%location
option.OUTOFRANGE_LAT='WARNFIX';
option.OUTOFRANGE_LON='WARNFIX';
option.OLD_DELAZ='FIX';
option.KM_DEPTH='WARN';
option.UNSET_ELEV='IGNORE';
option.UNSET_DEPTH='IGNORE';
%vsdata
option.INCONSISTENT_DEP_NPTS='WARNFIX';
option.INCONSISTENT_DEP_NCMP='WARNFIX';
option.BAD_SPECTRAL_DEP='ERROR';
option.MISSING_IND='WARNFIX';
option.EVEN_IND='FIX';
option.MULCMP_IND='ERROR';
option.INCONSISTENT_IND_NPTS='ERROR';
option.CMPLX_IND='ERROR';
option.NONMONOTONIC_IND='ERROR';
option.REPEAT_IND='ERROR';
option.INCONSISTENT_IND_B='FIX';
option.INCONSISTENT_IND_E='FIX';
option.INCONSISTENT_IND_DELTA='FIX';
option.CMPLX_DEP='ERROR';
option.OLD_DEP_STATS='FIX';
option.CMPLX_HEAD='ERROR';

% valid states
states={'ERROR' 'WARN' 'IGNORE' 'FIX' 'WARNFIX'};

% pull up and check over SEIZMO global
global SEIZMO; fields=fieldnames(option).';

% parse options
if(nargin==2)
    % check input is 'defaults'
    if(~any(strcmpi({'def' 'default' 'defaults'},varargin{:})))
        error('seizmo:checkparameters:unknownOption',...
            'Unknown option or bad option usage!');
    end
    if(setglobal)
        % clear SEIZMO settings
        if(isfield(SEIZMO,'CHECKPARAMETERS'))
            SEIZMO.CHECKPARAMETERS=[];
        end
    end
    return;
elseif(nargin>2)
    nin=nargin-1;
    % all inputs must be 'option','state' pairs
    if(mod(nin,2))
        error('seizmo:checkparameters:unpairedOption',...
            'Option missing a value!');
    end
    % must be valid
    varargin=upper(varargin);
    for i=varargin(1:2:end)
        if(~any(strcmp([{'ALL'} fields],i)))
            error('seizmo:checkparameters:unknownOption',...
                'Unknown option: %s',i{:});
        end
    end
    for i=varargin(2:2:end)
        if(~any(strcmp(i,states)))
            error('seizmo:checkparameters:unknownOptionState',...
                'Unknown option state: %s',i{:});
        end
    end
    % assign settings
    for i=1:2:nin
        if(strcmpi(varargin{i},'all'))
            for j=fields
                option.(j{:})=varargin{i+1};
            end
            if(setglobal)
                for j=fields
                    SEIZMO.CHECKPARAMETERS.(j{:})=varargin{i+1};
                end
            end
        else
            option.(varargin{i})=varargin{i+1};
            if(setglobal)
                SEIZMO.CHECKPARAMETERS.(varargin{i})=varargin{i+1};
            end
        end
    end
end

% pull in pre-existing settings too
if(isfield(SEIZMO,'CHECKPARAMETERS'))
    for i=fields
        if(isfield(SEIZMO.CHECKPARAMETERS,i))
            if(~any(strcmpi(SEIZMO.CHECKPARAMETERS.(i{:}),states)))
                % change improper regardless of setglobal
                warning('seizmo:checkparameters:badState',...
                    '%s in unknown state => changing to default!',i{:});
                SEIZMO.CHECKPARAMETERS.(i{:})=option.(i{:});
            else
                option.(i{:})=upper(SEIZMO.CHECKPARAMETERS.(i{:}));
            end
        end
    end
end

end
