function [data]=checkheader(data,varargin)
%CHECKHEADER    Check and fix header values of SEIZMO records
%
%    Usage: data=checkheader(data)
%           data=checkheader(data,...
%                            'option1','error|warn|ignore|fix|warnfix',...,
%                            'optionN','error|warn|ignore|fix|warnfix')
%
%    Description:
%     CHECKHEADER(DATA) does a number of consistency checks on SEIZMO
%     records in DATA.  Run CHECKOPERR() to get the current option states.
%     See the table below for a description of each option.  Run
%     checkoperr('defaults') to put CHECKHEADER into the default mode.
%
%     DATA=CHECKHEADER(DATA,'OPTION1','ERROR|WARN|IGNORE|FIX|WARNFIX',...,
%     'OPTIONN','ERROR|WARN|IGNORE|FIX|WARNFIX') modifies the state of
%     options OPTION1 to OPTIONN to the specified state for this run.  To
%     alter the state of these options for the duration of the Matlab/
%     Octave session, use CHECKOPERR.
%
%    Notes:
%     - use GET_CHECKHEADER_STATE and SET_CHECKHEADER_STATE to check and
%       change if CHECKHEADER is on or off (turning CHECKHEADER off will
%       speed up most functions but will allow inconsistencies to cause
%       trouble)
%     - Valid Options:
%       INVALID_NVHDR
%           CHECKS:     VALID VERSION FIELD WITHIN HEADER
%           FIX:        CHANGE NVHDR TO STRUCT .VERSION FIELD
%           DEFAULT:    ERROR
%       INVALID_IFTYPE
%           CHECKS:     VALID DATATYPE FIELD WITHIN HEADER
%           FIX:        CHANGE IFTYPE TO GENERAL XY DATATYPE
%           DEFAULT:    ERROR
%       NONTIME_IFTYPE
%           CHECKS:     DATASET HAS NON-TIMESERIES/XY DATATYPES
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       NONSPECTRAL_IFTYPE
%           CHECKS:     DATASET HAS NON-SPECTRAL DATATYPES
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       NONXYZ_IFTYPE
%           CHECKS:     DATASET HAS NON-XYZ DATATYPES
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       TIME_IFTYPE
%           CHECKS:     DATASET HAS TIMESERIES/XY DATATYPES
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       SPECTRAL_IFTYPE
%           CHECKS:     DATASET HAS SPECTRAL DATATYPES
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       XYZ_IFTYPE
%           CHECKS:     DATASET HAS XYZ DATATYPES
%           FIX:        NONE
%           DEFAULT:    IGNORE
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
%       FALSE_LEVEN
%           CHECKS:     DATASET HAS SOME UNEVENLY SAMPLED RECORDS
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       MULTIPLE_LEVEN
%           CHECKS:     DATASET HAS SOME EVENLY AND SOME UNEVENLY SAMPLED
%           FIX:        NONE
%           DEFAULT:    WARN
%       UNSET_DELTA
%           CHECKS:     SAMPLE SPACING FIELD IS UNSET
%           FIX:        SET DELTA TO 1
%           DEFAULT:    ERROR
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
%       UNSET_B
%           CHECKS:     RELATIVE START TIME FIELD IS UNSET
%           FIX:        SET B TO ZERO
%           DEFAULT:    ERROR
%       MULTIPLE_B
%           CHECKS:     DATASET HAS MULTIPLE RELATIVE START TIMES
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       UNSET_NPTS
%           CHECKS:     NUMBER OF POINTS FIELD IS UNSET
%           FIX:        SET NPTS TO ZERO
%           DEFAULT:    ERROR
%       MULTIPLE_NPTS
%           CHECKS:     DATASET HAS MULTIPLE NUMBER OF POINTS
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       NEGATIVE_NPTS
%           CHECKS:     NEGATIVE NUMBER OF POINTS
%           FIX:        SET NPTS TO ZERO
%           DEFAULT:    ERROR
%       ZERO_NPTS
%           CHECKS:     ZERO NUMBER OF POINTS
%           FIX:        NONE
%           DEFAULT:    IGNORE
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
%       UNSET_REFTIME
%           CHECKS:     REFTIME FIELDS (NZ*) ARE DEFINED
%           FIX:        SET NZ* TO ZERO (NZJDAY IS SET TO 1)
%           DEFAULT:    WARN
%       NONINTEGER_REFTIME
%           CHECKS:     REFTIME FIELDS (NZ*) ARE IN INTEGERS
%           FIX:        SET NZ* TO ZERO (NZJDAY IS SET TO 1)
%           DEFAULT:    WARN
%       OUTOFRANGE_REFTIME
%           CHECKS:     REFTIME FIELDS (NZ*) ARE IN DEFINED RANGES
%           FIX:        SET NZ* TO ZERO (NZJDAY IS SET TO 1)
%           DEFAULT:    WARN
%       MULTIPLE_REFTIME
%           CHECKS:     REFTIME FIELDS (NZ*) ARE ALL THE SAME
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       INCONSISTENT_E
%           CHECKS:     E IS CONSISTENT WITH B & DELTA FOR EVENLY SAMPLED
%           FIX:        CHANGE E TO MATCH B & DELTA
%           DEFAULT:    FIX
%       INACCURATE_TIMING
%           CHECKS:     TIME PRECISION DEGRADATION
%           FIX:        NONE
%           DEFAULT:    WARN
%       UNSET_ST_LATLON
%           CHECKS:     STLA/STLO ARE DEFINED
%           FIX:        NONE
%           DEFAULT:    WARN
%       UNSET_EV_LATLON
%           CHECKS:     EVLA/EVLO ARE DEFINED
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       OUTOFRANGE_LAT
%           CHECKS:     STLA/EVLA ARE IN RANGE -90 TO 90
%           FIX:        WRAP STLA/EVLA TO -90 TO 90 (SHIFTS STLO/EVLO TOO)
%           DEFAULT:    WARNFIX
%       OUTOFRANGE_LON
%           CHECKS:     STLO/EVLO ARE IN RANGE -180 TO 180
%           FIX:        WRAP STLO/EVLO TO -180 TO 180
%           DEFAULT:    WARNFIX
%       OLD_DELAZ
%           CHECKS:     GCARC/AZ/BAZ/DIST (IF LCALDA TRUE)
%           FIX:        UPDATE GCARC/AZ/BAZ/DIST
%           DEFAULT:    FIX
%       UNSET_DELAZ
%           CHECKS:     DELAZ FIELDS ARE DEFINED
%           FIX:        NONE
%           DEFAULT:    IGNORE
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
%       NAN_DEP
%           CHECKS:     DEP HAS NANS
%           FIX:        REMOVE RECORD
%           DEFAULT:    ERROR
%       INF_DEP
%           CHECKS:     DEP HAS INFINITE VALUES
%           FIX:        REMOVE RECORD
%           DEFAULT:    ERROR
%       REPEAT_DEP
%           CHECKS:     DEP HAS SUCCESSIVE REPEATED VALUES
%           FIX:        NONE
%           DEFAULT:    IGNORE
%       OLD_DEP_STATS
%           CHECKS:     DEPMAX/DEPMEN/DEPMIN
%           FIX:        UPDATE DEP* STATS
%           DEFAULT:    FIX
%       CMPLX_HEAD
%           CHECKS:     HEAD IS COMPLEX
%           FIX:        SET HEAD REAL
%           DEFAULT:    ERROR
%
%    Header changes: NVHDR, IFTYPE, IZTYPE, IDEP,
%                    DELTA, NPTS, NCMP, B, E, LEVEN, SDELTA, SB, NSPTS,
%                    NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC, NZMSEC
%                    DEPMEN, DEPMIN, DEPMAX,
%                    EVLA, EVLO, EVEL, EVDP, STLA, STLO, STEL, STDP,
%                    GCARC, AZ, BAZ, DIST
%
%    Examples:
%     % This may help point out header problems:
%     checkheader(data,'all','warn')
%
%    See also: LISTHEADER, CHANGEHEADER, GETHEADER, READHEADER, READDATA,
%              SEIZMOCHECK, CHECKOPERR, CHECKPARAMETERS

%     Version History:
%        Feb. 21, 2008 - initial version
%        Feb. 23, 2008 - renamed to chkhdr, some improvements to
%                        handle unevenly spaced records
%        Feb. 25, 2008 - Couple bugfixes
%        May  12, 2008 - Use new dep* formula
%        June 15, 2008 - Updated documentation
%        June 20, 2008 - Updates header then checks
%        June 22, 2008 - Major revision - supports LCALDA field, better
%                        handling of unevenly sampled data
%        June 28, 2008 - .dep and .ind rather than .x and .t
%        Oct. 15, 2008 - broke into sections, fixed location bugs, now does
%                        location calculations very much like SAC
%        Oct. 26, 2008 - reorganize sections, recursion to allow sections
%                        to build off on another, doc update
%        Dec.  4, 2008 - total rewrite, dropped recursion for better
%                        sequencing of checks, significant variable sharing
%                        to reduce repeated calls, more checks, support for
%                        turning this function on/off through SEIZMO global
%        Apr.  7, 2009 - fixed LOVROK handling (not checked here anymore),
%                        try/catch for quicker on/off flag check
%        Apr. 23, 2009 - fix seizmocheck for octave, move usage up
%        May  10, 2009 - add support for expanded idep set
%        June  8, 2009 - couple fixes for non-column vector data
%        June 23, 2009 - fixed enum checking bug that only returned one
%                        record per bad enum value
%        Sep. 18, 2009 - added warning for ghassan about datasets with
%                        multiple samplerates
%        Sep. 25, 2009 - fixed/updated multi-cmp code, allow iztype<1e-3
%        Oct.  5, 2009 - major revision
%        Oct.  6, 2009 - drop ARRAYFUN/LOGICAL functions (R14SP1 fix),
%                        nonzero_iztype fix
%        Oct.  7, 2009 - minor changes to a few messages
%        Oct. 13, 2009 - fix ncmp check for spectral records
%        Oct. 15, 2009 - fix for inaccurate_timing (big speed jump)
%        Oct. 16, 2009 - added complex data checks
%        Nov. 13, 2009 - update for geodetic to geographic
%        Dec.  8, 2009 - new options NAN_DEP & INF_DEP
%        Jan. 29, 2010 - forgot to predefine destruction arrays, minor
%                        changing of messages, proper SEIZMO handling,
%                        fixed EVEN_IND bug, VERSIONINFO caching
%        Feb.  2, 2010 - updates VERSIONINFO cache (requires calling
%                        function to update its variables), error if empty
%                        dataset is created
%        Feb.  3, 2010 - no longer remove records w/ nans or infs - they
%                        are converted to zero (if fix/warnfix)
%        Feb. 16, 2010 - OUTOFRANGE_LAT now shifts (ev/st)lo as necessary
%        Feb. 28, 2010 - added REPEAT_DEP, all FIX with no solution are now
%                        WARNFIX equivalent
%        Mar. 15, 2010 - allow NZMSEC to be 1000
%        May   5, 2010 - allow icounts for idep field
%        May  10, 2010 - added NONTIME_IFTYPE, NONSPECTRAL_IFTYPE,
%                        NONXYZ_IFTYPE, UNSET_B, UNSET_NPTS, UNSET_DELTA,
%                        UNSET_ST_LATLON, UNSET_EV_LATLON, FALSE_LEVEN,
%                        XYZ_IFTYPE, MULTIPLE_REFTIME, MULTIPLE_NPTS,
%                        MULTIPLE_B.  Edited several defaults.
%        June 10, 2010 - fixed bug in MULTIPLE_B
%        Aug. 21, 2010 - update for getheader no longer returning undefined
%                        values specific to a filetype, KM_DEPTH fixed
%        Sep. 16, 2010 - added UNSET_DELAZ
%        Feb. 11, 2011 - minor doc update
%        Feb. 25, 2011 - fix nonzero_iztype warning
%        Nov.  2, 2011 - added TIME_IFTYPE & SPECTRAL_IFTYPE, bugfix for
%                        NONXYZ_IFTYPE
%        Dec.  1, 2011 - slight adjustment to MULTIPLE_REFTIME msgs
%        Jan. 30, 2012 - better getheader usage
%        Feb. 14, 2013 - use strcmpi for consistency
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2013 at 22:10 GMT

% todo:

% check input
if(nargin<1)
    error('seizmo:checkheader:notEnoughInputs',...
        'Not enough input arguments.');
end

% check SEIZMO global for quick exit
if(seizmodebug); disp('debug checkheader 1!'); end
global SEIZMO
try
    if(~SEIZMO.CHECKHEADER.ON); return; end
catch
    % checks data below
end
if(seizmodebug); disp('debug checkheader 2!'); end

% check data structure & grab header setup
[h,vi]=versioninfo(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);
oldversioninfocache=versioninfo_cache(true);

% attempt header check
try
    % parse options
    option=checkparameters(false,varargin{:});

    % get header fields
    [nvhdr,delta,npts,ncmp,nz,b,e,sb,sdelta,nspts,st,ev,delaz,dep,...
        iftype,iztype,idep,leven,lcalda]=getheader(data,...
        'nvhdr','delta','npts','ncmp','nz','b','e','sb',...
        'sdelta','nspts','st','ev','delaz','dep','iftype id',...
        'iztype id','idep id','leven lgc','lcalda lgc');

    % logicals
    xyz=strcmpi(iftype,'ixyz');
    spectral=(strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph'));

    % fudge factor
    try
        fudge=SEIZMO.CHECKHEADER.FUDGE;
    catch
        fudge=1e-6;
    end

    %%%%%%%%%%%%%%%%
    % BEGIN CHECKS %
    %%%%%%%%%%%%%%%%

    % check for invalid header version
    if(~strcmp(option.INVALID_NVHDR,'IGNORE'))
        nvhdr=invalid_nvhdr(option.INVALID_NVHDR,h,vi,nvhdr);
    end

    % START ENUMS
    % check for invalid data type
    if(~strcmp(option.INVALID_IFTYPE,'IGNORE'))
        iftype=invalid_iftype(option.INVALID_IFTYPE,iftype);
    end
    
    % check for nontime data types
    if(~strcmp(option.NONTIME_IFTYPE,'IGNORE'))
        nontime_iftype(option.NONTIME_IFTYPE,iftype);
    end
    
    % check for nonspectral data types
    if(~strcmp(option.NONSPECTRAL_IFTYPE,'IGNORE'))
        nonspectral_iftype(option.NONSPECTRAL_IFTYPE,iftype);
    end
    
    % check for nonxyz data types
    if(~strcmp(option.NONXYZ_IFTYPE,'IGNORE'))
        nonxyz_iftype(option.NONXYZ_IFTYPE,iftype);
    end
    
    % check for time data types
    if(~strcmp(option.TIME_IFTYPE,'IGNORE'))
        time_iftype(option.TIME_IFTYPE,iftype);
    end
    
    % check for spectral data types
    if(~strcmp(option.SPECTRAL_IFTYPE,'IGNORE'))
        spectral_iftype(option.SPECTRAL_IFTYPE,iftype);
    end
    
    % check for xyz data types
    if(~strcmp(option.XYZ_IFTYPE,'IGNORE'))
        xyz_iftype(option.XYZ_IFTYPE,iftype);
    end

    % check for multiple data types
    if(~strcmp(option.MULTIPLE_IFTYPE,'IGNORE'))
        multiple_iftype(option.MULTIPLE_IFTYPE,iftype);
    end

    % check for undefined leven
    if(~strcmp(option.UNSET_LEVEN,'IGNORE'))
        leven=unset_leven(option.UNSET_LEVEN,leven);
    end

    % check for invalid uneven
    if(~strcmp(option.INVALID_UNEVEN,'IGNORE'))
        invalid_uneven(option.INVALID_UNEVEN,spectral | xyz,leven);
    end

    % check for invalid reference time type
    if(~strcmp(option.INVALID_IZTYPE,'IGNORE'))
        iztype=invalid_iztype(option.INVALID_IZTYPE,iztype);
    end

    % check for invalid reftime fields
    if(~strcmp(option.UNSET_REFTIME,'IGNORE'))
        nz=unset_reftime(option.UNSET_REFTIME,nz);
    end
    
    % check for invalid reftime fields
    if(~strcmp(option.NONINTEGER_REFTIME,'IGNORE'))
        nz=noninteger_reftime(option.NONINTEGER_REFTIME,nz);
    end

    % check for invalid reftime fields
    if(~strcmp(option.OUTOFRANGE_REFTIME,'IGNORE'))
        nz=outofrange_reftime(option.OUTOFRANGE_REFTIME,nz);
    end
    
    % check for multiple reftime
    if(~strcmp(option.MULTIPLE_REFTIME,'IGNORE'))
        multiple_reftime(option.MULTIPLE_REFTIME,nz);
    end

    % check for nonzero iztype
    if(~strcmp(option.NONZERO_IZTYPE,'IGNORE'))
        SEIZMO.VERSIONINFO.USECACHE=false;
        [iztype]=nonzero_iztype(...
            option.NONZERO_IZTYPE,spectral,data,nz,iztype);
        SEIZMO.VERSIONINFO.USECACHE=true;
        SEIZMO.VERSIONINFO.H=h;
        SEIZMO.VERSIONINFO.IDX=vi;
    end

    % check for multiple reference time types
    if(~strcmp(option.MULTIPLE_IZTYPE,'IGNORE'))
        multiple_iztype(option.MULTIPLE_IZTYPE,iztype);
    end

    % check for invalid dependent component type
    if(~strcmp(option.INVALID_IDEP,'IGNORE'))
        idep=invalid_idep(option.INVALID_IDEP,idep);
    end

    % check for multiple dependent component types
    if(~strcmp(option.MULTIPLE_IDEP,'IGNORE'))
        multiple_idep(option.MULTIPLE_IDEP,idep);
    end
    % END ENUMS

    % START TIMING
    % check for unset b
    if(~strcmp(option.UNSET_B,'IGNORE'))
        b=unset_b(option.UNSET_B,b);
    end
    
    % check for inconsistent e
    if(~strcmp(option.INCONSISTENT_E,'IGNORE'))
        e=inconsistent_e(...
            option.INCONSISTENT_E,leven,spectral,b,delta,npts,e);
    end

    % check for inaccurate timing
    if(~strcmp(option.INACCURATE_TIMING,'IGNORE'))
        inaccurate_timing(option.INACCURATE_TIMING,b,delta,e,fudge);
    end
    % END TIMING

    % START SPECTRAL
    % check for bad spectral b
    if(~strcmp(option.BAD_SPECTRAL_B,'IGNORE'))
        b=bad_spectral_b(option.BAD_SPECTRAL_B,spectral,b);
    end

    % check for bad spectral sdelta
    if(~strcmp(option.BAD_SPECTRAL_SDELTA,'IGNORE'))
        sdelta=bad_spectral_sdelta(...
            option.BAD_SPECTRAL_SDELTA,spectral,e,sdelta);
    end

    % check for bad spectral e
    if(~strcmp(option.BAD_SPECTRAL_E,'IGNORE'))
        e=bad_spectral_e(option.BAD_SPECTRAL_E,spectral,sdelta,e,fudge);
    end

    % check for bad spectral npts
    if(~strcmp(option.BAD_SPECTRAL_NPTS,'IGNORE'))
        bad_spectral_npts(option.BAD_SPECTRAL_NPTS,spectral,npts);
    end

    % check for bad spectral nspts
    if(~strcmp(option.BAD_SPECTRAL_NSPTS,'IGNORE'))
        nspts=bad_spectral_nspts(...
            option.BAD_SPECTRAL_NSPTS,spectral,npts,nspts);
    end

    % check for bad spectral delta
    if(~strcmp(option.BAD_SPECTRAL_DELTA,'IGNORE'))
        delta=bad_spectral_delta(...
            option.BAD_SPECTRAL_DELTA,spectral,e,npts,delta,fudge);
    end

    % check for bad spectral sb
    if(~strcmp(option.BAD_SPECTRAL_SB,'IGNORE'))
        sb=bad_spectral_sb(option.BAD_SPECTRAL_SB,spectral,sb);
    end
    % END SPECTRAL

    % START LOCATION
    % check for unset station lat/lon
    if(~strcmp(option.UNSET_ST_LATLON,'IGNORE'))
        unset_st_latlon(option.UNSET_ST_LATLON,st);
    end
    
    % check for unset event lat/lon
    if(~strcmp(option.UNSET_EV_LATLON,'IGNORE'))
        unset_ev_latlon(option.UNSET_EV_LATLON,ev);
    end
    
    % check for bad lat
    if(~strcmp(option.OUTOFRANGE_LAT,'IGNORE'))
        [ev,st]=outofrange_lat(option.OUTOFRANGE_LAT,ev,st);
    end

    % check for bad lon
    if(~strcmp(option.OUTOFRANGE_LON,'IGNORE'))
        [ev,st]=outofrange_lon(option.OUTOFRANGE_LON,ev,st);
    end

    % check for undef depth
    if(~strcmp(option.UNSET_DEPTH,'IGNORE'))
        [ev,st]=unset_depth(option.UNSET_DEPTH,ev,st);
    end

    % check for undef elev
    if(~strcmp(option.UNSET_ELEV,'IGNORE'))
        [ev,st]=unset_elev(option.UNSET_ELEV,ev,st);
    end

    % check for km depths
    if(~strcmp(option.KM_DEPTH,'IGNORE'))
        ev=km_depth(option.KM_DEPTH,ev);
    end

    % check for old delaz
    if(~strcmp(option.OLD_DELAZ,'IGNORE'))
        delaz=old_delaz(option.OLD_DELAZ,lcalda,st,ev,delaz);
    end
    
    % check for unset delaz
    if(~strcmp(option.UNSET_DELAZ,'IGNORE'))
        unset_delaz(option.UNSET_DELAZ,delaz);
    end
    % END LOCATION

    % START VSDATA
    % check for mismatch between ncmp and size of dep
    if(~strcmp(option.INCONSISTENT_DEP_NCMP,'IGNORE'))
        ncmp=inconsistent_dep_ncmp(...
            option.INCONSISTENT_DEP_NCMP,spectral,data,ncmp);
    end

    % check for multiple component records
    if(~strcmp(option.MULCMP_DEP,'IGNORE'))
        mulcmp_dep(option.MULCMP_DEP,ncmp);
    end

    % check for multiple component records using a non-multi-cmp filetype
    if(~strcmp(option.INVALID_MULCMP_DEP,'IGNORE'))
        % NOTE THAT THIS UPDATES H & VI BOTH HERE AND IN THE SEIZMO CACHE
        % - ANY FUNCTION THAT CALLS CHECKHEADER WILL THUS NEED TO UPDATE
        %   ITS LOCAL H & VI TO MATCH THE CACHE!
        % - IGNORE THE UNUSED VARIABLE MLINT MESSAGE B/C WE MAY NEED H & VI
        %   FOR A TO-BE-ADDED TEST BELOW THIS ONE
        [data,h,vi,nvhdr]=invalid_mulcmp_dep(...
            option.INVALID_MULCMP_DEP,data,h,vi,nvhdr,ncmp);
    end

    % check for negative ncmp
    if(~strcmp(option.NEGATIVE_NCMP,'IGNORE'))
        ncmp=negative_ncmp(option.NEGATIVE_NCMP,ncmp);
    end

    % check for no cmp records
    if(~strcmp(option.ZERO_NCMP,'IGNORE'))
        zero_ncmp(option.ZERO_NCMP,ncmp);
    end

    % check for mismatch between ncmp and size of dep for spectral
    if(~strcmp(option.BAD_SPECTRAL_DEP,'IGNORE'))
        bad_spectral_dep(option.BAD_SPECTRAL_DEP,spectral,data);
    end

    % check for missing ind
    if(~strcmp(option.MISSING_IND,'IGNORE'))
        leven=missing_ind(option.MISSING_IND,leven,data);
    end
    
    % check for false leven
    if(~strcmp(option.FALSE_LEVEN,'IGNORE'))
        false_leven(option.FALSE_LEVEN,leven);
    end
    
    % check for dataset with mix of evenly & unevenly sampled records
    if(~strcmp(option.MULTIPLE_LEVEN,'IGNORE'))
        multiple_leven(option.MULTIPLE_LEVEN,leven);
    end

    % clear even with ind
    if(~strcmp(option.EVEN_IND,'IGNORE'))
        data=even_ind(option.EVEN_IND,leven,data);
    end

    % check for multiple independent datum
    if(~strcmp(option.MULCMP_IND,'IGNORE'))
        mulcmp_ind(option.MULCMP_IND,leven,data);
    end

    % check for mismatch between dep size and ind size
    if(~strcmp(option.INCONSISTENT_IND_NPTS,'IGNORE'))
        data=inconsistent_ind_npts(...
            option.INCONSISTENT_IND_NPTS,leven,data);
    end

    % check for complex ind
    if(~strcmp(option.CMPLX_IND,'IGNORE'))
        data=cmplx_ind(option.CMPLX_IND,data);
    end

    % check for nonmonotonic ind
    if(~strcmp(option.NONMONOTONIC_IND,'IGNORE'))
        data=nonmonotonic_ind(option.NONMONOTONIC_IND,leven,data);
    end

    % check for repeat ind
    if(~strcmp(option.REPEAT_IND,'IGNORE'))
        data=repeat_ind(option.REPEAT_IND,leven,data);
    end

    % check for mismatch between npts and size of dep
    if(~strcmp(option.INCONSISTENT_DEP_NPTS,'IGNORE'))
        npts=inconsistent_dep_npts(...
            option.INCONSISTENT_DEP_NPTS,data,npts);
    end
    
    % check for unset npts
    if(~strcmp(option.UNSET_NPTS,'IGNORE'))
        npts=unset_npts(option.UNSET_NPTS,npts);
    end
    % check for multiple npts
    if(~strcmp(option.MULTIPLE_NPTS,'IGNORE'))
        multiple_npts(option.MULTIPLE_NPTS,npts);
    end
    
    % check for negative npts
    if(~strcmp(option.NEGATIVE_NPTS,'IGNORE'))
        npts=negative_npts(option.NEGATIVE_NPTS,npts);
    end

    % check for zero npts
    if(~strcmp(option.ZERO_NPTS,'IGNORE'))
        zero_npts(option.ZERO_NPTS,npts);
    end

    % check for mismatch between b and ind
    if(~strcmp(option.INCONSISTENT_IND_B,'IGNORE'))
        b=inconsistent_ind_b(...
            option.INCONSISTENT_IND_B,leven,data,b);
    end

    % check for mismatch between e and ind
    if(~strcmp(option.INCONSISTENT_IND_E,'IGNORE'))
        e=inconsistent_ind_e(...
            option.INCONSISTENT_IND_E,leven,data,e);
    end
    
    % multiple b
    if(~strcmp(option.MULTIPLE_B,'IGNORE'))
        multiple_b(option.MULTIPLE_B,b);
    end

    % check for mismatch between delta and ind
    if(~strcmp(option.INCONSISTENT_IND_DELTA,'IGNORE'))
        delta=inconsistent_ind_delta(...
            option.INCONSISTENT_IND_DELTA,leven,data,delta);
    end
    
    % check for unset delta
    if(~strcmp(option.UNSET_DELTA,'IGNORE'))
        delta=unset_delta(option.UNSET_DELTA,delta);
    end

    % check for negative delta
    if(~strcmp(option.NEGATIVE_DELTA,'IGNORE'))
        negative_delta(option.NEGATIVE_DELTA,delta);
    end

    % check for zero delta
    if(~strcmp(option.ZERO_DELTA,'IGNORE'))
        zero_delta(option.ZERO_DELTA,delta);
    end

    % check for multiple sample rates
    if(~strcmp(option.MULTIPLE_DELTA,'IGNORE'))
        multiple_delta(option.MULTIPLE_DELTA,delta);
    end

    % check for complex dep
    if(~strcmp(option.CMPLX_DEP,'IGNORE'))
        data=cmplx_dep(option.CMPLX_DEP,data);
    end

    % check for nan dep
    if(~strcmp(option.NAN_DEP,'IGNORE'))
        data=nan_dep(option.NAN_DEP,data);
    end

    % check for inf dep
    if(~strcmp(option.INF_DEP,'IGNORE'))
        data=inf_dep(option.INF_DEP,data);
    end
    
    % check for repeat dep
    if(~strcmp(option.REPEAT_DEP,'IGNORE'))
        data=repeat_dep(option.REPEAT_DEP,data);
    end

    % check for old dep stats
    if(~strcmp(option.OLD_DEP_STATS,'IGNORE'))
        dep=old_dep_stats(option.OLD_DEP_STATS,data,dep);
    end

    % check for complex head
    if(~strcmp(option.CMPLX_HEAD,'IGNORE'))
        data=cmplx_head(option.CMPLX_HEAD,data);
    end
    % END VSDATA

    %{
    % check for 
    if(~strcmp(option.,'IGNORE'))
        (option.,);
    end
    %}

    % update header
    data=changeheader(data,'nvhdr',nvhdr,'iftype',iftype,...
        'iztype',iztype,'idep',idep,'leven',leven,'delta',delta,...
        'npts',npts,'ncmp',ncmp,'nz',nz,'b',b,'e',e,'sb',sb,...
        'sdelta',sdelta,'nspts',nspts,'st',st,'ev',ev,'delaz',delaz,...
        'dep',dep);

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
	versioninfo_cache(oldversioninfocache);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
	versioninfo_cache(oldversioninfocache);
    
    % rethrow error
    error(lasterror);
end

end

function [nvhdr]=invalid_nvhdr(opt,h,vi,nvhdr)
bad=find([h(vi).version].'~=nvhdr);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:NVHDRbad';
    report.message=['Version info corrupted! NVHDR for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nis inconsistent with the version stored in the struct!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            nvhdr=[h(vi).version].';
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting invalid NVHDR to match struct!');
            nvhdr=[h(vi).version].';
    end
end

end

function [iftype]=invalid_iftype(opt,iftype)
validftype={'itime' 'irlim' 'iamph' 'ixy' 'ixyz'};
bad=find(~ismember(iftype,validftype));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:IFTYPEbad';
    report.message=['IFTYPE unknown/undefined for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nMust be one of the following:\n' sprintf('%s ',validftype{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            iftype(bad)={'ixy'};
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting invalid IFTYPE to ''ixy''!');
            iftype(bad)={'ixy'};
    end
end

end

function nontime_iftype(opt,iftype)
validftype={'itime' 'ixy'};
bad=find(~ismember(iftype,validftype));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:nontimeIFTYPE';
    report.message=['IFTYPE for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nMust be one of the following:\n' sprintf('%s ',validftype{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A TIMESERIES DATA TYPE.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A TIMESERIES DATA TYPE.']);
    end
end

end

function nonspectral_iftype(opt,iftype)
validftype={'irlim' 'iamph'};
bad=find(~ismember(iftype,validftype));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:nonspectralIFTYPE';
    report.message=['IFTYPE for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nMust be one of the following:\n' sprintf('%s ',validftype{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A SPECTRAL DATA TYPE.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A SPECTRAL DATA TYPE.']);
    end
end

end

function nonxyz_iftype(opt,iftype)
validftype={'ixyz'};
bad=find(~ismember(iftype,validftype));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:nonxyzIFTYPE';
    report.message=['IFTYPE for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nMust be one of the following:\n' sprintf('%s ',validftype{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO THE XYZ DATA TYPE.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO THE XYZ DATA TYPE.']);
    end
end

end

function time_iftype(opt,iftype)
invalidftype={'itime' 'ixy'};
bad=find(ismember(iftype,invalidftype));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:timeIFTYPE';
    report.message=['IFTYPE for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nCannot be one of the following:\n' ...
        sprintf('%s ',invalidftype{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A NON-TIME DATA TYPE.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A NON-TIME DATA TYPE.']);
    end
end

end

function spectral_iftype(opt,iftype)
invalidftype={'iamph' 'irlim'};
bad=find(ismember(iftype,invalidftype));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:spectralIFTYPE';
    report.message=['IFTYPE for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nCannot be one of the following:\n' ...
        sprintf('%s ',invalidftype{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A NON-SPECTRAL DATA TYPE.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A NON-SPECTRAL DATA TYPE.']);
    end
end

end

function xyz_iftype(opt,iftype)
invalidftype={'ixyz'};
bad=find(ismember(iftype,invalidftype));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:xyzIFTYPE';
    report.message=['IFTYPE for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nCannot be one of the following:\n' ...
        sprintf('%s ',invalidftype{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A NON-XYZ DATA TYPE.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR INCORRECT DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A NON-XYZ DATA TYPE.']);
    end
end

end

function multiple_iftype(opt,iftype)
bad=unique(iftype);
if(~isscalar(bad))
    report.identifier='seizmo:checkheader:multiIFTYPE';
    report.message=['Dataset has multiple datatypes:\n' ...
        sprintf('%s ',bad{:}) ...
        '\nThis will slow operations down and may cause problems!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A SINGLE FILETYPE.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE DATA TYPES!\n' ...
                'CONSIDER CONVERTING ALL RECORDS IN\n' ...
                'YOUR DATASET TO A SINGLE FILETYPE.']);
    end
end

end

function [leven]=unset_leven(opt,leven)
validleven={'true' 'false'};
bad=find(~ismember(leven,validleven));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:unsetLEVEN';
    report.message=['LEVEN must be TRUE or FALSE for record(s):\n'...
        sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            leven(bad)={'true'};
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting invalid LEVEN to TRUE!');
            leven(bad)={'true'};
    end
end

end

function false_leven(opt,leven)
bad=find(~ismember(leven,'true'));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:falseLEVEN';
    report.message=['LEVEN must be TRUE for record(s):\n'...
        sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            disp(['NO AUTOFIX FOR UNEVENLY SAMPLED DATA.\nCONSIDER ' ...
                'USING interpolate TO CONVERT TO EVENLY SAMPLED']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR UNEVENLY SAMPLED DATA.\nCONSIDER ' ...
                'USING interpolate TO CONVERT TO EVENLY SAMPLED']);
    end
end

end

function []=invalid_uneven(opt,reqeven,leven)
bad=find(reqeven & strcmpi(leven,'false'));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:badUneven';
    report.message=['LEVEN must be TRUE for Spectral/XYZ record(s):\n'...
        sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR UNEVEN SPECTRAL/XYZ! FIX IT MANUALLY.');
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR UNEVEN SPECTRAL/XYZ! FIX IT MANUALLY.');
    end
end

end

function [iztype]=invalid_iztype(opt,iztype)
validztype={'iunkn' 'ib' 'iday' 'io' 'ia' 'it0' 'it1' 'it2' 'it3' 'it4' ...
    'it5' 'it6' 'it7' 'it8' 'it9'};
bad=find(~ismember(iztype,validztype));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:IZTYPEbad';
    report.message=['IZTYPE unknown/undefined for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nMust be one of the following:\n' sprintf('%s ',validztype{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            iztype(bad)={'iunkn'};
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting invalid IZTYPE to ''iunkn''!');
            iztype(bad)={'iunkn'};
    end
end

end

function [nz]=noninteger_reftime(opt,nz)
bad=nz~=fix(nz) & ~(isnan(nz) | isinf(nz));
badrec=find(sum(bad,2)~=0);
if(~isempty(badrec))
    report.identifier='seizmo:checkheader:nonintREF';
    report.message=['Reference time fields for record(s):\n' ...
        sprintf('%d ',badrec) ...
        'stored as noninteger values!  This will cause problems.'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            nz=fix(nz(bad));
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Truncating decimal portion of NZ* fields!');
            nz=fix(nz(bad));
    end
end

end

function [nz]=unset_reftime(opt,nz)
bad=isnan(nz) | isinf(nz);
badrec=find(sum(bad,2)~=0);
bad=find(bad);
column=1+fix((bad-1)/size(bad,1));
if(~isempty(badrec))
    report.identifier='seizmo:checkheader:undefREF';
    report.message=['Reference time fields for record(s):\n' ...
        sprintf('%d ',badrec) ...
        '\nhave undefined values!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            nz(bad(column~=2))=0;
            nz(bad(column==2))=1;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting NZ* fields to 0 (and NZJDAY to 1)!');
            nz(bad(column~=2))=0;
            nz(bad(column==2))=1;
    end
end

end

function [nz]=outofrange_reftime(opt,nz)
bad1=false(size(nz,1),1);
bad2=nz(:,2)<1 | nz(:,2)>366 | (nz(:,2)==366 & ~isleapyear(nz(:,1)));
bad3=nz(:,3)<0 | nz(:,3)>23;
bad4=nz(:,4)<0 | nz(:,4)>59;
bad5=nz(:,5)<0 | nz(:,5)>60; % allow for leapsecond & precision error
bad6=nz(:,6)<0 | nz(:,6)>1000; % allow nzmsec to be 1000
bad=find(bad1 | bad2 | bad3 | bad4 | bad5 | bad6);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:outofrangeRef';
    report.message=['Record(s):\n' sprintf('%d ',bad)...
        '\nReference time fields NZ* not in defined range!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR OUT OF RANGE REFERENCE TIME!\n' ...
                'CONSIDER USING fixtimes.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR OUT OF RANGE REFERENCE TIME!\n' ...
                'CONSIDER USING fixtimes.']);
    end
end

end

function multiple_reftime(opt,nz)
bad=unique(nz,'rows');
if(size(bad,1)~=1)
    report.identifier='seizmo:checkheader:multiREFTIME';
    report.message=['Dataset has multiple reference times! Consider\n' ...
        'using SYNCHRONIZE to fix this! Current Reference Times:\n' ...
        sprintf('%04d:%03d@%02d:%02d:%02d.%03d ',bad')];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR MULTIPLE REFERENCE TIMES!');
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR MULTIPLE REFERENCE TIMES!');
    end
end
end

function [iztype]=nonzero_iztype(opt,spectral,data,nz,iztype)
day=strcmpi('iday',iztype);
badday=find(day & sum(nz(:,3:6)~=0,2)~=0);
unkn=strcmpi('iunkn',iztype) | ~strncmpi('i',iztype,1);
other=find(~day & ~unkn);
b=strcmpi('ib',iztype);
if(any(spectral & b)); iztype(spectral & b)={'isb'}; end
if(isscalar(unique(iztype(other))))
    iz=unique(iztype(other));
    izt=getheader(data(other),iz{1}(2:end));
    bad=other(abs(izt)>0.001);
else
    bad=true(numel(other),1);
    for i=1:numel(other)
        bad(i)=abs(getheader(data(other(i)),...
            iztype{other(i)}(2:end)))>0.001;
    end
end
if(any(spectral & b)); iztype(spectral & b)={'ib'}; end
bad=[find(bad(:)); badday(:)];
if(~isempty(bad))
    report.identifier='seizmo:checkheader:nonzeroIZTYPE';
    report.message=['Header inconsistency found! Record(s):\n' ...
        sprintf('%d ',bad)...
        '\nNZHOUR, NZMIN, NZSEC, NZMSEC are not zero for IZTYPE=IDAY\n' ...
        'or IZTYPE corresponds to a header field and that field is\n' ...
        'not zero (must be within 0.001s)!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            iztype(bad)={'iunkn'};
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting IZTYPE to ''iunkn'' for nonzero!');
            iztype(bad)={'iunkn'};
    end
end

end

function multiple_iztype(opt,iztype)
bad=unique(iztype);
if(~isscalar(bad))
    report.identifier='seizmo:checkheader:multiIZTYPE';
    report.message=['Dataset has multiple reference time types:\n' ...
        sprintf('%s ',bad{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE REFERENCE TIME TYPES! ' ...
                'CONSIDER SETTING ALL TO ''iunkn'' ' ...
                'OR USING synchronize.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE REFERENCE TIME TYPES! ' ...
                'CONSIDER SETTING ALL TO ''iunkn'' ' ...
                'OR USING synchronize.']);
    end
end

end

function [idep]=invalid_idep(opt,idep)
validdep={'iunkn' 'idisp' 'ivel' 'iacc' 'ivolts' 'iabsmnt' ...
        'iabsity' 'iabseler' 'iabserk' 'iabsnap' 'iabsackl' 'iabspop' ...
        'ijerk' 'isnap' 'icrackle' 'ipop' 'icounts'};
bad=find(~ismember(idep,validdep));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:IDEPbad';
    % only show first 5 (the rest are kinda crazy)
    report.message=['IDEP unknown/undefined for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nMust be one of the following:\n' sprintf('%s ',validdep{1:5})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            idep(bad)={'iunkn'};
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting invalid IDEP to ''iunkn''!');
            idep(bad)={'iunkn'};
    end
end

end

function multiple_idep(opt,idep)
bad=unique(idep);
if(~isscalar(bad))
    report.identifier='seizmo:checkheader:multiIDEP';
    report.message=['Dataset has multiple dependent component types:\n' ...
        sprintf('%s ',bad{:})];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE DEPENDENT COMPONENT TYPES! ' ...
                'CONSIDER SETTING ALL TO ''iunkn'' ' ...
                'OR USING integrate & differentiate.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE DEPENDENT COMPONENT TYPES! ' ...
                'CONSIDER SETTING ALL TO ''iunkn'' ' ...
                'OR USING integrate & differentiate.']);
    end
end

end

function multiple_leven(opt,leven)
bad=unique(leven);
if(~isscalar(bad))
    report.identifier='seizmo:checkheader:multiLEVEN';
    report.message=['Dataset has some evenly & unevenly sampled '...
        'record(s)!\n' ...
        'This will slow some operations down and may cause problems!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE SAMPLING TYPES!\n' ...
                'CONSIDER USING interpolate TO EVENLY SAMPLE UNEVEN.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE SAMPLING TYPES!\n' ...
                'CONSIDER USING interpolate TO EVENLY SAMPLE UNEVEN.']);
    end
end

end

function multiple_delta(opt,delta)
bad=unique(delta);
if(~isscalar(bad))
    report.identifier='seizmo:checkheader:multiDELTA';
    report.message=['Dataset has multiple sample intervals:\n' ...
        sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE SAMPLE INTERVALS! ' ...
                'CONSIDER USING syncrates.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE SAMPLE INTERVALS! ' ...
                'CONSIDER USING syncrates.']);
    end
end

end

function multiple_b(opt,b)
bad=unique(b);
if(~isscalar(bad))
    report.identifier='seizmo:checkheader:multiB';
    report.message=['Dataset has multiple begin times:\n' ...
        sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE BEGIN TIMES! ' ...
                'CONSIDER USING cut & interpolate.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE SAMPLE INTERVALS! ' ...
                'CONSIDER USING cut & interpolate.']);
    end
end

end

function multiple_npts(opt,npts)
bad=unique(npts);
if(~isscalar(bad))
    report.identifier='seizmo:checkheader:multiNPTS';
    report.message=['Dataset has multiple number of points:\n' ...
        sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE NUMBER OF POINTS! ' ...
                'CONSIDER USING cut.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR MULTIPLE NUMBER OF POINTS! ' ...
                'CONSIDER USING cut.']);
    end
end

end

function negative_delta(opt,delta)
bad=find(delta<0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:negativeDELTA';
    report.message=['DELTA for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nis NEGATIVE (time decreases in the record).\n' ...
        'This will cause problems!'];
    report.message='';
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR NEGATIVE SAMPLE INTERVALS!\n' ...
                'CONSIDER USING reverse & changing B, E, DELTA']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR NEGATIVE SAMPLE INTERVALS!\n' ...
                'CONSIDER USING reverse & changing B, E, DELTA']);
    end
end

end

function zero_delta(opt,delta)
bad=find(delta==0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:zeroDELTA';
    report.message=['DELTA for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nis ZERO (all points have the same time in the record).\n' ...
        'This will likely cause problems!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR ZERO DELTA! '...
                'TRY TO FIGURE OUT THE CAUSE.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR ZERO DELTA! '...
                'TRY TO FIGURE OUT THE CAUSE.']);
    end
end

end

function [npts]=negative_npts(opt,npts)
bad=find(npts<0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:negativeNPTS';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nNPTS is negative!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            npts(npts<0)=0;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting negative NPTS to 0!');
            npts(npts<0)=0;
    end
end

end

function zero_npts(opt,npts)
bad=find(npts==0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:zeroNPTS';
    report.message=['Record(s):\n' sprintf('%d ',bad) '\nNPTS is zero!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR ZERO NPTS! TRY TO FIGURE OUT THE CAUSE.');
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR ZERO NPTS! TRY TO FIGURE OUT THE CAUSE.');
    end
end

end

function mulcmp_dep(opt,ncmp)
bad=find(ncmp>1);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:NCMPgt1';
    report.message=['Record(s):\n' sprintf('%d ',bad) '\nNCMP>1!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR NCMP>1! CONSIDER USING splitrecords.');
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR NCMP>1! CONSIDER USING splitrecords.');
    end
end

end

function [data,h,vi,nvhdr]=invalid_mulcmp_dep(opt,data,h,vi,nvhdr,ncmp)
global SEIZMO
bad=find(ncmp>1 & ~getsubfield(h(vi),'mulcmp','valid').');
if(~isempty(bad))
    report.identifier='seizmo:checkheader:versNotMulCmp';
    report.message=['Record(s):\n' sprintf('%d ',bad)...
            '\nDo not support multiple components!\n'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            mulcmp=getsubfield(h(vi(bad)),'mulcmp');
            nvhdr(bad)=[mulcmp.altver];
            [data(bad).filetype]=deal(mulcmp.alttype);
            [data(bad).version]=deal(mulcmp.altver);
            SEIZMO.VERSIONINFO.USECACHE=false;
            [h,vi]=versioninfo(data);
            SEIZMO.VERSIONINFO.USECACHE=true;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Changing to a multi-component version!');
            mulcmp=getsubfield(h(vi(bad)),'mulcmp');
            nvhdr(bad)=[mulcmp.altver];
            [data(bad).filetype]=deal(mulcmp.alttype);
            [data(bad).version]=deal(mulcmp.altver);
            SEIZMO.VERSIONINFO.USECACHE=false;
            [h,vi]=versioninfo(data);
            SEIZMO.VERSIONINFO.USECACHE=true;
    end
end

end

function [ncmp]=negative_ncmp(opt,ncmp)
bad=find(ncmp<0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:negativeNCMP';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nNCMP is negative!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            ncmp(ncmp<0)=0;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting negative NCMP to 0!');
            ncmp(ncmp<0)=0;
    end
end

end

function zero_ncmp(opt,ncmp)
bad=find(ncmp==0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:zeroNCMP';
    report.message=['Record(s):\n' sprintf('%d ',bad) '\nNCMP is zero!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR ZERO NCMP! TRY TO FIGURE OUT THE CAUSE.');
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR ZERO NCMP! TRY TO FIGURE OUT THE CAUSE.');
    end
end

end

function b=unset_b(opt,b)
bad=find(isnan(b) | isinf(b));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:undefB';
    report.message=['B field for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nhave undefined values!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            b(:)=0;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting B field to 0!');
            b(:)=0;
    end
end

end

function npts=unset_npts(opt,npts)
bad=find(isnan(npts) | isinf(npts));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:undefNPTS';
    report.message=['NPTS field for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nhave undefined values!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            npts(:)=0;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting NPTS field to 0!');
            npts(:)=0;
    end
end

end

function delta=unset_delta(opt,delta)
bad=find(isnan(delta) | isinf(delta));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:undefDELTA';
    report.message=['DELTA field for record(s):\n' ...
        sprintf('%d ',bad) ...
        '\nhave undefined values!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            delta(:)=1;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting DELTA field to 1!');
            delta(:)=1;
    end
end

end

function [e]=inconsistent_e(opt,leven,spectral,b,delta,npts,e)
bad=find(~strcmpi(leven,'false') & ~spectral & e~=b+(npts-1).*delta);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:inconsistentE';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nE is inconsistent with that expected from B,NPTS,DELTA!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            e(bad)=b(bad)+(npts(bad)-1).*delta(bad);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting E to match expected!');
            e(bad)=b(bad)+(npts(bad)-1).*delta(bad);
    end
end

end

function inaccurate_timing(opt,b,delta,e,fudge)
bad=find((b+delta-b-delta)./delta>fudge | (e+delta-e-delta)./delta>fudge);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:timingDegraded';
    report.message=['Record(s):\n' sprintf('%d ',bad)...
        '\nDELTA is too insignificant to be well resolved!\n'...
        'This can introduce significant numerical error!\n'...
        'It is recommended to adjust the relative timing so\n'...
        'that the sample interval can be resolved.  This may\n'...
        'require splitting the record into smaller segments\n'...
        'or possibly just changing the reference time to be\n'...
        'closer to the time range of the record.'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR INACCURATE TIMING!');
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR INACCURATE TIMING!');
    end
end

end

function [b]=bad_spectral_b(opt,spectral,b)
bad=find(spectral & b~=0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:badSpectralB';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nB for spectral records must be 0!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            b(bad)=0;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting spectral B to 0!');
            b(bad)=0;
    end
end

end

function [sdelta]=bad_spectral_sdelta(opt,spectral,e,sdelta)
bad=spectral & (isnan(sdelta) | isinf(sdelta));
badd=find(bad & (isnan(e) | isinf(e)));
bad=find(bad);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:badSpectralSDELTA';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nSDELTA for spectral records must be defined!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            sdelta(bad)=1./(2*e);
            if(any(badd)); sdelta(badd)=1; end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['==> Setting SDELTA based on E ' ...
                '(set to 1 if E is undef)!']);
            sdelta(bad)=1./(2*e);
            if(any(badd)); sdelta(badd)=1; end
    end
end

end

function [e]=bad_spectral_e(opt,spectral,sdelta,e,fudge)
bad=find(spectral & abs((e-1./(2*sdelta))./e)>fudge);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:badSpectralE';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nE for spectral records gives the Nyquist frequency\n'...
        'and must be consistent with SDELTA (E=1/(2*SDELTA))!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            e(bad)=1./(2*sdelta(bad));
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['==> Setting E based on SDELTA ' ...
                '(set to 0.5 if SDELTA is undef)!']);
            e(bad)=1./(2*sdelta(bad));
    end
end

end

function bad_spectral_npts(opt,spectral,npts)
bad=find(spectral & npts~=2.^ceil(log2(npts)));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:badSpectralNPTS';
    report.message=['Record(s):\n' sprintf('%d ',bad)...
        '\nNPTS field for spectral records must be a power of 2!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO FIX FOR SPECTRAL NPTS NOT A POWER OF 2!\n' ...
                'TRY TO FIGURE OUT THE CAUSE.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO FIX FOR SPECTRAL NPTS NOT A POWER OF 2!\n' ...
                'TRY TO FIGURE OUT THE CAUSE.']);
    end
end

end

function [nspts]=bad_spectral_nspts(opt,spectral,npts,nspts)
bad=find(spectral & (nspts<0 | nspts>npts));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nNSPTS field for spectral records gives the number of\n' ...
        'points in the time domain and must be >=0 and <=NPTS!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            nspts(bad)=npts(bad);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting invalid NSPTS to equal NPTS!');
            nspts(bad)=npts(bad);
    end
end

end

function [delta]=bad_spectral_delta(opt,spectral,e,npts,delta,fudge)
bad=find(spectral & abs((delta-2*e./npts)./delta)>fudge);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:badSpectralDELTA';
    report.message=['Record(s):\n' sprintf('%d ',bad)...
        '\nDELTA for spectral records gives the frequency\n'...
        'resolution and must be consistent with E and NPTS!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            delta(bad)=2*e(bad)./npts(bad);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting DELTA based on E & NPTS!');
            delta(bad)=2*e(bad)./npts(bad);
    end
end

end

function [sb]=bad_spectral_sb(opt,spectral,sb)
bad=find(spectral & (isnan(sb) | isinf(sb)));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:badSpectralSB';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nSB for spectral records must be defined!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            sb(bad)=0;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting SB to 0!');
            sb(bad)=0;
    end
end

end

function unset_st_latlon(opt,st)
bad=isnan(st(:,1:2)) | isinf(st(:,1:2));
bad=find(sum(bad,2)~=0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:undefSTLALO';
    report.message=['STLA/STLO are undefined for record(s):\n' ...
        sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR UNSET STLA/STLO.');
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR UNSET STLA/STLO.');
    end
end

end

function unset_ev_latlon(opt,ev)
bad=isnan(ev(:,1:2)) | isinf(ev(:,1:2));
bad=find(sum(bad,2)~=0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:undefEVLALO';
    report.message=['EVLA/EVLO are undefined for record(s):\n' ...
        sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR UNSET EVLA/EVLO.');
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('NO AUTOFIX FOR UNSET EVLA/EVLO.');
    end
end
end

function [ev,st]=outofrange_lat(opt,ev,st)
badev=abs(ev(:,1))>90 & ...
    ~(isnan(ev(:,1)) | isinf(ev(:,1)));
badst=abs(st(:,1))>90 & ...
    ~(isnan(st(:,1)) | isinf(st(:,1)));
bad=find(badev | badst);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:outOfRangeLat';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nEVLA/STLA must be in range -90 to 90!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            % unwrap lats and shift lons appropriately
            [ev(badev,1),ev(badev,2)]=fixlatlon(ev(badev,1),ev(badev,2));
            [st(badst,1),st(badst,2)]=fixlatlon(st(badst,1),st(badst,2));
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Unwrapping out of range latitudes!');
            % unwrap lats and shift lons appropriately
            [ev(badev,1),ev(badev,2)]=fixlatlon(ev(badev,1),ev(badev,2));
            [st(badst,1),st(badst,2)]=fixlatlon(st(badst,1),st(badst,2));
    end
end

end

function [ev,st]=outofrange_lon(opt,ev,st)
badev=abs(ev(:,2))>180 & ...
    ~(isnan(ev(:,2)) | isinf(ev(:,2)));
badst=abs(st(:,2))>180 & ...
    ~(isnan(st(:,2)) | isinf(st(:,2)));
bad=find(badev | badst);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:outOfRangeLon';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nEVLO/STLO must be in range -180 to 180!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            ev(badev,2)=lonmod(ev(badev,2),360);
            st(badst,2)=lonmod(st(badst,2),360);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Unwrapping out of range longitudes!');
            ev(badev,2)=lonmod(ev(badev,2),360);
            st(badst,2)=lonmod(st(badst,2),360);
    end
end

end

function [ev,st]=unset_depth(opt,ev,st)
badev=isnan(ev(:,4)) | isinf(ev(:,4));
badst=isnan(st(:,4)) | isinf(st(:,4));
bad=find(badst | badev);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:unsetDepth';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nEVDP/STDP must be defined!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            ev(badev,4)=0;
            st(badst,4)=0;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting undefined depths to 0!');
            ev(badev,4)=0;
            st(badst,4)=0;
    end
end

end

function [ev,st]=unset_elev(opt,ev,st)
badev=isnan(ev(:,3)) | isinf(ev(:,3));
badst=isnan(st(:,3)) | isinf(st(:,3));
bad=find(badst | badev);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:unsetElevation';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nEVEL/STEL must be defined!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            ev(badev,3)=0;
            st(badst,3)=0;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting undefined elevations to 0!');
            ev(badev,3)=0;
            st(badst,3)=0;
    end
end

end

function [ev]=km_depth(opt,ev)
bad=find(ev(:,4)<1000 & ev(:,4)>0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:possibleKmDepth';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nEVDP must be in meters but appears to be in km!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            ev(bad,4)=ev(bad,4)*1000;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Converting kilometer depth to meters!');
            ev(bad,4)=ev(bad,4)*1000;
    end
end

end

function [delaz]=old_delaz(opt,lcalda,st,ev,delaz)

% get defined st/ev
def=0==sum(isnan(st(:,1:2)) | isinf(st(:,1:2)),2) ...
    & 0==sum(isnan(ev(:,1:2)) | isinf(ev(:,1:2)),2);

% get lcalda (don't mess with those with lcalda == 'false')
chk=~strcmpi(lcalda,'false');
newdelaz=nan(size(delaz));
ok=chk & def;
if(any(ok))
    % get geocentric latitude
    geocevla=geographic2geocentriclat(ev(ok,1));
    geocstla=geographic2geocentriclat(st(ok,1));
    
    % get gcarc, az, baz based on sphere (great-circle-arc)
    [newdelaz(ok,1),newdelaz(ok,2),newdelaz(ok,3)]=...
        sphericalinv(geocevla,ev(ok,2),geocstla,st(ok,2));
    
    % get km dist based on ellipsoid (geodesic)
    newdelaz(ok,4)=vincentyinv(ev(ok,1),ev(ok,2),st(ok,1),st(ok,2));
end

bad=find(sum(delaz~=newdelaz,2)~=0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:oldDelAz';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nGCARC, AZ, BAZ and/or DIST inconsistent w/ STLA/LO & EVLA/LO!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            % replace all to be checked
            delaz(chk,:)=newdelaz(chk,:);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Updating DELAZ values to match locations!');
            delaz(chk,:)=newdelaz(chk,:);
    end
end

end

function unset_delaz(opt,delaz)
baddelaz=sum(isnan(delaz) | isinf(delaz),2)>0;
bad=find(baddelaz);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:unsetDelAz';
    report.message=['DELAZ must be defined!\n' ...
        'Record(s):\n' sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR UNSET DELAZ!\n' ...
                'SETTING THE LCALDA FIELD TO TRUE MAY FIX THIS.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTOFIX FOR UNSET DELAZ!\n' ...
                'SETTING THE LCALDA FIELD TO TRUE MAY FIX THIS.']);
    end
end
end

function [npts]=inconsistent_dep_npts(opt,data,npts)
ok=find([data.hasdata].');
nrows=nan(size(npts));
for i=ok.'; nrows(i)=size(data(i).dep,1); end
bad=ok(nrows(ok)~=npts(ok));
%[nrows]=arrayfun(@(x)size(x.dep,1),data(ok));
%bad=nrows~=npts(ok);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:NPTSinconsistent';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nNPTS does not match data!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            npts(bad)=nrows(bad);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting NPTS to match data!');
            npts(bad)=nrows(bad);
    end
end

end

function [ncmp]=inconsistent_dep_ncmp(opt,spectral,data,ncmp)
ok=find([data.hasdata].');
nrows=nan(size(ncmp)); ncols=nrows;
for i=ok.'; [nrows(i),ncols(i)]=size(data(i).dep); end
%[nrows,ncols]=arrayfun(@(x)size(x.dep),data(ok));
ncols(spectral)=ncols(spectral)./2;
bad=ok(ncols(ok)~=ncmp(ok));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:NCMPinconsistent';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nNCMP does not match data!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            ncmp(bad)=ncols(bad);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting NCMP to match data!');
            ncmp(bad)=ncols(bad);
    end
end

end

function bad_spectral_dep(opt,spectral,data)
ok=find([data.hasdata].' & spectral);
nrows=nan(size(ok)); ncols=nrows;
for i=1:numel(ok); [nrows(i),ncols(i)]=size(data(ok(i)).dep); end
%[nrows,ncols]=arrayfun(@(x)size(x.dep),data(ok));
%bad=mod(ncols,2) | nrows~=2.^ceil(log2(nrows));
bad=ok(mod(ncols,2) | nrows~=2.^ceil(log2(nrows)));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:badDepSize';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
        '\nSpectral records must have even number of columns' ...
        '\nand NPTS equal to a power of two!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO FIX FOR IMPROPERLY SIZED SPECTRAL DATA!\n' ...
                'THIS WILL HAVE TO BE DELT WITH MANUALLY.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO FIX FOR IMPROPERLY SIZED SPECTRAL DATA!\n' ...
                'THIS WILL HAVE TO BE DELT WITH MANUALLY.']);
    end
end

end

function [leven]=missing_ind(opt,leven,data)
ok=find(strcmpi(leven,'false') & [data.hasdata].');
nrows=nan(size(ok)); ncols=nrows; nirows=nrows; nicols=nrows;
for i=1:numel(ok)
    [nrows(i),ncols(i)]=size(data(ok(i)).dep);
    [nirows(i),nicols(i)]=size(data(ok(i)).ind);
end
%[nrows,ncols]=arrayfun(@(x)size(x.dep),data(ok));
%[nirows,nicols]=arrayfun(@(x)size(x.ind),data(ok));
bad=ok(~(nirows.*nicols) & (nrows.*ncols));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:noIND';
    report.message=['Record(s):\n' sprintf('%d ',bad)...
            '\nLEVEN set FALSE requires an independent component\n'...
            'dataset but there is none.'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            leven(bad)=true;
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting LEVEN to TRUE!');
            leven(bad)=true;
    end
end

end

function [data]=even_ind(opt,leven,data)
ok=find(~strcmpi(leven,'false') & [data.hasdata].');
bad=false(size(ok));
for i=1:numel(ok); bad(i)=~isempty(data(ok(i)).ind); end
bad=ok(bad);
%bad=ok(~arrayfun(@(x)isempty(x.ind),data(ok)));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:unnecessaryIND';
    report.message=['Record(s):\n' sprintf('%d ',bad)...
            '\nLEVEN set TRUE but independent component data present.'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            for i=bad.'
                data(i).ind=[];
            end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Clearing .ind!');
            for i=bad.'
                data(i).ind=[];
            end
    end
end

end

function mulcmp_ind(opt,leven,data)
ok=find(strcmpi(leven,'false') & [data.hasdata].');
nirows=nan(size(ok)); nicols=nirows;
for i=1:numel(ok); [nirows(i),nicols(i)]=size(data(ok(i)).ind); end
%[nirows,nicols]=arrayfun(@(x)size(x.ind),data(ok));
bad=ok(nicols>1);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:badNumIndCmp';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
                '\nToo many independent components (only 1 allowed)!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO FIX FOR 2+ INDEPENDENT COMPONENTS!\n' ...
                'THIS WILL HAVE TO BE DELT WITH MANUALLY.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO FIX FOR 2+ INDEPENDENT COMPONENTS!\n' ...
                'THIS WILL HAVE TO BE DELT WITH MANUALLY.']);
    end
end

end

function [data]=inconsistent_ind_npts(opt,leven,data)
ok=find(strcmpi(leven,'false') & [data.hasdata].');
nrows=nan(size(ok)); nirows=nrows;
for i=1:numel(ok)
    nrows(i)=size(data(ok(i)).dep,1);
    nirows(i)=size(data(ok(i)).ind,1);
end
%nrows=arrayfun(@(x)size(x.dep,1),data(ok));
%nirows=arrayfun(@(x)size(x.ind,1),data(ok));
bad=ok(nrows~=nirows);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:inconsistentNPTS';
    report.message=['Record(s):\n' sprintf('%d ',bad)...
                '\nNPTS inconsistent between IND and DEP data!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            for i=bad.';
                if(nirows(i)>nrows(i))
                    data(i).ind=data(i).ind(1:nrows(i));
                else
                    data(i).dep=data(i).dep(1:nirows(i),:);
                end
            end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Truncating longer dataset to length of shorter!');
            for i=bad.';
                if(nirows(i)>nrows(i))
                    data(i).ind=data(i).ind(1:nrows(i));
                else
                    data(i).dep=data(i).dep(1:nirows(i),:);
                end
            end
    end
end

end

function [data]=nonmonotonic_ind(opt,leven,data)
ok=find(strcmpi(leven,'false') & [data.hasdata].');
bad=false(size(ok));
for i=1:numel(ok)
    if(~isempty(data(ok(i)).ind)); bad(i)=min(diff(data(ok(i)).ind))<0; end
end
bad=ok(bad);
%bad=ok(arrayfun(@(x)double(min(diff(x.ind))),data(ok))<0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:nonmonotonicIND';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
                '\nIndependent component must be monotonic!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            for i=bad'
                [data(i).ind,idx]=sort(data(i).ind);
                data(i).dep=data(i).dep(idx,:);
            end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Sorting IND and DEP!');
            for i=bad'
                [data(i).ind,idx]=sort(data(i).ind);
                data(i).dep=data(i).dep(idx,:);
            end
    end
end

end

function [data]=repeat_ind(opt,leven,data)
ok=find(strcmpi(leven,'false') & [data.hasdata].');
bad=false(size(ok));
for i=1:numel(ok); bad(i)=any(diff(sort(data(ok(i)).ind))==0); end
bad=ok(bad);
%bad=ok(arrayfun(@(x)any(diff(sort(x.ind))==0),data(ok)));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:repeatIND';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
                '\nIndependent component must not have repeated values!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            for i=bad.'
                [data(i).ind,idx]=unique(data(i).ind);
                data(i).dep=data(i).dep(idx,:);
            end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Sorting and dropping repeated values!');
            for i=bad.'
                [data(i).ind,idx]=unique(data(i).ind);
                data(i).dep=data(i).dep(idx,:);
            end
    end
end

end

function [b]=inconsistent_ind_b(opt,leven,data,b)
newb=nan(size(b));
ok=find(strcmpi(leven,'false') & [data.hasdata].');
for i=ok.'
    if(~isempty(data(i).ind)); newb(i)=double(data(i).ind(1)); end
end
%empty=arrayfun(@(x)isempty(x.ind),data(ok));
%newb(ok(~empty))=arrayfun(@(x)double(x.ind(1)),data(ok(~empty)));
bad=ok(b(ok)~=newb(ok));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:';
    report.message='';
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            b(bad)=newb(bad);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting B to match data!');
            b(bad)=newb(bad);
    end
end

end

function [e]=inconsistent_ind_e(opt,leven,data,e)
newe=nan(size(e));
ok=find(strcmpi(leven,'false') & [data.hasdata].');
for i=ok.'
    if(~isempty(data(i).ind)); newe(i)=double(data(i).ind(end)); end
end
%empty=arrayfun(@(x)isempty(x.ind),data(ok));
%newe(ok(~empty))=arrayfun(@(x)double(x.ind(end)),data(ok(~empty)));
bad=ok(e(ok)~=newe(ok));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:oldE';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
                '\nE does not match data!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            e(bad)=newe(bad);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting E to match data!');
            e(bad)=newe(bad);
    end
end

end

function [delta]=inconsistent_ind_delta(opt,leven,data,delta)
ok=find(strcmpi(leven,'false') & [data.hasdata].');
bad=false(size(ok)); newd=nan(size(delta));
for i=1:numel(ok)
    if(numel(data(ok(i)).ind)>1)
        bad(i)=true;
        newd(ok(i))=double(data(ok(i)).ind(end)-data(ok(i)).ind(1)) ...
            /(size(data(ok(i)).ind,1)-1);
    end
end
ok=ok(bad);
%ok=ok(arrayfun(@(x)numel(x.ind)>1,data(ok)));
%newd=arrayfun(@(x)double(x.ind(end)-x.ind(1))/(size(x.ind,1)-1),data(ok));
bad=ok(delta(ok)~=newd(ok));
if(~isempty(bad))
    report.identifier='seizmo:checkheader:oldDelta';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
                '\nDELTA does not match data!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            delta(bad)=newd(bad);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Setting DELTA to match data!');
            delta(bad)=newd(bad);
    end
end

end

function [dep]=old_dep_stats(opt,data,dep)
newdep=nan(numel(data),3);
ok=find([data.hasdata].');
for i=ok.'
    if(~isempty(data(i).dep))
        newdep(i,:)=[double(min(data(i).dep(:))) ...
            nanmean(double(data(i).dep(:))) ...
            double(max(data(i).dep(:)))];
    end
end
%empty=arrayfun(@(x)isempty(x.dep),data(ok));
%newdep(ok(~empty),:)=[ ...
%    arrayfun(@(x)double(min(x.dep(:))),data(ok(~empty))) ...
%    arrayfun(@(x)nanmean(double(x.dep(:))),data(ok(~empty))) ...
%    arrayfun(@(x)double(max(x.dep(:))),data(ok(~empty)))];
bad=ok(sum(dep(ok,:)~=newdep(ok,:),2)~=0);
if(~isempty(bad))
    report.identifier='seizmo:checkheader:oldDepStats';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
                '\nDEP stats do not match data!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            dep(bad,:)=newdep(bad,:);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Updating DEP stats to match data!');
            dep(bad,:)=newdep(bad,:);
    end
end

end

function [data]=cmplx_head(opt,data)
nrecs=numel(data);
bad=false(nrecs,1);
for i=1:nrecs
    bad(i)=~isreal(data(i).head);
end
bad=find(bad).';
if(~isempty(bad))
    report.identifier='seizmo:checkheader:cmplxHeader';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
                '\nHeader is complex (not real valued)!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            for i=bad
                data(i).head=real(data(i).head);
            end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Stripping imaginary values from header!');
            for i=bad
                data(i).head=real(data(i).head);
            end
    end
end

end

function [data]=cmplx_ind(opt,data)
nrecs=numel(data);
bad=false(nrecs,1);
for i=1:nrecs
    bad(i)=~isreal(data(i).ind);
end
bad=find(bad).';
if(~isempty(bad))
    report.identifier='seizmo:checkheader:cmplxIND';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
                '\nIND is complex (not real valued)!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            for i=bad
                data(i).ind=real(data(i).ind);
            end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Stripping imaginary values from IND!');
            for i=bad
                data(i).ind=real(data(i).ind);
            end
    end
end

end

function [data]=cmplx_dep(opt,data)
nrecs=numel(data);
bad=false(nrecs,1);
for i=1:nrecs
    bad(i)=~isreal(data(i).dep);
end
bad=find(bad).';
if(~isempty(bad))
    report.identifier='seizmo:checkheader:cmplxDEP';
    report.message=['Record(s):\n' sprintf('%d ',bad) ...
                '\nDEP is complex (not real valued)!'];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            for i=bad
                data(i).dep=real(data(i).dep);
            end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Stripping imaginary values from DEP!');
            for i=bad
                data(i).dep=real(data(i).dep);
            end
    end
end

end

function [data]=nan_dep(opt,data)
nrecs=numel(data);
bad=false(nrecs,1);
for i=1:nrecs
    bad(i)=any(isnan(data(i).dep(:)));
end
bad=find(bad).';
if(~isempty(bad))
    report.identifier='seizmo:checkheader:nanDEP';
    report.message=['Dependent data has NaNs!' ...
        '\nRecord(s):\n' sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            for i=1:nrecs
                nans=isnan(data(i).dep(:));
                data(i).dep(nans)=0;
            end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Changing NaN values to 0 (zero)!');
            for i=1:nrecs
                nans=isnan(data(i).dep(:));
                data(i).dep(nans)=0;
            end
    end
end

end

function [data]=inf_dep(opt,data)
nrecs=numel(data);
bad=false(nrecs,1);
for i=1:nrecs
    bad(i)=any(isinf(data(i).dep(:)));
end
bad=find(bad).';
if(~isempty(bad))
    report.identifier='seizmo:checkheader:infDEP';
    report.message=['Dependent data has infinite values!' ...
        '\nRecord(s):\n' sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            for i=1:nrecs
                infs=isinf(data(i).dep(:));
                data(i).dep(infs)=0;
            end
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp('==> Changing Inf values to 0 (zero)!');
            for i=1:nrecs
                infs=isinf(data(i).dep(:));
                data(i).dep(infs)=0;
            end
    end
end

end

function [data]=repeat_dep(opt,data)
nrecs=numel(data);
bad=false(nrecs,1);
for i=1:nrecs
    % double any for multi-cmp records
    bad(i)=any(any(~diff(data(i).dep,1,1)));
end
bad=find(bad).';
if(~isempty(bad))
    report.identifier='seizmo:checkheader:repeatDEP';
    report.message=['Dependent data has successive repeated values!\n' ...
        'Record(s):\n' sprintf('%d ',bad)];
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            warning(report.identifier,report.message);
            disp(['NO AUTO-FIX FOR SUCCESSIVE REPEATED VALUES!\n' ...
                'THIS MAY INDICATE DATA CLIPPING, FILLING, DROPOUTS,\n' ...
                'OR A MERGING ISSUE.  IT ALSO COULD JUST BE RANDOM\n' ...
                'CHANCE (ACTUALLY THIS IS USUALLY THE CASE).\n' ...
                'THIS WILL HAVE TO BE CHECKED & DELT WITH MANUALLY.']);
        case 'WARNFIX'
            warning(report.identifier,report.message);
            disp(['NO AUTO-FIX FOR SUCCESSIVE REPEATED VALUES!\n' ...
                'THIS MAY INDICATE DATA CLIPPING, FILLING, DROPOUTS,\n' ...
                'OR A MERGING ISSUE.  IT ALSO COULD JUST BE RANDOM\n' ...
                'CHANCE (ACTUALLY THIS IS USUALLY THE CASE).\n' ...
                'THIS WILL HAVE TO BE CHECKED & DELT WITH MANUALLY.']);
    end
end

end

%{
function []=(opt,)
bad=find();
if(~isempty(bad))
    report.identifier='seizmo:checkheader:';
    report.message='';
    switch opt
        case 'ERROR'
            error(report.identifier,report.message);
        case 'WARN'
            warning(report.identifier,report.message);
        case 'FIX'
            
        case 'WARNFIX'
            warning(report.identifier,report.message);
            
    end
end

end
%}
