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
%     Setting SETGLOBAL to FALSE will cause OPTION1 to OPTIONN to not
%     affect subsequent calls to CHECKHEADER.  See CHECKHEADER or the
%     source code below for available options.
%
%    Notes:
%
%    Examples:
%     Make all subsequent calls to CHECKHEADER only check for the version
%     stored in the header of records:
%      checkparameters(true,'all','ignore','invalid_nvhdr','error');
%
%    See also: CHECKHEADER, CHECKOPERR

%     Version History:
%        Sep. 29, 2009 - initial version taken from CHECKOPERR
%        Oct.  5, 2009 - added option table to notes, new option EVEN_IND
%        Oct. 16, 2009 - added complex data checks
%        Oct. 21, 2009 - dropped rmfield usage (slow)
%        Dec.  1, 2009 - fixed bug where global settings overrode command
%                        line settings when SETGLOBAL was set FALSE
%        Dec.  8, 2009 - new options NAN_DEP & INF_DEP
%        Jan. 28, 2010 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2010 at 00:55 GMT

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
option.NAN_DEP='ERROR';
option.INF_DEP='ERROR';
option.OLD_DEP_STATS='FIX';
option.CMPLX_HEAD='ERROR';

% make a defaults backup
defaults=option;

% get option names
fields=fieldnames(option).';

% valid states
states={'ERROR' 'WARN' 'IGNORE' 'FIX' 'WARNFIX'};

% pull up SEIZMO global
global SEIZMO

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

% parse options
if(nargin==2)
    % check input is 'defaults'
    if(~any(strcmpi({'def' 'default' 'defaults'},varargin{:})))
        error('seizmo:checkparameters:unknownOption',...
            'Unknown option or bad option usage!');
    end
    option=defaults;
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
        if(strcmp(varargin{i},'ALL'))
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

end
