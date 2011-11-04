function [varargout]=checkoperr(varargin)
%CHECKOPERR    Controls behavior of SEIZMO CHECKHEADER function
%
%    Usage:    checkoperr()
%              checkoperr('defaults')
%              checkoperr('option1','error|warn|ignore|fix|warnfix',...,
%                         'optionN','error|warn|ignore|fix|warnfix')
%              s=checkoperr(...)
%
%    Description:
%     CHECKOPERR provides access (through display, output, or modification)
%     to the action options of CHECKHEADER.  Modifying these options alters
%     the behavior of CHECKHEADER for the entire session or until another
%     CHECKOPERR call changes the setting.  All options may be set to one
%     of five values: 'error', 'warn', 'ignore', 'fix', and 'warnfix'.  The
%     first three values ('error' 'warn' and 'ignore') are always available
%     (for better or worse) while the last two ('fix' and 'warnfix') may
%     not be implemented (a warning message will usually indicate why that
%     is so).
%
%     CHECKOPERR() displays the current CHECKHEADER settings.
%
%     CHECKOPERR('DEFAULTS') clears any previous settings to assure that
%     CHECKHEADER will behave in the default manner.
%
%     CHECKOPERR('ALL','ERROR|WARN|IGNORE|FIX|WARNFIX') changes all
%     action options to the specified value.  Useful when chained with more
%     specific settings to focus on specialized header problems.
%
%     CHECKOPERR('OPTION','ERROR|WARN|IGNORE|FIX|WARNFIX') alters a
%     specific option to the specified state.  See the table in the Notes
%     section of CHECKHEADER for valid options.
%
%     CHECKOPERR('OPTION1','ERROR|WARN|IGNORE|FIX|WARNFIX',...,
%                'OPTIONN','ERROR|WARN|IGNORE|FIX|WARNFIX') changes
%                multiple options in a single command.
%
%     S=CHECKOPERR(...) outputs the current option settings in a struct.
%
%    Notes:
%
%    Examples:
%     % Make CHECKHEADER only require that the version field is valid:
%     checkoperr('all','ignore','invalid_nvhdr','error')
%
%    See also: CHECKHEADER, CHECKPARAMETERS, BINOPERR

%     Version History:
%        June  4, 2009 - initial version
%        Sep. 28, 2009 - major revision for upcoming checkheader rewrite
%        Sep. 29, 2009 - uses CHECKPARAMETERS now
%        Oct.  5, 2009 - doc update
%        Oct. 16, 2009 - added complex data checks
%        Dec.  8, 2009 - new options NAN_DEP & INF_DEP
%        Jan. 28, 2010 - doc update
%        Feb. 11, 2011 - fprintf fix
%        Nov.  2, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  2, 2011 at 00:55 GMT

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
    for i=fields; fprintf('%25s = %s\n',i{:},option.(i{:})); end
end

% output
if(nargout); varargout{1}=option; end

end
