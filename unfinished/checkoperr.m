function []=checkoperr(varargin)
%CHECKOPERR    Controls behavior of SEIZMO CHECKHEADER function
%
%    Usage:    checkoperr()
%              checkoperr('defaults')
%              checkoperr('option1','error'|'warn'|'ignore'|'fix',...,
%                         'optionN','error'|'warn'|'ignore'|'fix')
%
%    Description: Allows changing how CHECKHEADER behaves when header
%     inconsistencies are encountered.  This modifies the behavior of
%     CHECKHEADER until the session is finished or another call to
%     CHECKOPERR undoes it.  All options can be set to one of four values:
%     'error', 'warn', 'ignore', and 'fix'.  The first three values are
%     always available while the last ('fix') may not be.  See individual
%     options to find out if and how fixes are implemented.
%
%     CHECKOPERR() displays the current CHECKHEADER error settings.
%
%     CHECKOPERR('defaults') clears any previous settings to assure that
%     CHECKHEADER will behave in the default manner.
%
%     CHECKOPERR('all','error'|'warn'|'ignore'|'fix')
%     CHECKOPERR('version','error'|'warn'|'ignore'|'fix')
%     CHECKOPERR('enums','error'|'warn'|'ignore'|'fix')
%     CHECKOPERR('leven','error'|'warn'|'ignore'|'fix')
%     CHECKOPERR('delta','error'|'warn'|'ignore'|'fix')
%     CHECKOPERR('npts','error'|'warn'|'ignore'|'fix')
%     CHECKOPERR('timing','error'|'warn'|'ignore'|'fix')
%     CHECKOPERR('vsdata','error'|'warn'|'ignore'|'fix')
%     CHECKOPERR('location','error'|'warn'|'ignore'|'fix')
%
%    Notes:
%     - multiple options may be strung together in a single command
%
%    Examples:
%
%    See also: checkheader, binoperr

%     Version History:
%        June  4, 2009 - initial version
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  4, 2009 at 02:00 GMT

% todo:

% default options
option.ALL='WARN';
option.VERSION='ERROR';
option.ENUMS='WARN';
option.LEVEN='ERROR';
option.DELTA='WARN';
option.NPTS='WARN';
option.TIMING='FIX';
option.VSDATA='FIX';
option.LOCATION='FIX';

% available states
states={'ERROR' 'WARN' 'IGNORE' 'FIX'};

% pull up and check over SEIZMO global
global SEIZMO; fields=fieldnames(option).';
if(isfield(SEIZMO,'CHECKOPERR'))
    for i=fields
        if(isfield(SEIZMO.CHECKOPERR,i))
            if(~any(strcmpi(SEIZMO.CHECKOPERR.(i{:}),states)))
                warning('seizmo:checkoperr:badState',...
                    '%s in unknown state => changing to default!',i{:});
                SEIZMO.CHECKOPERR.(i{:})=option.(i{:});
            else
                option.(i{:})=upper(SEIZMO.CHECKOPERR.(i{:}));
            end
        end
    end
end

% no inputs = display options
if(nargin==0)
    for i=fields; disp(sprintf('%12s = %s',i{:},option.(i{:}))); end
% one input = 'defaults'
elseif(nargin==1)
    % check input is 'defaults'
    if(~strcmpi('defaults',varargin{:}))
        error('seizmo:checkoperr:unknownOption',...
            'Unknown option or bad option usage!');
    end
    % clear SEIZMO settings
    if(isfield(SEIZMO,'CHECKOPERR'))
        SEIZMO=rmfield(SEIZMO,'CHECKOPERR');
    end
% set options
else
    % all inputs must be 'option','state' pairs
    if(mod(nargin,2))
        error('seizmo:checkoperr:unpairedOption','Option missing a value!');
    end
    % must be valid
    varargin=upper(varargin);
    for i=varargin(1:2:end)
        if(~isfield(option,i))
            error('seizmo:checkoperr:unknownOption',...
                'Unknown option: %s',i{:});
        end
    end
    for i=varargin(2:2:end)
        if(~any(strcmpi(i,states)))
            error('seizmo:checkoperr:unknownOptionState',...
                'Unknown option state: %s',i{:});
        end
    end
    % assign settings
    for i=1:2:nargin
        SEIZMO.CHECKOPERR.(varargin{i})=varargin{i+1};
    end
end

end