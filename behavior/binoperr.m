function []=binoperr(varargin)
%BINOPERR    Controls behavior of SEIZMO binary functions
%
%    Usage:    binoperr()
%              binoperr('defaults')
%              binoperr('option1',state,...,
%                       'optionN',state)
%
%    Description:
%     Allows for changing the behavior of how binary functions (addrecords,
%     subtractrecords, multiplyrecords, dividerecords) handle records with
%     some unequal aspects.  In particular, this function will modify the
%     behavior of all of these functions until the session is finished or a
%     subsequent call to BINOPERR undoes it.
%
%     BINOPERR() displays the current binary operator error settings.
%
%     BINOPERR('defaults') clears any previous settings and causes all
%     SEIZMO binary functions to use their default settings.
%     
%     BINOPERR('npts','error'|'warn'|'truncate'|'pad'|'ignore') sets how
%     binary operations handle records with different numbers of points.
%     If the option is set to 'warn' or 'ignore', the number of points in
%     the records is not altered - which will likely cause an error during
%     the operation.  If the option is set to 'truncate', the number of
%     points in the records being operated on will be equal to that with
%     the least.  Option 'pad' will make the records being operated on have
%     number of points equal to that with the most (note that padding is
%     done with zeros).  By default 'npts' is set to 'error'.
%     
%     BINOPERR('delta','error'|'warn'|'ignore') sets how binary operations
%     handle records with different sample rates.  If the option is set to
%     'warn' or 'ignore', the records are operated on point for point
%     (ignoring timing).  The resultant records' sample rates are
%     determined by the parent of their header fields (set by the option 
%     'newhdr' in the binary function call).  By default 'delta' is set to
%     'error'.
%     
%     BINOPERR('begin','error'|'warn'|'ignore') sets how binary operations
%     handle records with different begin times.  If the option is set to
%     'warn' or 'ignore', the resultant records' begin times are 
%     determined by the parent of their header fields (set by option 
%     'newhdr' in the binary function call).  By default 'begin' is set to 
%     'warn'.
%     
%     BINOPERR('ref','error'|'warn'|'ignore') sets how binary operations 
%     handle records with different reference times.  If the option is set
%     to 'warn' or 'ignore', the resultant records' reference times are 
%     determined by the parent of their header fields (set by option 
%     'newhdr' in the binary function call).  By default 'ref' is set to 
%     'warn'.
%     
%     BINOPERR('ncmp','error'|'warn'|'truncate'|'pad'|'ignore') sets how
%     binary operations handle records with different number of components.
%     If the option is set to 'warn' or 'ignore', the number of components
%     in the records is not altered - which will likely lead to an error.
%     If the option is set to 'truncate', the number of components in the 
%     records being operated on will be equal to that with the least.
%     Option 'pad' will make the number of components for records in the
%     operation equal to that of the record with the most (note that
%     padding is done with zeros).  By default 'ncmp' is set to 'error'.
%     
%     BINOPERR('leven','error'|'warn'|'ignore') sets how binary operations
%     handle unevenly sampled records.  If the option is set to 'warn' or 
%     'ignore', the functions behave as normal working point for point 
%     (basically ignoring timing).  The resultant records' leven fields are
%     determined by the parent of their header fields (set by option 
%     'newhdr' in the binary function call).  By default 'leven' is set to 
%     'error'.
%     
%     BINOPERR('iftype','error'|'warn'|'ignore') sets how binary operations
%     handle records of different types.  If the option is set to 'warn' or
%     'ignore', the records are just worked on point for point.  The 
%     resultant records' iftype fields are determined by the parent of 
%     their header fields (set by option 'newhdr' in the binary function
%     call).  By default 'iftype' is set to 'warn'.
%
%    Notes:
%     - multiple options may be strung together in a single command
%
%    Examples:
%     % Turn off warnings for different timing:
%     binoperr('begin','ignore','ref','ignore')
%
%     % Show binary error settings:
%     binoperr
%
%    See also: ADDRECORDS, DIVIDERECORDS, MULTIPLYRECORDS, SUBTRACTRECORDS
%              MULTIFUN, CHECKOPERR

%     Version History:
%        June 20, 2008 - initial version
%        June 28, 2008 - doc update
%        Oct.  6, 2008 - doc update, use new SEIZMO layout
%        Nov. 22, 2008 - update for new name schema
%        Apr. 23, 2009 - move usage up
%        June  4, 2009 - minor doc fixes
%        June 28, 2009 - added 'truncate' and 'pad' states to 'npts' and
%                        'ncmp' options, some code reworking to accomadate
%        Aug. 21, 2009 - changed IFTYPE from ERROR to WARN to allow working
%                        with mixed xy and timeseries data
%        Oct. 21, 2009 - dropped rmfield usage (slow)
%        Aug. 16, 2010 - fprintf fix
%        Jan. 28, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2012 at 07:25 GMT

% todo:

% default options
option.NPTS='ERROR';
option.DELTA='ERROR';
option.BEGIN='WARN';
option.REF='WARN';
option.NCMP='ERROR';
option.LEVEN='ERROR';
option.IFTYPE='WARN';

% available states
valid.NPTS={'ERROR' 'WARN' 'TRUNCATE' 'PAD' 'IGNORE'};
valid.NCMP={'ERROR' 'WARN' 'TRUNCATE' 'PAD' 'IGNORE'};
valid.REF={'ERROR' 'WARN' 'IGNORE'};
valid.DELTA={'ERROR' 'WARN' 'IGNORE'};
valid.BEGIN={'ERROR' 'WARN' 'IGNORE'};
valid.LEVEN={'ERROR' 'WARN' 'IGNORE'};
valid.IFTYPE={'ERROR' 'WARN' 'IGNORE'};

% pull up and check over SEIZMO global
global SEIZMO; fields=fieldnames(option).';
if(isfield(SEIZMO,'BINOPERR'))
    for i=fields
        j=i{:};
        if(isfield(SEIZMO.BINOPERR,j))
            if(~any(strcmpi(SEIZMO.BINOPERR.(j),valid.(j))))
                warning('seizmo:binoperr:badState',...
                    '%s in unknown state => changing to default!',j);
                SEIZMO.BINOPERR.(i{:})=option.(j);
            else
                option.(j)=upper(SEIZMO.BINOPERR.(j));
            end
        end
    end
end

% no inputs = display options
if(nargin==0)
    for i=fields; fprintf('%12s = %s\n',i{:},option.(i{:})); end
% one input = 'defaults'
elseif(nargin==1)
    % check input is 'defaults'
    if(~strcmpi({'defaults' 'default' 'def'},varargin{:}))
        error('seizmo:binoperr:unknownOption',...
            'Unknown option or bad option usage!');
    end
    % clear SEIZMO settings
    if(isfield(SEIZMO,'BINOPERR'))
        SEIZMO.BINOPERR=[];
    end
% set options
else
    % all inputs must be 'option','state' pairs
    if(mod(nargin,2))
        error('seizmo:binoperr:unpairedOption','Option missing a value!');
    end
    % must be valid
    varargin=upper(varargin);
    for i=1:2:numel(varargin)
        if(~isfield(option,varargin{i}))
            error('seizmo:binoperr:unknownOption',...
                'Unknown option: %s',varargin{i});
        end
        
        if(~any(strcmpi(valid.(varargin{i}),varargin{i+1})))
            error('seizmo:binoperr:unknownOptionState',...
                'Option: %s\nUnknown state: %s',...
                varargin{i},varargin{i+1});
        end
    end
    % assign settings
    for i=1:2:nargin
        SEIZMO.BINOPERR.(varargin{i})=varargin{i+1};
    end
end

end
