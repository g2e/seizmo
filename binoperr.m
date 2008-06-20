function []=binoperr(varargin)
%BINOPERR    Controls behavior of SAClab binary functions
%
%    Description: Allows for changing the behavior of how binary functions 
%     (addf, divf, mulf, divf) handle records with some unequal aspects.
%     In particular, this function will modify the behavior of all of these
%     functions until the Matlab session is finished or another call to
%     BINOPERR undoes it.
%
%     BINOPERR() displays the current binary operator error settings.
%
%     BINOPERR('defaults') clears any previous settings and causes all
%     SAClab binary functions to use their default settings.
%     
%     BINOPERR('npts','error|warn|ignore') sets how binary operations
%     handle records with different numbers of points.  If the option is 
%     set to 'warn' or 'ignore', the number of points in the resultant 
%     records will be equal to that of the shortest record.  Note that 
%     binary operations work point by point (ignore timing).  By default 
%     'npts' is set to 'error'.
%     
%     BINOPERR('delta','error|warn|ignore') sets how binary operations
%     handle records with different sample rates.  If the option is set to
%     'warn' or 'ignore', the records are just operated on point for point
%     (basically ignoring timing).  The resultant records' sample rates are
%     determined by the parent of their header fields (set by the option 
%     'newhdr' in the binary function call).  By default 'delta' is set to
%     'error'.
%     
%     BINOPERR('begin','error|warn|ignore') sets how binary operations
%     handle records with different begin times.  If the option is set to
%     'warn' or 'ignore', the resultant records' begin times are 
%     determined by the parent of their header fields (set by option 
%     'newhdr' in the binary function call).  By default 'begin' is set to 
%     'warn'.
%     
%     BINOPERR('ref','error|warn|ignore') sets how binary operations handle
%     records with different reference times.  If the option is set to 
%     'warn' or 'ignore', the resultant records' reference times are 
%     determined by the parent of their header fields (set by option 
%     'newhdr' in the binary function call).  By default 'ref' is set to 
%     'warn'.
%     
%     BINOPERR('ncmp','error|warn|ignore') sets how binary operations 
%     handle records with different numbers of components.  If the option
%     is set to 'warn' or 'ignore', the number of components in the
%     resultant records will be equal to that of the record with the least.  
%     Note that components are operated on according to their order in the 
%     record so that the first components always go together.  By default 
%     'ncmp' is set to 'error'.
%     
%     BINOPERR('leven','error|warn|ignore') sets how binary operations
%     handle unevenly sampled records.  If the option is set to 'warn' or 
%     'ignore', the functions behave as normal working point for point 
%     (basically ignoring timing).  The resultant records' leven fields are
%     determined by the parent of their header fields (set by option 
%     'newhdr' in the binary function call).  By default 'leven' is set to 
%     'error'.
%     
%     BINOPERR('iftype','error|warn|ignore') sets how binary operations
%     handle records of different types.  If the option is set to 'warn' or
%     'ignore', the records are just worked on point for point.  The 
%     resultant records' iftype fields are determined by the parent of 
%     their header fields (set by option 'newhdr' in the binary function
%     call).  By default 'iftype' is set to 'error'.
%
%    Notes:
%     - options may be strung together
%
%    System requirements: Matlab 7
%
%    Data requirements: Inputs must all be strings.
%
%    Header Changes: NONE
%
%    Usage: binoperr('option1','error|warn|ignore',...,
%                    'optionN','error|warn|ignore')
%
%    See also: addf, divf, mulf, subf

%     Version History:
%        June 20, 2008 - Initial Version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 20, 2008 at 04:25 GMT

% default options
option.NPTS='ERROR';
option.DELTA='ERROR';
option.BEGIN='WARN';
option.REF='WARN';
option.NCMP='ERROR';
option.LEVEN='ERROR';
option.IFTYPE='ERROR';

% available states
states={'ERROR' 'WARN' 'IGNORE'};

% pull up and check over SACLAB global
global SACLAB; fields=fieldnames(option).';
for i=fields 
    if(isfield(SACLAB,i))
        if(~any(strcmpi(SACLAB.(i{:}),states)))
            warning('SAClab:binoperr:badState',...
                '%s in unknown state => changing to default',i{:})
            SACLAB.(i{:})=option.(i{:});
        else
            option.(i{:})=upper(SACLAB.(i{:}));
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
        error('SAClab:binoperr:unknownOption',...
            'Unknown option or bad option usage!')
    end
    % clear SACLAB settings
    for i=fields; if(isfield(SACLAB,i)); SACLAB=rmfield(SACLAB,i); end; end
% set options
else
    % all inputs must be 'option','state' pairs
    if(mod(nargin,2))
        error('SAClab:binoperr:unpairedOption','Option missing a value!')
    end
    % must be valid
    varargin=upper(varargin);
    for i=varargin(1:2:end)
        if(~isfield(option,i))
            error('SAClab:binoperr:unknownOption',...
                'Unknown option: %s',i{:})
        end
    end
    for i=varargin(2:2:end)
        if(~any(strcmpi(i,states)))
            error('SAClab:binoperr:unknownOptionState',...
                'Unknown option state: %s',i{:})
        end
    end
    % assign settings
    for i=1:2:nargin
        SACLAB.(varargin{i})=varargin{i+1};
    end
end 

end
