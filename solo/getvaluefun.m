function [value]=getvaluefun(data,func,type,scalar)
%GETVALUEFUN    Applies a function to SEIZMO records and returns the value
%
%    Usage:    value=getvaluefun(data,func)
%              value=getvaluefun(data,func,cmptype)
%              value=getvaluefun(data,func,cmptype,scalarOutput)
%
%    Description:
%     VALUE=GETVALUEFUN(DATA,FUNC) applies the function handle FUNC to the
%     dependent component of records in SEIZMO struct DATA.  The function
%     FUNC is expected to return a scalar value.  The values are returned
%     in VALUE (a Nx1 column vector where N is the number of records in
%     DATA).
%
%     VALUE=GETVALUEFUN(DATA,FUNC,CMPTYPE) specifies whether FUNC is
%     applied to the dependent or independent data using CMPTYPE.  CMPTYPE
%     must be either 'dep' or 'ind'.  The default value is 'dep'.  An empty
%     value ([]) also will specify the default value.
%
%     VALUE=GETVALUEFUN(DATA,FUNC,CMPTYPE,SCALAROUTPUT) allows returning
%     non-scalar output in VALUE when SCALAROUTPUT is set FALSE.  VALUE
%     will be a Nx1 cell array in this case.  The default SCALAROUTPUT is
%     TRUE.
%
%    Notes:
%     - Both .dep & .ind data fields are always passed in double precision.
%
%    Examples:
%     % Get dependent data medians:
%     medians=getvaluefun(data,@median);
%
%     % If you have multi-cmp records the following may be necessary:
%     medians=getvaluefun(data,@(x)median(x(:)));
%     % or (if you want the median for each component):
%     medians=getvaluefun(data,@median,[],false);
%
%     % Get RMS:
%     rms=getvaluefun(data,@(x)sqrt(mean(x.^2));
%
%     % Get robust RMS:
%     rms=getvaluefun(data,@(x)sqrt(median(x.^2));
%
%     % Get maximum amplitude of multi-cmp records assuming each component
%     % is orthogonal (was an old function called getnorm):
%     normalizers=getvaluefun(data,@(x)max(sqrt(sum(x.^2,2))));
%
%    See also: SOLOFUN, SLIDINGFUN, MULTIFUN

%     Version History:
%        Mar. 18, 2010 - initial version
%        Mar. 20, 2010 - fixed bug that added extra time point
%        Mar. 26, 2010 - added example to copy GETNORM (now deprecated)
%        May  10, 2010 - better SCALAROUTPUT implementation
%        Jan.  6, 2011 - nargchk fix, seizmofun/solofun rename,
%                        recordfun/multifun rename
%        Jan. 12, 2012 - minor doc update
%        Mar. 13, 2012 - seizmocheck fix, use getheader improvements
%        June 11, 2012 - fix rms & robust rms examples
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 18:25 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt to apply function
try
    % verbosity
    verbose=seizmoverbose;
    
    % defaults
    if(nargin<3 || isempty(type)); type='dep'; end
    if(nargin<4 || isempty(scalar)); scalar=true; end

    % check other arguments
    if(~isa(func,'function_handle'))
        error('seizmo:getvaluefun:badInput',...
            'FUNC must be a function handle!');
    elseif(~ischar(type) || size(type,1)~=1 || ...
            ~any(strcmpi(type,{'dep' 'ind'})))
        error('seizmo:getvaluefun:badInput',...
            'CMPTYPE must be one of the following:\n''DEP'' or ''IND''');
    elseif(~isscalar(scalar) || ~islogical(scalar))
        error('seizmo:getvaluefun:badInput',...
            'SCALAROUTPUT must be TRUE or FALSE!');
    end

    % number of records
    nrecs=numel(data);

    % preallocate value
    value=cell(nrecs,1);

    % which cmp are we working on?
    type=lower(type);
    switch type
        case 'dep'
            % detail message
            if(verbose)
                disp('Getting Value From Dependent Data of Record(s)');
                print_time_left(0,nrecs);
            end
            
            % apply function
            for i=1:nrecs
                value{i}=func(double(data(i).dep));
                
                % detail message
                if(verbose); print_time_left(i,nrecs); end
            end
        case 'ind'
            % fill .ind for evenly spaced arrays
            leven=~strcmpi(getheader(data,'leven lgc'),'false');
            
            % are any are evenly spaced?
            if(any(leven))
                % pull header values
                [b,delta,npts]=getheader(data(leven),'b','delta','npts');
                
                % loop over even, add .ind
                idx=find(leven);
                for i=1:sum(leven)
                    data(idx(i)).ind=...
                        b(i)+(0:delta(i):delta(i)*(npts(i)-1)).';
                end
            end
            
            % detail message
            if(verbose)
                disp('Getting Value From Independent Data of Record(s)');
                print_time_left(0,nrecs);
            end

            % apply function
            for i=1:nrecs
                value{i}=func(double(data(i).ind));
                
                % detail message
                if(verbose); print_time_left(i,nrecs); end
            end
    end
    
    % error if not scalar
    if(scalar)
        if(any(cellfun('prodofsize',value)~=1))
            error('seizmo:getvaluefun:nonScalarOutput',...
                ['Non-scalar output when scalar output expected!\n' ...
                'Please set SCALAROUTPUT option to FALSE!']);
        else
            value=cell2mat(value);
        end
    end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);

    % rethrow error
    error(lasterror);
end

end
