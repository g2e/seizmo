function [dep,idx1,ind,idx2,store,npts]=records2mat(data)
%RECORDS2MAT    Combines SEIZMO data records into a single numeric matrix
%
%    Usage:    [dep,idx1,ind,idx2,store,npts]=records2mat(data)
%
%    Description:
%     [DEP,IDX1,IND,IDX2,STORE,NPTS]=RECORDS2MAT(DATA) combines SEIZMO
%     records stored in DATA into numeric arrays DEP & IND (components are
%     in separate columns).  IDX1 & IDX2 indicate which columns belong to
%     which records in DATA.  STORE is a cell array of data class strings,
%     one for each record.  NPTS gives the original npts in each record
%     (records are padded with zeros when combined).  This function is
%     useful for providing easy access to functions not in the SEIZMO
%     toolbox.  Use MAT2RECORDS to redistribute the records back into DATA
%     if necessary.
%    
%    Notes:
%
%    Header changes: see CHECKHEADER
%
%    Examples:
%     % Get interquartile range of records and assign to a header field:
%     dep=records2mat(data);
%     data=changeheader(data,'user3',iqr(dep));
%
%    See also: MAT2RECORDS, BSEIZMO, GETVALUEFUN, SOLOFUN, MULTIFUN,
%              SLIDINGFUN

%     Version History:
%        Feb. 16, 2008 - initial version
%        Feb. 21, 2008 - minor doc update
%        Feb. 23, 2008 - minor doc update
%        Feb. 28, 2008 - seischk support
%        Mar.  4, 2008 - minor doc update
%        June 15, 2008 - doc update
%        June 20, 2008 - minor doc update
%        June 29, 2008 - doc update, .dep rather than .x,
%                        dataless support
%        Nov. 22, 2008 - update for new name schema (now COMBINERECORDS)
%                        now outputs independent component
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        June 12, 2009 - minor doc update
%        June 25, 2009 - name change from COMBINERECORDS to RECORDS2MAT
%        Jan. 30, 2010 - proper SEIZMO handling, seizmoverbose support
%        Feb.  3, 2010 - versioninfo caching
%        Feb. 11, 2011 - mass nargchk fix, see also section update, 
%                        dropped versioninfo caching
%        Mar. 24, 2012 - minor doc update
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:05 GMT

% todo:

% check input
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt matrix creation
try
    % check headers
    data=checkheader(data);
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);

    % get header fields
    [b,npts,delta]=getheader(data,'b','npts','delta');
    
    % detail message
    if(verbose)
        disp('Exporting Record(s) as Matrices');
        print_time_left(0,nrecs);
    end

    % loop through records
    store=cell(nrecs,1);
    npts=zeros(nrecs,1);
    ncol=zeros(nrecs,1);
    for i=1:nrecs
        store{i}=class(data(i).dep);
        [npts(i),ncol(i)]=size(data(i).dep);
        if(npts(i)*ncol(i)==0); npts(i)=0; ncol(i)=0; end
    end
    leven=true(nrecs,1);
    if(isfield(data,'ind'))
        for i=1:nrecs
            if(~isempty(data(i).ind))
                leven(i)=false;
            end
        end
    end

    % preallocate dep matrix (as double precision)
    col2=cumsum(ncol); col=[1; col2+1];
    dep=zeros(max(npts),col2(end));
    idx1=zeros(1,col2(end));

    % preallocate ind matrix
    ind=zeros(max(npts),nrecs);
    idx2=1:nrecs;

    % loop through records
    for i=1:nrecs
        idx1(col(i):col2(i))=i;
        dep(1:npts(i),col(i):col2(i))=data(i).dep;
        if(leven(i))
            ind(1:npts(i),i)=b(i)+(0:delta(i):delta(i)*(npts(i)-1));
        else
            ind(1:npts(i),i)=data(i).ind;
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
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
