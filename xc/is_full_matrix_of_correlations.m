function [lgc,missing]=is_full_matrix_of_correlations(xc)
%IS_FULL_MATRIX_OF_CORRELATIONS    Are all correlation pairs present?
%
%    Usage:    lgc=is_full_matrix_of_correlations(xc)
%              [lgc,missing]=is_full_matrix_of_correlations(xc)
%
%    Description:
%     LGC=IS_FULL_MATRIX_OF_CORRELATIONS(XC) returns TRUE if the SEIZMO
%     dataset XC contains all possible pairings between the stations.  This
%     requirement includes autocorrelations.  Correlations that are
%     redundant (can be made by REVERSE_CORRELATIONS on records in XC) are
%     not required.
%
%     [LGC,REV]=IS_FULL_MATRIX_OF_CORRELATIONS(XC) returns the indices of
%     the correlations that need to be reversed to complete the full
%     correlation matrix (using REVERSE_CORRELATIONS) if and only if LGC is
%     TRUE (see the next Usage form if LGC is FALSE).  REV is a list of
%     indices corresponding to correlograms in XC.
%
%     [LGC,MISSING]=IS_FULL_MATRIX_OF_CORRELATIONS(XC) returns the missing
%     station pairs needed to allow for a full correlation matrix if LGC is
%     FALSE.  This list does not include redundant station pairings.
%     MISSING is a Nx2 cell array of full station-component names as
%     strings ('NET.STN.HOLE.CMP') for identifying the N missing station
%     pairs (e.g., MISSING(1,:) gives the 1st missing pair).
%
%    Notes:
%
%    Examples:
%     % Both return TRUE but the first requires reversing some records:
%     [lgc,rev]=is_full_matrix_of_correlations(correlate(data,'mcxc'))
%     [lgc,rev]=is_full_matrix_of_correlations(correlate(data,data,'mcxc'))
%
%    See also: CORRELATE, REVERSE_CORRELATIONS, ROTATE_CORRELATIONS, ISXC,
%              SPLIT_AUTO_CORRELATIONS, HORZ_CORRELATIONS_SETS,
%              NAME_CORRELATIONS, NO_REDUNDANT_CORRELATIONS

%     Version History:
%        Sep. 19, 2013 - initial version
%        Jan. 15, 2014 - fixed warning id
%        June 12, 2014 - handle fd i/o
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 12, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check structure
error(seizmocheck(xc,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% safely check headers & grab necessary info
try
    % check headers
    xc=checkheader(xc,...
        'MULCMP_DEP','ERROR',...
        'XYZ_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    
    % necessary header info
    [kuser,mnm,snm]=getheader(xc,'kuser','kt','kname');
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% number of correlograms
nxc=numel(xc);

% require all correlograms
if(~all(strcmp(kuser(:,1),'MASTER') & strcmp(kuser(:,2),'SLAVE')))
    error('seizmo:is_full_matrix_of_correlations:badInput',...
        'Some records appear not to be correlations!');
end

% full station names
mnm=lower(strcat(mnm(:,1),'.',mnm(:,2),'.',mnm(:,3),'.',mnm(:,4)));
snm=lower(strcat(snm(:,1),'.',snm(:,2),'.',snm(:,3),'.',snm(:,4)));

% use unique to determine the matrix members and find the indices
[rnm,x2ridx,r2xidx]=unique([snm;mnm],'first');

% number of underlying records
nrecs=numel(rnm);

% double over indexing so that we now have row/column indices
r2xidx=reshape(r2xidx,nxc,2);
r2xlidx=sub2ind([nrecs nrecs],r2xidx(:,1),r2xidx(:,2)); % linear

% create logical matrix for indicating which correlations are present
% - defaults is none present
% - rows are slave index, columns are master index
xcmat=false(nrecs);

% indicate which are present by setting their indice to true
xcmat(r2xlidx)=true;

% allow the reverse correlation to count
xcmatr=xcmat | xcmat';

% are all correlations now accounted for?
lgc=all(xcmatr(:));

% shortcut if only one output
if(nargout==1); return; end

% second output
if(lgc) % return correlations that need to be reversed
    % transposed negative of correlation logical matrix gives
    % the pairs to be reversed...
    lind=find(~xcmat');
    
    % now match the linear indices to an input correlation record
    % - note: returns the last record in the case of multiples
    [tf,missing]=ismember(lind,r2xlidx);
else % return the station pair strings in cells
    % find the indices of those missing
    [s,m]=find(~xcmatr);
    
    % require slave index >= master to eliminate redundancies
    [m,s]=deal(m(s>=m),s(s>=m));
    
    % create output cell array
    missing=[rnm(m) rnm(s)];
end

end
