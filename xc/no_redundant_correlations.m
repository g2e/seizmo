function [xc]=no_redundant_correlations(xc)
%NO_REDUNDANT_CORRELATIONS    Removes redundant correlations from dataset
%
%    Usage:    xc=no_redundant_correlations(xc)
%
%    Description:
%     XC=NO_REDUNDANT_CORRELATIONS(XC) removes any correlograms that are
%     redundant because they can be created from a REVERSE_CORRELATIONS
%     call on another correlogram in the same dataset.  Redundancies are
%     decided based on header info, in particular on the master & slave
%     name fields.
%
%    Notes:
%     - This does not check the relative nor absolute timing and so will
%       potentially eliminate correlograms that could be stacked and
%       correlograms that represent different lag time ranges (positive &
%       negative "lag component" correlations come to mind).
%     - This does not check component orientation, just the component code.
%
%    Examples:
%     % The following two commands should give equivalent datasets (the
%     % correlograms may be ordered slightly differently though):
%     xc1=correlate(data,'mcxc');
%     xc2=no_redundant_correlations(correlate(data,data,'mcxc'));
%
%    See also: CORRELATE, REVERSE_CORRELATIONS, ROTATE_CORRELATIONS, ISXC,
%              SPLIT_AUTO_CORRELATIONS, HORZ_CORRELATIONS_SETS,
%              NAME_CORRELATIONS, IS_FULL_MATRIX_OF_CORRELATIONS

%     Version History:
%        Sep.  9, 2013 - initial version
%        Sep. 19, 2013 - drop relative timing constraint, optimized
%                        checking, improved docs
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
nrecs=numel(xc);

% require all correlograms
if(~all(strcmp(kuser(:,1),'MASTER') & strcmp(kuser(:,2),'SLAVE')))
    error('seizmo:no_redundant_correlations:badInput',...
        'Some records appear not to be correlations!');
end

% full station names
mnm=lower(strcat(mnm(:,1),'.',mnm(:,2),'.',mnm(:,3),'.',mnm(:,4)));
snm=lower(strcat(snm(:,1),'.',snm(:,2),'.',snm(:,3),'.',snm(:,4)));

% division into 2 set names to handle master<=>slave ambiguity
setname1=strcat(mnm,'_',snm);
setname2=strcat(snm,'_',mnm);
[setname,r2sidx,s2ridx]=unique([setname1;setname2],'first');

% number of sets
nsets=numel(setname);

% double over set to record indexing so the 2 sets
% each record is in are easy to visualize/index/link
s2ridx=reshape(s2ridx,nrecs,2);

% annihilation loop
keep=false(nrecs,1);
settaken=false(nsets,1);
for i=1:nrecs
    % is this set taken?
    if(~settaken(s2ridx(i,1)))
        % nope, so keep this record ...
        keep(i)=true;
        
        % ... and nobody else in this set
        settaken(s2ridx(i,:))=true;
    end
end

% annihilate redundancies
xc(~keep)=[];

end
