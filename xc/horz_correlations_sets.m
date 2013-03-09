function [idx1,idx2,idx3,reverse]=horz_correlations_sets(xc)
%HORZ_CORRELATIONS_SETS   Returns indices for horiz. correlation sets
%
%    Usage:    [idx1,idx2,idx3]=horz_correlations_sets(xc)
%
%    Description:
%     [IDX1,IDX2,IDX3]=HORZ_CORRELATION_SETS(XC) groups horizontal
%     correlograms into rotatible sets for ROTATE_CORRELATIONS.  This is
%     similar to the functions HORZPAIRS & FINDTRIPLETS which find the
%     rotatible data for the functions ROTATE & ROTATE3.  XC is a SEIZMO
%     struct of correlograms created by CORRELATE.  Sets are identified by
%     having common stations, horizontal orientation, and common lag times
%     and will always come in sets of 3 (for autocorr) or 4.  Components
%     are ordered as either:
%          (1) NN, NE, EN, EE
%          (2) RR, RT, TR, TT.
%     depending on the input.  Any other orientation is not useful and will
%     generate an error.  IDX1 gives the indices of the correlograms in XC
%     that are in a set.  IDX2 gives the set indices (use max(IDX2) to get
%     the number of sets). IDX3 gives the component indice (1, 2, 3, or 4).
%
%    Notes:
%
%    Examples:
%     % Remove vertical correlations and non-set horizontal correlations:
%     xc=xc(horz_correlation_sets(xc));
%
%    See also: CORRELATE, SPLIT_HORZ_CORRELATIONS, ROTATE_CORRELATIONS,
%              REVERSE_CORRELATIONS, SPLIT_AUTO_CORRELATIONS, HORZPAIRS,
%              FINDTRIPLETS



% todo:
% - auto = 3 files
% - cross = 4 files

% check nargin
error(nargchk(1,1,nargin));

% check structure
error(seizmocheck(xc));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    xc=checkheader(xc,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% header info
[kuser,scmp,b,e,st,ev,npts,delta,stnm,evnm,ecmp,gcp,az]=getheader(xc,...
    'kuser','cmp','b','e','st','ev','npts','delta','kname','kt','user',...
    'gcp','az');
mi=ecmp(:,1); si=ecmp(:,2); ecmp=ecmp(:,3:4);
mc=char(evnm(:,4)); sc=char(stnm(:,4));

% error if non-correlations
if(~all(strcmp(kuser(:,1),'MASTER') & strcmp(kuser(:,2),'SLAVE')))
    error('seizmo:horz_correlation_sets:badInput',...
        'Some records appear not to be correlations!');
end

% retain only horizontal correlations
horz=ecmp(:,1)==90 & scmp(:,1)==90;
if(sum(horz)<4); [idx1,idx2,idx3,reverse]=deal([]); return; end



end
