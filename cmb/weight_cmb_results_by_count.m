function [results]=weight_cmb_results_by_count(results)
%WEIGHT_CMB_RESULTS_BY_COUNT    Weight cmb results for profile measurements
%
%    Usage:    results=weight_cmb_results_by_count(results)
%
%    Description:
%     RESULTS=WEIGHT_CMB_RESULTS_BY_COUNT(RESULTS) adds a .weights field
%     to RESULTS that is used by SLOWDECAYPROFILES when calculating the
%     rayparameter and decay constants.  Weights are determined by a simple
%     station count.  This emphasizes the data from stations that are more
%     common while de-emphasizing data from stations that are less common
%     (like due to station glitches).  For example if RESULTS has 25
%     elements (say from running CMB_2ND_PASS), a station found in only 3
%     of those would be given a weight of 3.  Weights are normalized to sum
%     to 1 (for each RESULTS element).  Outliers are not included in the
%     count.
%
%    Notes:
%     - Uses the station codes in the KNAME group field for discrimination.
%     - If you rerun CMB_OUTLIERS after WEIGHT_CMB_RESULTS_BY_COUNT then
%       SLOWDECAYPROFILES will probably throw an error.  Re-run this
%       function after CMB_OUTLIERS to account for changes in the outliers.
%
%    Examples:
%     % The main intention of this function is to help smooth dispersion
%     % curves that are disrupted by stations not included in every
%     % frequency.  This is my usual way of including this function:
%     result2=cmb_2nd_pass(results,5);
%     for i=1:size(results2,2)
%         result2(:,i)=weight_cmb_results_by_count(results2(:,i));
%     end
%     pf=slowdecayprofiles(results2);
%
%    See also: SLOWDECAYPROFILES, CMB_2ND_PASS, CMB_OUTLIERS

%     Version History:
%        Mar.  6, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2012 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check results struct
error(check_cmb_results(results));

% clear out empty results
g=find(~cellfun('isempty',{results.useralign}));
nr=numel(results(g));

% loop over results and get names
kname=cell(nr,1);
for i=1:nr
    n=getheader(results(g(i)).useralign.data(...
        ~results(g(i)).outliers.bad),'kname');
    kname{i,1}=strcat(n(:,1),'.',n(:,2),'.',n(:,3),'.',n(:,4));
end

% get the count
[names,i,j]=unique(cat(1,kname{:}));
n=histc(j,1:max(j));

% set weights
for i=1:nr
    [tf,loc]=ismember(kname{i,1},names);
    results(g(i)).weights=zeros(size(kname));
    results(g(i)).weights(tf)=n(loc(tf))./sum(n(loc(tf)));
end

end
