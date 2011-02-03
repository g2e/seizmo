function [results]=cmb_clustering(results,odir,figdir)
%CMB_CLUSTERING    Cluster analysis of core-diffracted data
%
%    Usage:    results=cmb_clustering(results)
%              results=cmb_clustering(results,odir)
%              results=cmb_clustering(results,odir,figdir)
%
%    Description:
%     RESULTS=CMB_CLUSTERING(RESULTS) provides a graphical interface for
%     clustering the core-diffracted wave analysis RESULTS struct.  This is
%     intended to aid in removing noisy waveforms and isolating dissimilar
%     waveform sets from one another.  RESULTS is the output from either
%     CMB_1ST_PASS, CMB_OUTLIERS or CMB_2ND_PASS.  Plots and the modified
%     RESULTS struct are saved to the current directory.
%
%     RESULTS=CMB_CLUSTERING(RESULTS,ODIR) sets the output directory
%     where the figures and RESULTS struct is saved.  By default ODIR is
%     '.' (the current directory.
%
%     RESULTS=CMB_CLUSTERING(RESULTS,ODIR,FIGDIR) allows saving figures to
%     a different directory than ODIR (where the RESULTS struct is saved).
%
%    Notes:
%     - Redoing the alignment (from CMB_1ST_PASS) is probably necessary if
%       you specify that some records need their ground units changed.
%     - While outliers specified using CMB_OUTLIERS are ignored here, they
%       are not disturbed -- so they remain outliers even though their
%       cluster may have changed.  Running CMB_OUTLIERS will reset the
%       outliers each time it is ran.
%     - Clustering is reset each time CMB_CLUSTERING is ran on a RESULTS
%       struct.
%
%    Examples:
%     % Typical alignment and refinement workflow:
%     results=cmb_1st_pass;
%     results=cmb_clustering(results);
%     results=cmb_outliers(results);
%
%    See also: PREP_CMB_DATA, CMB_1ST_PASS, CMB_2ND_PASS, SLOWDECAYPAIRS,
%              SLOWDECAYPROFILES, MAP_CMB_PROFILES, CMB_OUTLIERS

%     Version History:
%        Jan. 15, 2011 - initial version
%        Jan. 16, 2011 - fix for results standardization
%        Jan. 18, 2011 - .time field
%        Jan. 29, 2011 - prepend datetime to output names
%        Jan. 31, 2011 - odir & figdir inputs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 31, 2011 at 13:35 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check results
error(check_cmb_results(results));

% default odir & figdir
if(nargin<2 || isempty(odir)); odir='.'; end
if(nargin<3 || isempty(figdir)); figdir=odir; end

% check odir & figdir
if(~isstring(odir))
    error('seizmo:cmb_clustering:badInput',...
        'ODIR must be a string!');
elseif(~isstring(figdir))
    error('seizmo:cmb_clustering:badInput',...
        'FIGDIR must be a string!');
end

% make sure odir/figdir exists (create it if it does not)
[ok,msg,msgid]=mkdir(odir);
if(~ok)
    warning(msgid,msg);
    error('seizmo:cmb_clustering:pathBad',...
        'Cannot create directory: %s',odir);
end
[ok,msg,msgid]=mkdir(figdir);
if(~ok)
    warning(msgid,msg);
    error('seizmo:cmb_clustering:pathBad',...
        'Cannot create directory: %s',figdir);
end

% loop over each event
for i=1:numel(results)
    % display name
    disp(results(i).runname);
    
    % abandon events we skipped
    if(isempty(results(i).useralign)); continue; end
    
    % cluster analysis
    [results(i).usercluster,ax]=usercluster(results(i).useralign.data,...
        results(i).useralign.xc.cg(:,:,1),...
        [],[],[],[],'normstyle','single');
    if(any(ishandle(ax)))
        fh=unique(cell2mat(get(ax(ishandle(ax)),'parent')));
        for j=1:numel(fh)
            saveas(fh(j),fullfile(figdir,[datestr(now,30) '_' ...
                results(i).runname '_usercluster_' num2str(j) '.fig']));
            close(fh(j));
        end
    end
    
    % advanced clustering
    [results(i).useralign.data,results(i).usercluster,...
        results(i).useralign.solution.arr,...
        results(i).useralign.solution.pol,...
        results(i).usercluster.units]=adjustclusters(...
        results(i).useralign.data,results(i).usercluster,...
        results(i).useralign.solution.arr,...
        results(i).useralign.solution.pol);
    
    % time
    results(i).time=datestr(now);

    % save results
    tmp=results(i);
    save(fullfile(figdir,[datestr(now,30) '_' results(i).runname ...
        '_clustering_results.mat']),'-struct','tmp');
end

end
