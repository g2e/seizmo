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
%     '.' (the current directory.  You may set ODIR to FALSE for no
%     written output.
%
%     RESULTS=CMB_CLUSTERING(RESULTS,ODIR,FIGDIR) allows saving figures to
%     a different directory than ODIR (where the RESULTS struct is saved).
%     You may set FIGDIR to FALSE to skip saving figures.
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
%        Feb. 11, 2011 - results now output to odir not figdir
%        Apr.  1, 2011 - update .time field of skipped
%        Mar.  1, 2012 - octave ascii save workaround
%        Mar.  5, 2012 - allow no written output
%        Mar. 11, 2013 - directory input (reads indir/*.mat), selection
%                        list, advanced clustering commented out
%        Jan. 27, 2014 - abs path fix & reduced filesep/fullfile calls
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 13:35 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% directory separator
fs=filesep;

% handle directory input
if(isstring(results))
    if(~isabspath(results)); results=[pwd fs results]; end
    if(isdir(results))
        files=xdir([results fs '*.mat']);
        clear results;
        for i=1:numel(files)
            results(i)=load([files(i).path files(i).name]);
        end
    end
end

% check results
error(check_cmb_results(results));

% default odir & figdir
if(nargin<2 || isempty(odir)); odir='.'; end
if(islogical(odir) && isscalar(odir) && odir); odir='.'; end
if(nargin<3 || isempty(figdir)); figdir=odir; end
if(islogical(figdir) && isscalar(figdir) && figdir); figdir='.'; end

% check odir & figdir
if(~isstring(odir) && ~(islogical(odir) && isscalar(odir)))
    error('seizmo:cmb_clustering:badInput',...
        'ODIR must be a string or TRUE/FALSE!');
elseif(~isstring(figdir) && ~(islogical(figdir) && isscalar(figdir)))
    error('seizmo:cmb_clustering:badInput',...
        'FIGDIR must be a string or TRUE/FALSE!');
end
if(islogical(odir)); out=odir; else out=true; end
if(islogical(figdir)); figout=odir; else figout=true; end

% make sure odir/figdir exists (create it if it does not)
if(out)
    [ok,msg,msgid]=mkdir(odir);
    if(~ok)
        warning(msgid,msg);
        error('seizmo:cmb_clustering:pathBad',...
            'Cannot create directory: %s',odir);
    end
elseif(~out && ~nargout)
    error('seizmo:cmb_clustering:badInput',...
        'Output variable must be assigned when no written output!');
end
if(figout)
    [ok,msg,msgid]=mkdir(figdir);
    if(~ok)
        warning(msgid,msg);
        error('seizmo:cmb_clustering:pathBad',...
            'Cannot create directory: %s',figdir);
    end
end

% select events
datelist=char({results.runname}.');
s=listdlg('PromptString','Select events:',...
          'InitialValue',1:numel(results),...
          'ListSize',[170 300],...
          'ListString',datelist);

% error if none selected
if(isempty(s))
    error('seizmo:cmb_clustering:noDirsSelected',...
        'No earthquakes selected!');
end

% loop over each event
for i=1:numel(s)
    % display name
    disp(results(s(i)).runname);
    
    % time (for skipped)
    results(s(i)).time=datestr(now);
    
    % abandon events we skipped
    if(isempty(results(s(i)).useralign)); continue; end
    
    % cluster analysis
    [results(s(i)).usercluster,ax]=usercluster(...
        results(s(i)).useralign.data,...
        results(s(i)).useralign.xc.cg(:,:,1),...
        [],[],[],[],'normstyle','single');
    if(any(ishandle(ax)))
        fh=unique(cell2mat(get(ax(ishandle(ax)),'parent')));
        for j=1:numel(fh)
            if(figout)
                saveas(fh(j),[figdir fs datestr(now,30) '_' ...
                    results(s(i)).runname '_usercluster_' num2str(j) ...
                    '.fig']);
            end
            close(fh(j));
        end
    end
    
    % advanced clustering
    %[results(s(i)).useralign.data,results(s(i)).usercluster,...
    %    results(s(i)).useralign.solution.arr,...
    %    results(s(i)).useralign.solution.pol,...
    %    results(s(i)).usercluster.units]=adjustclusters(...
    %    results(s(i)).useralign.data,results(s(i)).usercluster,...
    %    results(s(i)).useralign.solution.arr,...
    %    results(s(i)).useralign.solution.pol);
    
    % time
    results(s(i)).time=datestr(now);

    % save results
    if(out)
        tmp=results(s(i));
        if(isoctave)
            save([odir fs datestr(now,30) '_' results(s(i)).runname ...
                '_clustering_results.mat'],'-7','-struct','tmp');
        else % matlab
            save([odir fs datestr(now,30) '_' results(s(i)).runname ...
                '_clustering_results.mat'],'-struct','tmp');
        end
    end
end

end
