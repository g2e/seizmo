function [ax]=investigate_cmb_results(results,observable)
%INVESTIGATE_CMB_RESULTS    Plots results output from cmb functions
%
%    Usage:    investigate_cmb_results(results)
%              investigate_cmb_results(results,observable)
%              ax=investigate_cmb_results(...)
%
%    Description:
%     INVESTIGATE_CMB_RESULTS(RESULTS) creates a figure showing all
%     observables as well as rayparameter and decay constant.  This is
%     great for diagnosing data problems.
%
%     INVESTIGATE_CMB_RESULTS(RESULTS,OBSERVABLE) plots
%
%     AX=INVESTIGATE_CMB_RESULTS(...) returns the handle(s) of the axes
%     plotted in.  If no observable was specified 
%
%    Notes:
%     - Outliers are set to NaN in the observables plots.
%
%    Examples:
%     % Diagnose dispersion curve issues for an event:
%     results2=cmb_2nd_pass(results);
%     investigate_cmb_results(results2(:,1));
%
%    See also: CMB_1ST_PASS, CMB_2ND_PASS, CMB_OUTLIERS, CMB_CLUSTERING,
%              PLOT_CMB_PDF, PLOT_CMB_MEASUREMENTS

%     Version History:
%        Mar.  5, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  5, 2012 at 13:35 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check results struct
error(check_cmb_results(results));

% check observable
if(nargin>1 && ~isempty(observable))
    valid={'arr' 'err' 'amp' 'snr' 'all'};
    if(~ischar(observable) || size(observable,1)~=1 ...
            || ~any(strcmpi(observable,valid)))
        error('seizmo:investigate_cmb_results:badInput',...
            ['OBSERVABLE must be ''arr'', ''arrerr'', ''amp'', ' ...
            '''snr'', or ''all''!']);
    end
else
    observable='all';
end

% create figure
fh=figure('color','k','defaultaxesxcolor','w','defaultaxesycolor','w',...
    'defaulttextcolor','w');

% axes
if(strcmpi(observable,'all'))
    ax=makesubplots(2,3,1:6,'parent',fh);
    ax=reshape(ax,3,2).';
else
    ax=axes('color','none','parent',fh);
end

% clear out empty results
bad=cellfun('isempty',{results.useralign});
oi=find(~bad);
results(bad)=[];
nr=numel(results);

% single or multi-event?
% - single event has stations sorted by distance
% - multi-event has stations sorted by name
dirs={results.dirname};
if(numel(unique(dirs))>1)
    mevflag=true;
else
    mevflag=false;
end

% get observables for each result
[d,kname,x]=deal(cell(nr,1));
for i=1:numel(results)
    % station names & distances
    [n,d{i,1}]=getheader(results(i).useralign.data,'kname','gcarc');
    kname{i,1}=strcat(n(:,1),'.',n(:,2),'.',n(:,3),'.',n(:,4));
    
    % which observable
    switch lower(observable)
        case 'arr'
            % deviation from the expected travel time
            % rather than the total travel time
            verbose=seizmoverbose(false);
            x{i,1}=results(i).useralign.solution.arr ...
                - (getarrival(results(i).useralign.data,...
                {results(i).phase(1) results(i).phase}) ...
                - getheader(results(i).useralign.data,'o'));
            seizmoverbose(verbose);
            x{i,1}(results(i).outliers.bad)=nan;
        case 'err'
            % arrival time error
            x{i,1}=results(i).useralign.solution.arrerr;
            x{i,1}(results(i).outliers.bad)=nan;
        case 'snr'
            % measured snr
            snr=results(i).usersnr.snr;
            snr=snr(snr>=results(i).usersnr.snrcut);
            snr(results(i).userwinnow.cut)=[];
            snr=snr(results(i).finalcut);
            x{i,1}=log10(snr);
            x{i,1}(results(i).outliers.bad)=nan;
        case 'amp'
            % amplitude
            x{i,1}=log10(results(i).useralign.solution.amp);
            x{i,1}(results(i).outliers.bad)=nan;
        case 'all'
            % all of the above
            verbose=seizmoverbose(false);
            x{i,1}=results(i).useralign.solution.arr ...
                - (getarrival(results(i).useralign.data,...
                {results(i).phase(1) results(i).phase}) ...
                - getheader(results(i).useralign.data,'o'));
            seizmoverbose(verbose);
            x{i,1}(results(i).outliers.bad)=nan;
            x{i,2}=results(i).useralign.solution.arrerr;
            x{i,2}(results(i).outliers.bad)=nan;
            snr=results(i).usersnr.snr;
            snr=snr(snr>=results(i).usersnr.snrcut);
            snr(results(i).userwinnow.cut)=[];
            snr=snr(results(i).finalcut);
            x{i,3}=log10(snr);
            x{i,3}(results(i).outliers.bad)=nan;
            x{i,4}=log10(results(i).useralign.solution.amp);
            x{i,4}(results(i).outliers.bad)=nan;
    end
end

% yaxis names and ordering
[s,i]=unique(cat(1,kname{:}));
if(~mevflag)
    d=cat(1,d{:}); d=d(i);
    [d,i]=sort(d,'descend'); s=s(i);
end

% plot all or single observable
if(strcmpi(observable,'all'))
    % make observable matrices
    o=nan(size(s,1),numel(bad),4);
    for i=1:numel(results)
        [tf,loc]=ismember(kname{i},s);
        for j=1:4; o(loc,oi(i),j)=x{i,j}; end
    end
    
    % plot observables
    h(1)=imagesc(o(:,:,1),'parent',ax(1));
    h(2)=imagesc(o(:,:,2),'parent',ax(2));
    h(3)=imagesc(o(:,:,3),'parent',ax(3));
    h(4)=imagesc(o(:,:,4),'parent',ax(4));
    set(ax(1:4),'color','none',...
        'ytick',(1:numel(s)),...
        'yticklabel',s,...
        'xtick',oi);
    colorbar('peer',ax(1));
    colorbar('peer',ax(2));
    colorbar('peer',ax(3));
    colorbar('peer',ax(4));
    set(h(1),'alphadata',~isnan(o(:,:,1)));
    set(h(2),'alphadata',~isnan(o(:,:,2)));
    set(h(3),'alphadata',~isnan(o(:,:,3)));
    set(h(4),'alphadata',~isnan(o(:,:,4)));
    title(ax(1),'arr');
    title(ax(2),'err');
    title(ax(3),'snr');
    title(ax(4),'amp');
    xlabel(ax(1),'index #');
    xlabel(ax(2),'index #');
    xlabel(ax(3),'index #');
    xlabel(ax(4),'index #');
    
    % plot profile info
    pf=slowdecayprofiles(results,[],[],false);
    plot_cmb_pdf(pf,'cslow',[],[],ax(5));
    plot_cmb_pdf(pf,'cdecay',[],[],ax(6));
else
    % make observable matrix
    o=nan(size(s,1),numel(bad));
    for i=1:numel(results)
        [tf,loc]=ismember(kname{i},s);
        o(loc,oi(i))=x{i};
    end
    
    % plot observable
    h=imagesc(o,'parent',ax);
    set(ax,'color','none',...
        'ytick',(1:numel(s)),...
        'yticklabel',s,...
        'xtick',oi);
    colorbar('peer',ax);
    set(h,'alphadata',~isnan(o));
    title(ax,observable);
    xlabel(ax,'index #');
end

end
