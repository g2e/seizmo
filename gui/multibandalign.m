function [info]=multibandalign(data,varargin)
%MULTIBANDALIGN    Aligns signal at multiple frequency bands
%
%    Usage:    results=multibandalign(data)
%              results=multibandalign(data,'option',value,...)
%
%    Description:
%     RESULTS=MULTIBANDALIGN(DATA) presents an interface to align the data
%     given in SEIZMO struct DATA over a series of 25 frequency bands from
%     50s to 5s.  Details on the default frequency bands can be found in
%     the Notes section below.  To alter the frequency bands see the next
%     usage form.  The pre-alignment processing includes a SNR-based data
%     removal and distance based data removal.  Alignment uses USERALIGN.
%     All pertinent info is saved to the RESULTS struct, which includes the
%     USERSNR, USERWINNOW, & USERALIGN output.  All figures are saved using
%     a name called 'multibandalign_date&time_figname.fig.
%
%     RESULTS=MULTIBANDALIGN(DATA,BANK,RUNNAME,'OPTION',VALUE,...) allows
%     altering certain parameters.  Those pertinent to MULTIBANDALIGN:
%       runname         - Name used for saving figures
%       initwin         - initial window limits (default is [-150 200])
%       noisewin        - noise window (default is [-125 25])
%       signalwin       - signal window (default is [35 135])
%       snrmethod       - SNR method (default is 'peak2rms')
%       snrcutoff       - Default cutoff (3) for removing low SNR records
%       filterbank      - filter bank (use FILTER_BANK)
%       filtertype      - filter type (default is 'bp')
%       filterstyle     - filter style (default is 'butter')
%       filterorder     - filter order (default is 4)
%       filterpasses    - number of passes (default is 1)
%       auto            - fully automated mode (false)
%       zxing           - zero-crossing windowing (true)
%       zxoffset        - auto mode: # zero xings from guess to win start
%       zxwidth         - auto mode: window is # zero crossings wide
%       winnow          - prompt user to winnow (false)
%       winnowyfield    - header field to winnow by (default is 'gcarc')
%       winnowlimits    - default winnow limits (default is [])
%       figdir          - directory to save figures to (default is '.')
%       lags            - number of periods to look from 0 for xc peaks
%     The remaining options are passed to USERALIGN.  See that function for
%     more details.
%
%    Notes:
%     - The default filter bank is created using:
%       filter_bank([1/50 1/5],'variable',0.2,0.1)
%
%    Examples:
%     % Use a filter bank with filters of 10mHz width and stepped at 5mHz:
%     linbank=filter_bank([1/50 1/5],'constant',.01,.005);
%     results=multibandalign(data,'bank',linbank,'runname','examplerun');
%
%     % Save figures to a separate directory and skip winnowing:
%     results=multibandalign(data,'figdir','figs','winnow',false,...
%                                 'runname','examplerun');
%
%    See also: USERALIGN, USERSNR, FILTER_BANK, IIRFILTER, CMB_2ND_PASS

%     Version History:
%        Mar. 25, 2010 - initial version
%        Sep. 21, 2010 - working version
%        Sep. 30, 2010 - added amplitude measurements, adjust output
%        Oct.  1, 2010 - added runname input
%        Nov.  5, 2010 - userwinnow added
%        Nov.  9, 2010 - prints filter corners for each band, save
%                        adjustment to initial window for subsequent bands,
%                        workaround for ill-conditioning of pre-alignment,
%                        better capture of erroring out and useralign,
%                        allow redo of a an entire filter band sequence,
%                        save snr windows for subsequent bands
%        Nov. 13, 2010 - update for userwindow arg change
%        Nov. 21, 2010 - force userwinnow single trace normstyle
%        Jan. 14, 2011 - added a commented bit about amplitudes, added fast
%                        mode, improved some menus, info is now col vector
%        Jan. 18, 2011 - multibandalign options added, winnowing memory,
%                        improved snr codes, bank & runname are optional,
%                        zero-crossing feature added
%        Jan. 23, 2011 - save window positions on iter 2+ of first band,
%                        phase input allows for a bit more
%        Jan. 26, 2011 - use 2 digit band number
%        Jan. 29, 2011 - alter filter bank from 1/80-1/8 to 1/50-1/5, fully
%                        automated mode, more params added, fix stack bug,
%                        fix phase bug, fixed alignment reuse bug
%        Jan. 30, 2011 - zxing option allows turning on/off zx codes,
%                        advance to next band if too few high snr or in
%                        winnow range to avoid infinite loop in auto mode
%        Jan. 31, 2011 - some doc fixes, figdir option
%        Feb.  1, 2011 - added autozx checks
%        Feb.  6, 2011 - make userwinnow optional, adjust zxing code to be
%                        more useful for the default filter bank
%        Feb. 26, 2011 - fix zeroxing bug
%        Mar.  1, 2011 - fix winnow not using its state, better examples
%        Mar.  7, 2011 - winnow is off by default (annoying)
%        Mar. 17, 2011 - only uses useralign_auto now (no poweruser mode),
%                        usermoveout added to main menu
%        Apr.  6, 2011 - no error if aligned plot closed, fix normstyle for
%                        editing initial window or moveout
%        Apr.  7, 2011 - msg for pre-aligning step (can be time-consuming)
%        Apr. 17, 2011 - lags option for correlate (input should be in # of
%                        periods rather than lags), initwin bugfix
%        Apr. 22, 2011 - added post-align deletion, delete useralign output
%                        if user did not like it
%        Apr. 23, 2011 - window presets are more appropriate to default
%                        filter set, fix bug where window settings lost if
%                        too few high-snr waveforms
%        May  19, 2011 - fixed error with .finalcut not always there
%        Mar.  1, 2012 - old plot_taupcurve is now plot_taupcurve_dt
%        Mar.  5, 2012 - allow no written output by setting figdir=false
%        Mar. 15, 2012 - fix for pick functions
%        Apr.  3, 2012 - use seizmocheck
%        Jan. 28, 2013 - deleted permutation of xc outputs, fixed phase
%                        parsing for findpicks, fixed figdir handling
%        Jan. 30, 2013 - update for new correlate (no spacing option)
%        Feb. 14, 2013 - bugfix for xc output squeeze
%        Feb. 27, 2013 - bugfix: winnow cut uninitialized, adjust default
%                        windowing parameters & zx parameters so full auto
%                        works much better
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 27, 2013 at 12:25 GMT

% todo:

% check nargin
if(nargin<1)
    error('seizmo:multibandalign:notEnoughInputs',...
        'Not enough input arguments.');
elseif(nargin>1 && ~mod(nargin,2))
    error('seizmo:multbandalign:optionMustBePaired',...
        'Options must be paired with a value!');
end

% check data structure
error(seizmocheck(data,'dep'));

% number of records
nrecs=numel(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt multi-band align
try
    % parse inputs
    %[bank,runname,snrcut,filt,phase,varargin]=parse_mba_params(varargin{:});
    [p,varargin]=parse_mba_params(varargin{:});
    
    % preallocate output
    info([])=struct('useralign',[],'filter',[],'initwin',[],...
        'usermoveout',[],'usersnr',[],'userwinnow',[]);
    
    % require o field is set to the same value (within +/- 2 millisec)
    o=getheader(data,'o');
    outc=cell2mat(getheader(data,'o utc'));
    if(any(abs(timediff(outc(1,:),outc))>0.002))
        error('seizmo:multibandalign:oFieldVaries',...
            'O field must correspond to one UTC time for all records!');
    end
    
    % build lag table (for pre-aligning each band using previous lags)
    goodarr=-o;
    goodidx=(1:nrecs)';
    lags=o(:,ones(1,nrecs))'-o(:,ones(1,nrecs));
    idx=ones(nrecs);
    
    % loop over each filter
    i=1;         % starting with first band
    last=[];     % for memory (last band aligned)
    while(i<=size(p.bank,1))
        % display filter info
        istr=num2str(i,'%02d');
        disp(['BAND ' istr ' FILTER CORNERS: ' ...
            num2str(1/p.bank(i,3)) 's  to  ' num2str(1/p.bank(i,2)) 's']);
        
        % filter the records
        data0=iirfilter(data,p.filt.type,p.filt.style,...
            'c',p.bank(i,2:3),'o',p.filt.order,'p',p.filt.passes);
        p.filt.corners=p.bank(i,2:3);
        info(i,1).filter=p.filt;
        
        % get new relative arrivals & prealign
        % - lags are absolute
        % - use 10^idx as weight (higher weight to more recent freqs)
        %   - limit to >10^-5 for stability (-5 may be changed)
        % - use latest "good" arrivals to get absolute times
        disp('Pre-aligning using info from previous alignment(s)');
        tt=ttalign(lags,10.^(max(idx-i,-5)),goodarr,1,goodidx);
        data0=timeshift(data0,-o-tt); % shift to origin then new times
        info(i,1).tt_start=tt;
        
        % initial window and assessment
        info(i).initwin=p.initwin;
        if(p.auto) % skipping main menu
            % no skipping filter bands
            skip=0;
            
            % preview plot
            ax=preview_plot(data0,info(i),p);
            if(ishandle(ax))
                if(p.figout)
                    saveas(get(ax,'parent'),fullfile(p.figdir,...
                        [datestr(now,30) '_' p.runname '_band_' istr ...
                        '_preview.fig']));
                end
                close(get(ax,'parent'));
            end
            data0=cut(data0,info(i).initwin(1),info(i).initwin(2));
        else % main menu
            [data0,info(i),skip,p,ax]=get_initial_window(data0,...
                info(i),p,i);
            if(ishandle(ax))
                if(p.figout)
                    saveas(get(ax,'parent'),fullfile(p.figdir,...
                        [datestr(now,30) '_' p.runname '_band_' istr ...
                        '_preview.fig']));
                end
                close(get(ax,'parent'));
            end
        end
        
        % handle skipping
        info(i).finalcut=[];
        if(skip); i=i+skip; continue; end
        
        % auto-window or user specified (usersnr interface)
        if(p.auto)
            % default settings
            if(i==1)
                info(i).usersnr.noisewin=p.noisewin;
                info(i).usersnr.signalwin=p.signalwin;
                info(i).usersnr.method=p.snrmethod;
            else
                info(i).usersnr.noisewin=info(last).usersnr.noisewin;
                info(i).usersnr.signalwin=info(last).usersnr.signalwin;
                info(i).usersnr.method=p.snrmethod;
            end
            
            % zx-based window
            if(p.zxing);
                zx(1)=zeroxing321(data0,info(i),1);
                zx(2)=zeroxing(data0,info(i),2);
                info(i).usersnr.signalwin=zx;
                %info(i).usersnr.signalwin=autozx(data0,info(i),p);
            end
            
            % plot windows
            info(i).usersnr.plottype=@plot0;
            ax=windows_plot(data0,info(i));
            if(ishandle(ax))
                if(p.figout)
                    saveas(get(ax,'parent'),fullfile(p.figdir,...
                        [datestr(now,30) '_' p.runname '_band_' istr ...
                        '_windows.fig']));
                end
                close(get(ax,'parent'));
            end
            
            % get snr
            snr=quicksnr(data0,info(i).usersnr.noisewin,...
                info(i).usersnr.signalwin,info(i).usersnr.method);
            info(i).usersnr.snr=snr;
            info(i).usersnr.snrcut=p.snrcut;
        else
            % user specified snr
            if(isempty(last))
                info(i).usersnr.noisewin=p.noisewin;
                info(i).usersnr.signalwin=p.signalwin;
                info(i).usersnr.method=p.snrmethod;
            elseif(last==i)
                % just update signal window if we used useralign
                if(~isempty(info(last).useralign) && ...
                        ~isempty(info(last).useralign.userwindow.limits))
                    info(i).usersnr.signalwin=...
                        info(last).useralign.userwindow.limits;
                end
            else % look for zero-crossing
                info(i).usersnr.noisewin=info(last).usersnr.noisewin;
                if(~isempty(info(last).useralign) && ...
                        ~isempty(info(last).useralign.userwindow.limits))
                    info(i).usersnr.signalwin=...
                        info(last).useralign.userwindow.limits;
                else
                    info(i).usersnr.signalwin=info(last).usersnr.signalwin;
                end
                info(i).usersnr.method=info(last).usersnr.method;
                
                % use nearby zero crossing for new window limits
                if(p.zxing)
                    zx(1)=zeroxing321(data0,info(i),1);
                    zx(2)=zeroxing(data0,info(i),2);
                    info(i).usersnr.signalwin=zx;
                end
            end
            [snr,info(i).usersnr,ax]=usersnr(data0,...
                info(i).usersnr.noisewin,info(i).usersnr.signalwin,...
                info(i).usersnr.method,'normstyle','single');
            info(i).usersnr.snr=snr;
            info(i).usersnr.snrcut=p.snrcut;
            if(ishandle(ax))
                if(p.figout)
                    saveas(get(ax,'parent'),fullfile(p.figdir,...
                        [datestr(now,30) '_' p.runname '_band_' istr ...
                        '_usersnr.fig']));
                end
                close(get(ax,'parent'));
            end
        end
        
        % trimming dataset to those above snr cutoff
        snridx=find(snr>=p.snrcut);
        data0=data0(snridx);
        snr=snr(snridx);
        nn=numel(data0);
        
        % skip to next band if less than 3
        if(nn<3)
            last=i;
            i=i+1;
            disp(['BAND ' istr ': Too Few High SNR data!']);
            continue;
        end
        
        % distance winnow only if manual
        if(~p.userwinnow || p.auto)
            % fake userwinnow
            info(i).userwinnow.yfield=p.winnowyfield;
            if(isempty(p.winnowlimits))
                info(i).userwinnow.limits=[];
                info(i).userwinnow.cut=[];
            else % winnow!
                % get winnow
                info(i).userwinnow.limits=p.winnowlimits;
                yfield=getheader(data0,p.winnowyfield);
                if(p.winnowlimits(1)>p.winnowlimits(2))
                    info(i).userwinnow.cut=find(...
                        yfield>=info(i).userwinnow.limits(1) ...
                        | yfield<=info(i).userwinnow.limits(2));
                else % normal case
                    info(i).userwinnow.cut=find(...
                        yfield<info(i).userwinnow.limits(1) ...
                        | yfield>info(i).userwinnow.limits(2));
                end
                
                % plot it
                ax=recordsection(data0,'yfield',p.winnowyfield,...
                    'normstyle','single');
                span=xlim(ax);
                hold(ax,'on');
                plot(ax,span,[p.winnowlimits(1) p.winnowlimits(1)],'g',...
                    'linewidth',4);
                plot(ax,span,[p.winnowlimits(2) p.winnowlimits(2)],'r',...
                    'linewidth',4);
                hold(ax,'off');
                if(ishandle(ax))
                    if(p.figout)
                        saveas(get(ax,'parent'),fullfile(p.figdir,...
                            [datestr(now,30) '_' p.runname '_band_' ...
                            istr '_winnow.fig']));
                    end
                    close(get(ax,'parent'));
                end
                
                % winnow
                data0(info(i).userwinnow.cut)=[];
                snr(info(i).userwinnow.cut)=[];
                snridx(info(i).userwinnow.cut)=[];
                nn=numel(data0);
            end
        else
            if(isempty(last));
                info(i).userwinnow.limits=p.winnowlimits;
            else
                try
                    if(isempty(info(last).userwinnow.limits))
                        info(i).userwinnow.limits=p.winnowlimits;
                    else
                        info(i).userwinnow.limits=...
                            info(last).userwinnow.limits;
                    end
                catch
                    info(i).userwinnow.limits=p.winnowlimits;
                end
            end
            [data0,info(i).userwinnow,ax]=userwinnow(data0,...
                info(i).userwinnow.limits,...
                'normstyle','single',...
                'yfield',p.winnowyfield);
            snr(info(i).userwinnow.cut)=[];
            snridx(info(i).userwinnow.cut)=[];
            if(ishandle(ax))
                if(p.figout)
                    saveas(get(ax,'parent'),fullfile(p.figdir,...
                        [datestr(now,30) '_' p.runname '_band_' istr ...
                        '_userwinnow.fig']));
                end
                close(get(ax,'parent'));
            end
            nn=numel(data0);
        end
        
        % skip to next band if less than 3
        if(nn<3)
            last=i;
            i=i+1;
            disp(['BAND ' istr ': Too few data in ' ...
                p.winnowyfield ' range!']);
            continue;
        end
        
        % quiet align
        [info(i).useralign,info(i).useralign.xc,...
            info(i).useralign.data]=useralign_quiet(data0,'absxc',false,...
            'estarr',zeros(nn,1),'snr',snr,'lags',p.lags/p.bank(i,3),...
            'window',info(i).usersnr.signalwin,varargin{:});
        
        % Get Amplitude Info
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % THIS IS THE MOST VISUALLY APPEALING
        % BUT ERRORS TEND TO BE OVERESTIMATED BY ~10%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [info(i).useralign.data,scale]=normalize(...
            info(i).useralign.data);
        info(i).useralign.solution.amp=scale;
        info(i).useralign.solution.amperr=scale./snr;
        
        % only ask if satisfied if not auto
        info(i).finalcut=true(numel(info(i).useralign.data),1);
        ax=plot2(info(i).useralign.data);
        if(p.auto)
            if(ishandle(ax))
                if(p.figout)
                    saveas(get(ax,'parent'),fullfile(p.figdir,...
                        [datestr(now,30) '_' p.runname '_band_' istr ...
                        '_aligned.fig']));
                end
                % this was missing before (probably for testing)
                close(get(ax,'parent'));
            end
        else
            satisfied=false;
            deleted=false(nn,1);
            while(~satisfied)
                choice=user_satisfied;
                if(ishandle(ax))
                    if(p.figout)
                        saveas(get(ax,'parent'),fullfile(p.figdir,...
                            [datestr(now,30) '_' p.runname '_band_' ...
                            istr '_aligned.fig']));
                    end
                    close(get(ax,'parent'));
                end
                switch choice
                    case 'yes'
                        satisfied=true;
                        break; % skip replotting
                    case 'no'
                        info(i).useralign=[]; % user didn't like it
                        last=i;
                        break; % gets caught after while loop
                    case 'delete'
                        [deleted,deleted,ax]=selectrecords(...
                            info(i).useralign.data,'delete','p2',deleted);
                        if(ishandle(ax)); close(get(ax,'parent')); end
                end
                ax=plot2(info(i).useralign.data(~deleted));
            end
            % catch no from above (redo band)
            if(~satisfied); continue; end
            
            % save/delete info
            info(i).finalcut=~deleted;
            info(i).useralign.solution.arr(deleted)=[];
            info(i).useralign.solution.arrerr(deleted)=[];
            info(i).useralign.solution.amp(deleted)=[];
            info(i).useralign.solution.amperr(deleted)=[];
            info(i).useralign.solution.pol(deleted)=[];
            info(i).useralign.solution.zmean(deleted)=[];
            info(i).useralign.solution.zstd(deleted)=[];
            info(i).useralign.data(deleted)=[];
            info(i).useralign.usermoveout.adjust(deleted)=[];
            bigdeleted=ndsquareform(deleted(:,ones(nn,1)) ...
                | deleted(:,ones(nn,1))');
            info(i).useralign.xc.cg(bigdeleted,:,:)=[];
            info(i).useralign.xc.lg(bigdeleted,:,:)=[];
            info(i).useralign.xc.pg(bigdeleted,:,:)=[];
            info(i).useralign.xc.zg(bigdeleted,:,:)=[];
            info(i).useralign.xc.wg(bigdeleted,:,:)=[];
            clear bigdeleted;
            
            % update variables
            nn=sum(~deleted);
            snridx=snridx(~deleted);
        end
        
        % use alignment for next band
        if(i<size(p.bank,1))
            goodarr=info(i).useralign.solution.arr;
            goodidx=snridx;
            lags(goodidx,goodidx)=goodarr(:,ones(1,nn))...
                -goodarr(:,ones(1,nn))';
            idx(goodidx,goodidx)=i+1;
        end
        
        % we aligned this band so lets use its info for later bands
        last=i;
        i=i+1;
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end


function [choice]=user_satisfied()
happy_user=false;
answer={'yes' 'no' 'delete'};
while(~happy_user)
    choice=menu('Satisfied w/ alignment?',...
        'YES','NO','DELETE SOME RECORDS');
    if(choice); choice=answer{choice}; happy_user=true; end
end
end


function [data,info,skip,p,ax]=get_initial_window(data,info,p,bandidx)

% default initial window
if(isempty(info.initwin))
    info.initwin=p.initwin;
end

% shift if phase given
if(~isempty(p.phase))
    [t,n]=findpicks(data,[p.phase(1) ',' p.phase],1);
    data=timeshift(data,-t,strcat('it',num2str(n)));
end

% additional moveout/time shifts
moveout=0;
shift=zeros(size(data));

% ask the user if they wish to continue
skip=0; happy_user=false; ax=-1;
while(~happy_user)
    % plot records (so user can decide to skip or not)
    if(~ishandle(ax))
        ax=preview_plot(data,info,p);
    end
    
    choice=menu('What to do with this filter band?',...
        'Align The Waveforms!',...
        ['Adjust SNR Cutoff (' num2str(p.snrcut) ')'],...
        'Adjust Initial Window',...
        'Adjust Moveout',...
        'Skip This Filter Band',...
        'Skip Remaining Filter Bands',...
        'Jump To Filter Band ...');
    switch choice
        case 1 % process
            happy_user=true;
        case 2 % snrcut
            tmp=inputdlg(...
                ['SNR CutOff? [' num2str(p.snrcut) ']:'],...
                'Enter SNR CutOff',1,{num2str(p.snrcut)});
            if(~isempty(tmp))
                try
                    tmp=str2double(tmp{:});
                    if(~isnan(tmp) && ~isinf(tmp))
                        p.snrcut=tmp;
                    end
                catch
                    % do not change snrcut
                end
            end
        case 3 % adjust window
            % get new initial window
            if(ishandle(ax)); close(get(ax,'parent')); end
            [win,win,ax]=userwindow(data,info.initwin,[],[],...
                'normstyle','single');
            if(ishandle(ax)); close(get(ax,'parent')); end
            
            % default initial window if none
            if(~isempty(win.limits))
                info.initwin=win.limits;
                p.initwin=win.limits;
            end
        case 4 % moveout
            % get new initial window
            if(ishandle(ax)); close(get(ax,'parent')); end
            [mvo,mvo,ax]=usermoveout(...
                cut(data,info.initwin(1),info.initwin(2)),...
                'normstyle','single');
            if(ishandle(ax)); close(get(ax,'parent')); end
            data=timeshift(data,mvo.shift);
            
            % get total moveout and time shifts
            shift=shift+mvo.shift;
            moveout=moveout+mvo.moveout;
        case 5 % skip this band
            happy_user=true;
            skip=1;
        case 6 % skip all
            happy_user=true;
            skip=inf;
        case 7 % skip to ...
            nb=size(p.bank,1);
            n=cellstr(num2str((1:nb)','%02d'));
            p1=cellstr(num2str(1./p.bank(:,3),'%5.1f'));
            p2=cellstr(num2str(1./p.bank(:,2),'%5.1f'));
            cur={''}; cur=cur(ones(nb,1),1); cur{bandidx}='(current)';
            skip=listdlg(...
                'liststring',...
                strcat({'BAND '},n,{':  '},p1,'-',p2,{'s '},cur),...
                'selectionmode','single',...
                'promptstring','Jump to bandpass:',...
                'initialvalue',bandidx,...
                'listsize',[240 400]);
            if(isempty(skip) || skip==bandidx)
                skip=0;
            else
                happy_user=true;
                skip=skip-bandidx;
            end
    end
end

% implement window
if(choice==1); data=cut(data,info.initwin(1),info.initwin(2)); end
info.usermoveout.moveout=moveout;
info.usermoveout.shift=shift;

end


function [p,varargin]=parse_mba_params(varargin)
%PARSE_MBA_PARAMS  Parses out multibandalign parameters

% basic checks on optional inputs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:multibandalign:badInput',...
        'All OPTIONs must be specified with a string!');
end

% defaults
p.initwin=[-150 200];
p.noisewin=[-150 0];
p.signalwin=[50 170];
p.bank=filter_bank([1/50 1/5],'variable',0.2,0.1);
p.runname=['multibandalign_' datestr(now,30)];
p.snrcut=5;
p.snrmethod='peak2rms';
p.filt.type='bandpass';
p.filt.style='butter';
p.filt.order=4;
p.filt.passes=1;
p.phase=[];   % pre-align phase field (only for core-diff currently)
p.auto=false; % auto multibandalign (uses zx window parameters if zxing)
p.zxing=true; % zero-crossing windowing flag (off uses signalwin or user)
p.zxwidth=6;  % number of zero crossings increments in window if auto
p.zxoffset=1; % number of zero crossings to window start from initial guess
p.userwinnow=false;
p.winnowyfield='gcarc';
p.winnowlimits=[];
p.figdir='.';
p.figout=true;
p.lags=[];

% valid filter types & styles (grabbed from iirdesign)
validtypes={'low' 'lo' 'l' 'lp' 'high' 'hi' 'h' 'hp' 'bandpass' 'pass' ...
    'bp' 'bandstop' 'stop' 'bs' 'notch' 'n' 'bandreject' 'reject' 'br'};
validstyles={'butter' 'butt' 'butterworth' 'bu' 'b' 'cheby1' 'cheb1' ...
    'chebychev1' 'c1' 'cheby2' 'cheb2' 'chebychev2' 'c2' 'ellip' 'el' ...
    'elliptical' 'e'};

% find multibandalign options
keep=true(nargin,1);
for i=1:2:nargin
    if(isempty(varargin{i+1})); continue; end
    switch lower(varargin{i})
        case {'fb' 'filterbank' 'bank'}
            if(size(varargin{i+1},2)~=3 || any(varargin{i+1}(:)<=0 ...
                    | isnan(varargin{i+1}(:)) ...
                    | isinf(varargin{i+1}(:))) ...
                    || any(varargin{i+1}(:,1)<=varargin{i+1}(:,2) ...
                    | varargin{i+1}(:,3)<=varargin{i+1}(:,1)))
                error('seizmo:multibandalign:badInput',...
                    'BANK must be in the format from FILTER_BANK!');
            end
            p.bank=varargin{i+1};
            keep(i:i+1)=false;
        case {'rn' 'run' 'name' 'runname'}
            if(~isstring(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'RUNNAME must be a string!');
            end
            p.runname=varargin{i+1};
            keep(i:i+1)=false;
        case {'sc' 'snrcutoff' 'snrcut'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}<0)
                error('seizmo:multibandalign:badInput',...
                    'SNRCUTOFF must be a real-valued scalar!');
            end
            p.snrcut=varargin{i+1};
            keep(i:i+1)=false;
        case {'ft' 'filtertype'}
            if(~isstring(varargin{i+1}) ...
                    || ~any(strcmpi(varargin{i+1},validtypes)))
                error('seizmo:multibandalign:badInput',...
                    'FILTERTYPE must be a filter type in IIRDESIGN!');
            end
            p.filt.type=varargin{i+1};
            keep(i:i+1)=false;
        case {'fs' 'filterstyle'}
            if(~isstring(varargin{i+1}) ...
                    || ~any(strcmpi(varargin{i+1},validstyles)))
                error('seizmo:multibandalign:badInput',...
                    'FILTERSTYLE must be a filter style in IIRDESIGN!');
            end
            p.filt.style=varargin{i+1};
            keep(i:i+1)=false;
        case {'fo' 'filterorder'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}~=fix(varargin{i+1}) ...
                    || varargin{i+1}<=0)
                error('seizmo:multibandalign:badInput',...
                    'FILTERORDER must be a scalar integer >0!');
            end
            p.filt.order=varargin{i+1};
            keep(i:i+1)=false;
        case {'fp' 'filterpasses' 'filterpass'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}~=fix(varargin{i+1}) ...
                    || ~any(abs(varargin{i+1})==[1 2]))
                error('seizmo:multibandalign:badInput',...
                    'FILTERPASSES must be 1 or 2!');
            end
            p.filt.passes=varargin{i+1};
            keep(i:i+1)=false;
        case {'phase'}
            if(~isstring(varargin{i+1}) ...
                    || ~any(strcmpi(varargin{i+1},{'Pdiff' 'Sdiff'})))
                error('seizmo:multibandalign:badInput',...
                    'PHASE must be a ''Pdiff'' or ''Sdiff''!');
            end
            p.phase=varargin{i+1};
            keep(i:i+1)=false;
        case 'auto'
            if(~isscalar(varargin{i+1}) || (~islogical(varargin{i+1}) ...
                    && ~isreal(varargin{i+1})))
                error('seizmo:multibandalign:badInput',...
                    'AUTO must be TRUE or FALSE!');
            end
            p.auto=varargin{i+1};
            keep(i:i+1)=false;
        case 'zxing'
            if(~isscalar(varargin{i+1}) || (~islogical(varargin{i+1}) ...
                    && ~isreal(varargin{i+1})))
                error('seizmo:multibandalign:badInput',...
                    'ZXING must be TRUE or FALSE!');
            end
            p.zxing=varargin{i+1};
            keep(i:i+1)=false;
        case 'zxwidth'
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}<=0 ...
                    || varargin{i+1}~=fix(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'ZXWIDTH must be a positive integer!');
            end
            p.zxwidth=varargin{i+1};
            keep(i:i+1)=false;
        case 'zxoffset'
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}~=fix(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'ZXOFFSET must be an integer!');
            end
            p.zxoffset=varargin{i+1};
            keep(i:i+1)=false;
        case 'initwin'
            if(~isequal(size(varargin{i+1}),[1 2]) ...
                    || ~isreal(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'INITWIN must be a 1x2 real-valued array!');
            end
            p.initwin=varargin{i+1};
            keep(i:i+1)=false;
        case 'signalwin'
            if(~isequal(size(varargin{i+1}),[1 2]) ...
                    || ~isreal(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'SIGNALWIN must be a 1x2 real-valued array!');
            end
            p.signalwin=varargin{i+1};
            keep(i:i+1)=false;
        case 'noisewin'
            if(~isequal(size(varargin{i+1}),[1 2]) ...
                    || ~isreal(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'NOISEWIN must be a 1x2 real-valued array!');
            end
            p.noisewin=varargin{i+1};
            keep(i:i+1)=false;
        case 'snrmethod'
            if(~isstring(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'SNRMETHOD must be a string!');
            end
            p.snrmethod=varargin{i+1};
            keep(i:i+1)=false;
        case {'winnow' 'userwinnow'}
            if(~isscalar(varargin{i+1}) || (~islogical(varargin{i+1}) ...
                    && ~isreal(varargin{i+1})))
                error('seizmo:multibandalign:badInput',...
                    'WINNOW must be TRUE or FALSE!');
            end
            p.userwinnow=varargin{i+1};
            keep(i:i+1)=false;
        case 'winnowlimits'
            if(~isequal(size(varargin{i+1}),[1 2]) ...
                    || ~isreal(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'WINNOWLIMITS must be a 1x2 real-valued array!');
            end
            p.winnowlimits=varargin{i+1};
            keep(i:i+1)=false;
        case 'winnowyfield'
            if(~isstring(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'WINNOWYFIELD must be a string!');
            end
            p.winnowyfield=varargin{i+1};
            keep(i:i+1)=false;
        case {'fdir' 'figdir'}
            if(~isstring(varargin{i+1}) && ~(islogical(varargin{i+1}) ...
                    && isscalar(varargin{i+1})))
                error('seizmo:multibandalign:badInput',...
                    'FIGDIR must be a string or TRUE/FALSE!');
            end
            p.figdir=varargin{i+1};
            keep(i:i+1)=false;
        case {'lag' 'lags'}
            if(numel(varargin{i+1})>2 || ~isnumeric(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'LAGS must be given as [MINLAG MAXLAG]!');
            end
            p.lags=varargin{i+1};
            keep(i:i+1)=false;
    end
end
varargin=varargin(keep);

% force fast to true if auto is true
if(p.auto); p.fast=true; end

% fix TRUE directory input
if(islogical(p.figdir))
    if(~p.figdir); p.figout=false; end
    p.figdir='.';
end

% create directory if it does not exist
% check that this does not fail
if(p.figout)
    [ok,msg,msgid]=mkdir(p.figdir);
    if(~ok)
        warning(msgid,msg);
        error('seizmo:multibandalign:pathBad',...
            'Cannot create directory: %s',p.figdir);
    end
end

end


function [zx]=zeroxing(data,info,idx)
% get sample spacing
delta=getheader(data,'delta');
delta=unique(delta);
if(~isscalar(delta))
   error('seizmo:multibandalign:badInput',...
       'DATA records must have the same sample rate!');
end

% get stack
dstack=stack(data,1/delta,[],info.initwin(1),info.initwin(2));

% get zero crossings
zx=zerocrossings(dstack);
zx=sort(zx{1}); % just in case

% find closest within range
sw=info.usersnr.signalwin;
neg=find(zx<sw(idx) & zx>=(sw(1)+sw(2))/2,1,'last');

% choose one
if(isempty(neg))
    % use previous
    zx=sw(idx);
else
    zx=zx(neg);
end

end


function [zx]=zeroxing321(data,info,idx)
% get sample spacing
delta=getheader(data,'delta');
delta=unique(delta);
if(~isscalar(delta))
   error('seizmo:multibandalign:badInput',...
       'DATA records must have the same sample rate!');
end

% get stack
dstack=stack(data,1/delta,[],info.initwin(1),info.initwin(2));

% get zero crossings
zx=zerocrossings(dstack);
zx=sort(zx{1}); % just in case

% find closest within range
sw=info.usersnr.signalwin;
neg=find(zx<sw(idx) & zx>=sw(idx)-(sw(2)-sw(1))/2,1,'last');
pos=find(zx>=sw(idx) & zx<=(sw(1)+sw(2))/2,1,'first');

% choose one
if(isempty(neg) && isempty(pos))
    % use previous
    zx=sw(idx);
elseif(isempty(neg))
    zx=zx(pos);
elseif(isempty(pos))
    zx=zx(neg);
else
    % nearest zero crossing
    if((sw(idx)-zx(neg))<=(zx(pos)-sw(idx)))
        zx=zx(neg);
    else
        zx=zx(pos);
    end
end

end


function [ax]=preview_plot(data,info,p)
% record section
ax=recordsection(cut(data,info.initwin(1),info.initwin(2)),...
    'normstyle','single','xlim',info.initwin,'title',...
    ['FILTER CORNERS: ' num2str(1/info.filter.corners(2)) ...
    's  to  ' num2str(1/info.filter.corners(1)) 's']);

% this all assumes the phase is shifted to 0
if(~isempty(p.phase))
    evdp=getheader(data(1),'evdp')/1000;
    tt=taupcurve('dep',evdp);
    idx=find(strcmp(p.phase,{tt.phase}));
    intrcpt=tt(idx).time(1)...
        -tt(idx).distance(1)*tt(idx).rayparameter(1);
    tt=taupcurve('dep',evdp,...
        'reddeg',1/tt(idx).rayparameter(1),'ph','ttall');
    hold(ax,'on');
    h=plot_taupcurve_dt(tt,-intrcpt,true,'parent',ax,'linewidth',5);
    movekids(h,'back');
    hold(ax,'off');
end
end


function [ax]=windows_plot(data,info)
% make title
ptitle={['SNR Estimation Method: ' upper(info.usersnr.method)]
    'Yellow Dashed Line --  Noise Window Start'
    '  Blue Dashed Line --  Noise Window End  '
    '  Green Solid Line -- Signal Window Start'
    '    Red Solid Line -- Signal Window End  '};

% plot records
ax=info.usersnr.plottype(data,'title',ptitle);

% add window limit markers
% - yellow/blue dashed == noise window
% - red/green == signal window
span=ylim(ax);
hold(ax,'on');
plot(ax,[info.usersnr.noisewin(1) info.usersnr.noisewin(1)],...
    span,'--y','linewidth',4);
plot(ax,[info.usersnr.noisewin(2) info.usersnr.noisewin(2)],...
    span,'--b','linewidth',4);
plot(ax,[info.usersnr.signalwin(1) info.usersnr.signalwin(1)],...
    span,'g','linewidth',4);
plot(ax,[info.usersnr.signalwin(2) info.usersnr.signalwin(2)],...
    span,'r','linewidth',4);
hold(ax,'off');
end

