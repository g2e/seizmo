function [info]=multibandalign(data,varargin)
%MULTIBANDALIGN    Aligns signal at multiple frequency bands
%
%    Usage:    results=multibandalign(data)
%              results=multibandalign(data,'option',value,...)
%
%    Description:
%     RESULTS=MULTIBANDALIGN(DATA) presents an interface to align the data
%     given in SEIZMO struct DATA over a series of 25 frequency bands from
%     80s to 8s.  Details on the default frequency bands can be found in
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
%       snrcutoff       - Default cutoff (10) for removing low SNR records
%       filterbank      - filter bank (use FILTER_BANK)
%       filtertype      - filter type (default is 'bp')
%       filterstyle     - filter style (default is 'butter')
%       filterorder     - filter order (default is 4)
%       filterpasses    - number of passes (default is 1)
%     The remaining options are passed to USERALIGN.  See that function for
%     more details.
%
%    Notes:
%     - The default filter bank is created using:
%       filter_bank([0.0125 0.125],'variable',0.2,0.1)
%
%    Examples:
%     % This is my typical usage form (for really nice quakes the upper
%     % limit of the filter bank can be raised to something like 0.2Hz):
%     bank=filter_bank([0.0125 0.125],'variable',0.2,0.1);
%     results=multibandalign(data,'bank',bank,'runname','examplerun');
%
%    See also: USERALIGN, USERSNR, FILTER_BANK, IIRFILTER

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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2011 at 12:25 GMT

% todo:

% check nargin
if(nargin<1)
    error('seizmo:multibandalign:notEnoughInputs',...
        'Not enough input arguments.');
elseif(nargin>1 && ~mod(nargin,2))
    error('seizmo:multbandalign:optionMustBePaired',...
        'Options must be paired with a value!');
end

% check data (dep)
versioninfo(data,'dep');

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
    [bank,runname,snrcut,filt,phase,varargin]=parse_mba_params(...
        varargin{:});
    
    % preallocate output
    info([])=struct('useralign',[],'filter',[],'initwin',[],...
        'usersnr',[],'userwinnow',[]);
    
    % require o field is set to the same value (within +/- 2 millisec)
    o=getheader(data,'o');
    outc=cell2mat(getheader(data,'o utc'));
    if(any(abs(timediff(outc(1,:),outc))>0.002))
        error('seizmo:multibandalign:oFieldVaries',...
            'O field must correspond to one UTC time for all records!');
    end
    
    % build lag table (for pre-aligning each band using previous lags)
    goodarr=-o;
    goodidx=1:nrecs;
    lags=o(:,ones(1,nrecs))'-o(:,ones(1,nrecs));
    idx=ones(nrecs);
    
    % loop over each filter
    i=1;       % starting with first band
    fast=true; % start in fast mode (usersnr & userwinnow only)
    last=[];   % for memory (last band aligned)
    while(i<=size(bank,1))
        % display filter info
        istr=num2str(i,'%02d');
        disp(['BAND ' istr ' FILTER CORNERS: ' ...
            num2str(1/bank(i,3)) 's  to  ' num2str(1/bank(i,2)) 's']);
        
        % get new relative arrivals & prealign
        % - lags are absolute
        % - use 10^idx as weight (higher weight to more recent freqs)
        %   - limit to >10^-5 for stability (-5 may be changed)
        % - use latest "good" arrivals to get absolute times
        tt=ttalign(lags,10.^(max(idx-i,-5)),goodarr,1,goodidx);
        data0=timeshift(data,-o-tt); % shift to origin then new times
        info(i,1).tt_start=tt;
        
        % filter the records
        data0=iirfilter(data0,filt.type,filt.style,'c',bank(i,2:3),...
            'o',filt.order,'p',filt.passes);
        filt.corners=bank(i,2:3);
        info(i).filter=filt;
        
        % initial window and assessment
        if(i>1); info(i).initwin=info(i-1).initwin; end
        [data0,info(i),skip,skipall,snrcut,fast,ax]=...
            get_initial_window(data0,info(i),snrcut,fast,phase);
        if(ishandle(ax))
            saveas(get(ax,'parent'),...
                [runname '_band_' istr '_preview.fig']);
            close(get(ax,'parent'));
        end
        
        % handle skipping
        if(skip); i=i+1; continue; end
        if(skipall); break; end
        
        % user specified snr
        if(isempty(last))
            info(i).usersnr.noisewin=[-300 -20];
            info(i).usersnr.signalwin=[-10 70];
            info(i).usersnr.method='peak2rms';
        elseif(last==i)
            % just update signal window if we used useralign
            if(~isempty(info(last).useralign) ...
                    && ~isempty(info(last).useralign.userwindow.limits))
                info(i).usersnr.signalwin=...
                    info(last).useralign.userwindow.limits;
            end
        else % look for zero-crossing
            info(i).usersnr.noisewin=info(last).usersnr.noisewin;
            if(~isempty(info(last).useralign) ...
                    && ~isempty(info(last).useralign.userwindow.limits))
                info(i).usersnr.signalwin=...
                    info(last).useralign.userwindow.limits;
            else
                info(i).usersnr.signalwin=info(last).usersnr.signalwin;
            end
            info(i).usersnr.method=info(last).usersnr.method;
            
            % use nearby zero crossing for new window limits
            zx=zeroxing(data0,info(i));
            info(i).usersnr.signalwin=[info(i).usersnr.signalwin(1) zx];
        end
        [snr,info(i).usersnr,ax]=usersnr(data0,...
            info(i).usersnr.noisewin,info(i).usersnr.signalwin,...
            info(i).usersnr.method,'normstyle','single');
        info(i).usersnr.snr=snr;
        info(i).usersnr.snrcut=snrcut;
        if(ishandle(ax))
            saveas(get(ax,'parent'),...
                [runname '_band_' istr '_usersnr.fig']);
            close(get(ax,'parent'));
        end
        
        % trimming dataset to those above snr cutoff
        snridx=find(snr>=snrcut);
        data0=data0(snridx);
        snr=snr(snridx);
        nn=numel(data0);
        
        % skip if less than 3
        if(nn<3)
            disp(['BAND ' istr ': Too Few High SNR data!']);
            continue;
        end
        
        % distance winnow
        if(isempty(last));
            info(i).userwinnow.limits=[];
        else
            info(i).userwinnow.limits=info(last).userwinnow.limits;
        end
        [data0,info(i).userwinnow,ax]=userwinnow(data0,...
            info(i).userwinnow.limits,...
            'normstyle','single',...
            'yfield','gcarc');
        snr(info(i).userwinnow.cut)=[];
        snridx(info(i).userwinnow.cut)=[];
        if(ishandle(ax))
            saveas(get(ax,'parent'),...
                [runname '_band_' istr '_userwinnow.fig']);
            close(get(ax,'parent'));
        end
        nn=numel(data0);
        
        % skip if less than 3
        if(nn<3)
            disp(['BAND ' istr ': Too few data in distance range!']);
            continue;
        end
        
        % fast vs control
        if(fast)
            % quiet align
            [info(i).useralign,info(i).useralign.xc,...
                info(i).useralign.data]=useralign_quiet(data0,...
                'spacing',1/(4*bank(i,3)),'absxc',false,...
                'estarr',zeros(nn,1),'snr',snr,...
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
            
            % plot alignment
            ax=plot2(info(i).useralign.data);
            if(ishandle(ax))
                saveas(get(ax,'parent'),...
                    [runname '_band_' istr '_aligned.fig']);
            end
            
            % ask if satisfied
            if(~user_satisfied)
                close(get(ax,'parent'));
                last=i;
                continue;
            else
                close(get(ax,'parent'));
            end
            
            % use alignment for next band
            if(i<size(bank,1))
                goodarr=info(i).useralign.solution.arr;
                goodidx=find(snridx);
                lags(goodidx,goodidx)=goodarr(:,ones(1,nn))...
                    -goodarr(:,ones(1,nn))';
                idx(goodidx,goodidx)=i+1;
            end
        else
            try
                % detailed alignment
                [info(i).useralign,info(i).useralign.xc,...
                    info(i).useralign.data]=useralign(data0,...
                    'spacing',1/(4*bank(i,3)),'absxc',false,...
                    'estarr',zeros(nn,1),'snr',snr,...
                    'window',info(i).usersnr.signalwin,varargin{:});
                if(any(ishandle(info(i).useralign.handles)))
                    ax=info(i).useralign.handles(...
                        ishandle(info(i).useralign.handles));
                    saveas(get(ax(1),'parent'),...
                        [runname '_band_' istr '_useralign.fig']);
                    close(get(ax(1),'parent'));
                end
            catch
                % ask if user is satisfied with breakage
                if(~user_satisfied)
                    continue;
                else
                    i=i+1;
                    continue;
                end
            end
            
            % amplitude analysis
            % - amplitude is the peak2peak
            % - amplitude "standard error" is given my amp/snr
            %   where snr is the peak2peak amplitude of the signal
            %   to the rms of the noise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THIS IS IF RMS OF THE SIGNAL IS USED FOR SNR
            %scale=getvaluefun(info(i).useralign.data,@(x)sqrt(mean(x.^2)));
            %info(i).useralign.data=divide(info(i).useralign.data,scale);
            %info(i).useralign.solution.amp=scale;
            %info(i).useralign.solution.amperr=scale./snr;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THIS IS IF PEAK2RMS IS USED FOR SNR
            %scale=(gh(info(i).useralign.data,'depmax')...
            %    -gh(info(i).useralign.data,'depmin'))/2;
            %info(i).useralign.data=divide(info(i).useralign.data,scale);
            %info(i).useralign.solution.amp=scale;
            %info(i).useralign.solution.amperr=scale./snr;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THIS IS THE MOST VISUALLY APPEALING
            % BUT ERRORS TEND TO BE OVERESTIMATED BY ~10%
            [info(i).useralign.data,scale]=normalize(...
                info(i).useralign.data);
            info(i).useralign.solution.amp=scale;
            info(i).useralign.solution.amperr=scale./snr;
            
            % ask to utilize for subsequent bands
            if(i<size(bank,1) && use_align)
                goodarr=info(i).useralign.solution.arr;
                goodidx=find(snridx);
                lags(goodidx,goodidx)=goodarr(:,ones(1,nn))...
                    -goodarr(:,ones(1,nn))';
                idx(goodidx,goodidx)=i+1;
            end
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


function [lgc]=use_align()
happy_user=false;
while(~happy_user)
    choice=menu('Use alignment for subsequent filters?','YES','NO');
    switch choice
        case 1
            lgc=true;
            happy_user=true;
        case 2
            lgc=false;
            happy_user=true;
    end
end
end


function [lgc]=user_satisfied()
happy_user=false;
while(~happy_user)
    choice=menu('Satisfied with the alignment?','YES','NO');
    switch choice
        case 1
            lgc=true;
            happy_user=true;
        case 2
            lgc=false;
            happy_user=true;
    end
end
end


function [data,info,skip,skipall,snrcut,fast,ax]=get_initial_window(...
    data,info,snrcut,fast,phase)

% default initial window
if(isempty(info.initwin))
    win=[-300 300];
else
    win=info.initwin;
end

% default modestr
if(fast)
    modestr='POWERUSER';
else
    modestr='FAST';
end

% shift if phase given
if(~isempty(phase))
    [t,n]=getarrival(data,phase);
    data=timeshift(data,-t,strcat('it',num2str(n)));
end

% ask the user if they wish to continue
skip=false; skipall=false;
happy_user=false; ax=-1;
while(~happy_user)
    % plot records (so user can decide to skip or not)
    if(~ishandle(ax))
        % record section
        ax=recordsection(cut(data,win(1),win(2)),...
            'normstyle','single','xlim',win,'title',...
            ['FILTER CORNERS: ' num2str(1/info.filter.corners(2)) ...
            's  to  ' num2str(1/info.filter.corners(1)) 's']);
        
        % this all assumes the phase is shifted to 0
        if(~isempty(phase))
            evdp=getheader(data(1),'evdp')/1000;
            tt=taupcurve('dep',evdp);
            idx=find(strcmp(phase,{tt.phase}));
            intrcpt=tt(idx).time(1)...
                -tt(idx).distance(1)*tt(idx).rayparameter(1);
            tt=taupcurve('dep',evdp,...
                         'reddeg',1/tt(idx).rayparameter(1),'ph','ttall');
            hold(ax,'on');
            h=plot_taupcurve(tt,-intrcpt,true,'parent',ax,'linewidth',5);
            movekids(h,'back');
            hold(ax,'off');
        end
    end

    choice=menu('What to do with this filter band?',...
        'Align The Waveforms!',...
        'Adjust Initial Window',...
        ['Adjust SNR Cutoff (' num2str(snrcut) ')'],...
        ['Change To ' modestr ' Mode'],...
        'Skip This Filter Band',...
        'Skip All Remaining Filter Bands');
    switch choice
        case 1 % process
            happy_user=true;
        case 2 % adjust window
            % get new initial window
            if(ishandle(ax)); close(get(ax,'parent')); end
            [win,win,ax]=userwindow(data,win);
            if(ishandle(ax)); close(get(ax,'parent')); end
            
            % default initial window if none
            if(isempty(win.limits))
                if(isempty(info.initwin))
                    win=[-300 300];
                else
                    win=info.initwin;
                end
            else
                win=win.limits;
            end
        case 3 % snrcut
            tmp=inputdlg(...
                ['SNR CutOff? [' num2str(snrcut) ']:'],...
                'Enter SNR CutOff',1,{num2str(snrcut)});
            if(~isempty(tmp))
                try
                    tmp=str2double(tmp{:});
                    if(~isnan(tmp) && ~isinf(tmp))
                        snrcut=tmp;
                    end
                catch
                    % do not change snrcut
                end
            end
        case 4 % switch mode
            if(fast)
                fast=false;
                modestr='FAST';
            else
                fast=true;
                modestr='POWERUSER';
            end
        case 5 % skip band
            happy_user=true;
            skip=true;
        case 6 % skip all
            happy_user=true;
            skipall=true;
    end
end

% implement window
info.initwin=win;
if(choice==1); data=cut(data,win(1),win(2)); end

end


function [bank,run,snrcut,filt,phase,varargin]=parse_mba_params(varargin)
%PARSE_MBA_PARAMS  Parses out multibandalign parameters

% basic checks on optional inputs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:multibandalign:badInput',...
        'All OPTIONs must be specified with a string!');
end

% defaults
bank=filter_bank([0.0125 0.125],'variable',0.2,0.1);
run=['multibandalign_' datestr(now,30)];
snrcut=3;
filt.type='bandpass';
filt.style='butter';
filt.order=4;
filt.passes=1;
phase=[];

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
            bank=varargin{i+1};
            keep(i:i+1)=false;
        case {'rn' 'run' 'name' 'runname'}
            if(~isstring(varargin{i+1}))
                error('seizmo:multibandalign:badInput',...
                    'RUNNAME must be a string!');
            end
            run=varargin{i+1};
            keep(i:i+1)=false;
        case {'sc' 'snrcutoff' 'snrcut'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}<0)
                error('seizmo:multibandalign:badInput',...
                    'SNRCUTOFF must be a real-valued scalar!');
            end
            snrcut=varargin{i+1};
            keep(i:i+1)=false;
        case {'ft' 'filtertype'}
            if(~isstring(varargin{i+1}) ...
                    || ~any(strcmpi(varargin{i+1},validtypes)))
                error('seizmo:multibandalign:badInput',...
                    'FILTERTYPE must be a filter type in IIRDESIGN!');
            end
            filt.type=varargin{i+1};
            keep(i:i+1)=false;
        case {'fs' 'filterstyle'}
            if(~isstring(varargin{i+1}) ...
                    || ~any(strcmpi(varargin{i+1},validstyles)))
                error('seizmo:multibandalign:badInput',...
                    'FILTERSTYLE must be a filter style in IIRDESIGN!');
            end
            filt.style=varargin{i+1};
            keep(i:i+1)=false;
        case {'fo' 'filterorder'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}~=fix(varargin{i+1}) ...
                    || varargin{i+1}<=0)
                error('seizmo:multibandalign:badInput',...
                    'FILTERORDER must be a scalar integer >0!');
            end
            filt.order=varargin{i+1};
            keep(i:i+1)=false;
        case {'fp' 'filterpasses' 'filterpass'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                    || varargin{i+1}~=fix(varargin{i+1}) ...
                    || ~any(abs(varargin{i+1})==[1 2]))
                error('seizmo:multibandalign:badInput',...
                    'FILTERPASSES must be 1 or 2!');
            end
            filt.passes=varargin{i+1};
            keep(i:i+1)=false;
        case {'phase'}
            if(~isstring(varargin{i+1}) ...
                    || ~any(strcmpi(varargin{i+1},{'Pdiff' 'Sdiff'})))
                error('seizmo:multibandalign:badInput',...
                    'PHASE must be a ''Pdiff'' or ''Sdiff''!');
            end
            phase=varargin{i+1};
            keep(i:i+1)=false;
    end
end
varargin=varargin(keep);

end


function [zx]=zeroxing(data,info)
% get sample spacing
delta=getheader(data,'delta');
delta=unique(delta);
if(~isscalar(delta))
   error('seizmo:multibandalign:badInput',...
       'DATA records must have the same sample rate!');
end

% get stack
dstack=stack(data,delta,[],info.initwin(1),info.initwin(2));

% get zero crossings
zx=zerocrossings(dstack);
zx{1}=sort(zx{1}); % just in case

% find closest within range
sw=info.usersnr.signalwin;
neg=find(zx{1}<sw(2) & zx{1}>=(sw(1)+sw(2))/2,1,'last');
pos=find(zx{1}>=sw(2) & zx{1}<=(sw(2)+(sw(1)+sw(2))/2),1,'first');

% choose one
if(isempty(neg) && isempty(pos))
    % use previous
    zx=sw(2);
elseif(isempty(neg))
    zx=zx{1}(pos);
elseif(isempty(pos))
    zx=zx{1}(neg);
else
    % move positive only if < 1/3 negative distance
    if((sw(2)-zx{1}(neg))/3<=(zx{1}(pos)-sw(2)))
        zx=zx{1}(neg);
    else
        zx=zx{1}(pos);
    end
end

end


