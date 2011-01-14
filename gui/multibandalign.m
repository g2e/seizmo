function [info]=multibandalign(data,bank,runname,varargin)
%MULTIBANDALIGN    Aligns signal at multiple frequency bands
%
%    Usage:    results=multibandalign(data,bank,runname)
%              results=multibandalign(data,bank,runname,'option',value,...)
%
%    Description:
%     RESULTS=MULTIBANDALIGN(DATA,BANK,RUNNAME) presents an interface
%     to align the data given in SEIZMO struct DATA at the frequency bands
%     specified in the filter bank BANK.  BANK should be formatted as
%     output from FILTER_BANK.  DATA should be records where the phase of
%     interest has already been isolated (not necessary but saves you from
%     having to do it for every frequency band).  You can do this using a
%     combination of USERMOVEOUT, USERWINDOW, and USERTAPER.  The output is
%     a struct containing the USERALIGN output from each aligned frequency
%     as well as SNR info and filter info.  All figures are automatically
%     saved to the current directory as .fig files.
%
%     RESULTS=MULTIBANDALIGN(DATA,BANK,RUNNAME,'OPTION',VALUE,...) passes
%     options to USERALIGN.  See that function for details.
%
%    Notes:
%
%    Examples:
%     % This is my typical usage form (for really nice quakes the upper
%     % limit of the filter bank can be raised to something like 0.2Hz):
%     bank=filter_bank([0.0125 0.125],'variable',0.2,0.1);
%     results=multibandalign(data,bank,'examplerun');
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 14, 2011 at 12:25 GMT

% todo:

% check nargin
if(nargin<3)
    error('seizmo:multibandalign:notEnoughInputs',...
        'Not enough input arguments.');
elseif(nargin>3 && ~mod(nargin,2))
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
    % check filter bank
    if(size(bank,2)~=3 ...
            || any(bank(:)<=0 | isnan(bank(:)) | isinf(bank(:))) ...
            || any(bank(:,1)<=bank(:,2) | bank(:,3)<=bank(:,1)))
        error('seizmo:multibandalign:badInput',...
            'BANK must be in the format from FILTER_BANK!');
    end
    
    % check runname
    if(~isstring(runname))
        error('seizmo:multibandalign:badInput',...
            'RUNNAME must be a string!');
    end
    
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
    
    % get additional filter info (default to 4th order 1-pass butterworth)
    filt=bandpass_parameters('butter',4,1);
    
    % get snr cut (default to 10)
    snrcut=get_snr_cutoff(10);
    
    % build lag table (for pre-aligning each band using previous lags)
    lags=o(:,ones(1,nrecs))'-o(:,ones(1,nrecs));
    idx=ones(nrecs);
    goodarr=-o;
    goodidx=1:nrecs;
    
    % loop over each filter
    i=1;       % starting with first band
    fast=true; % start in fast mode (usersnr & userwinnow only)
    while(i<=size(bank,1))
        % display filter info
        disp(['BAND ' num2str(i,'%02d') ' FILTER CORNERS: ' ...
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
        data0=iirfilter(data0,'bp',filt.style,'c',bank(i,2:3),...
            'o',filt.order,'p',filt.passes);
        filt.corners=bank(i,2:3);
        info(i).filter=filt;
        
        % initial window and assessment
        if(i>1); info(i).initwin=info(i-1).initwin; end
        [data0,info(i),skip,skipall,fast,ax]=...
            get_initial_window(data0,info(i),fast);
        
        % handle skipping
        istr=num2str(i);
        saveas(get(ax,'parent'),[runname '_band_' istr '_preview.fig']);
        close(get(ax,'parent'));
        if(skip); i=i+1; continue; end
        if(skipall); break; end
        
        % user estimated snr
        if(~exist('s','var'))
            s.noisewin=[-300 -20];
            s.signalwin=[-10 70];
        end
        [snr,s,ax]=usersnr(data0,s.noisewin,s.signalwin,'peak2rms',...
            'normstyle','single');
        info(i).usersnr=s;
        info(i).usersnr.snrcut=snrcut;
        info(i).usersnr.snr=snr;
        saveas(get(ax,'parent'),[runname '_band_' istr '_usersnr.fig']);
        close(get(ax,'parent'));
        
        % trimming dataset to snr significant
        snridx=snr>=snrcut;
        data0=data0(snridx);
        snr=snr(snridx);
        nn=numel(data0);
        
        % skip if less than 3
        if(nn<3)
            disp([runname ' band ' num2str(istr) ': Too Few High SNR data!']);
            continue;
        end
        
        % distance winnow
        [data0,info(i).userwinnow,ax]=userwinnow(data0,'normstyle','single');
        snr(info(i).userwinnow.cut)=[];
        tmp=find(snridx);
        snridx(tmp(info(i).userwinnow.cut))=false;
        if(ishandle(ax))
            saveas(get(ax,'parent'),[runname '_band_' istr '_userwinnow.fig']);
            close(get(ax,'parent'));
        end
        nn=numel(data0);
        
        % skip if less than 3
        if(nn<3)
            disp([runname ' band ' num2str(istr) ...
                ': Too few data in distance range!']);
            continue;
        end
        
        % fast vs control
        if(fast)
            % quiet align
            [info(i).useralign,info(i).useralign.xc,...
                info(i).useralign.data]=useralign_quiet(data0,...
                'spacing',1/(4*bank(i,3)),'absxc',false,...
                'estarr',zeros(nn,1),'snr',snr,'window',s.signalwin,...
                varargin{:});
            
            % Get Amplitude Info
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THIS IS THE MOST VISUALLY APPEALING
            % BUT ERRORS TEND TO BE OVERESTIMATED BY ~10%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [info(i).useralign.data,scale]=normalize(info(i).useralign.data);
            info(i).useralign.solution.amp=scale;
            info(i).useralign.solution.amperr=scale./snr;
            
            % plot alignment
            ax=plot2(info(i).useralign.data);
            if(ishandle(ax))
                saveas(get(ax,'parent'),[runname '_band_' istr '_aligned.fig']);
            end
            
            % ask if satisfied
            if(~user_satisfied)
                close(get(ax,'parent'));
                continue;
            else
                close(get(ax,'parent'));
            end
            
            % use alignment for next band?
            if(i<size(bank,1))
                goodarr=info(i).useralign.solution.arr;
                goodidx=find(snridx);
                lags(goodidx,goodidx)=goodarr(:,ones(1,nn))-goodarr(:,ones(1,nn))';
                idx(goodidx,goodidx)=i+1;
            end
            i=i+1;
        else
            % align
            try
                [info(i).useralign,info(i).useralign.xc,...
                    info(i).useralign.data]=useralign(data0,...
                    'spacing',1/(4*bank(i,3)),'absxc',false,...
                    'estarr',zeros(nn,1),'snr',snr,varargin{:});
                ax=info(i).useralign.handles(ishandle(info(i).useralign.handles));
                saveas(get(ax(1),'parent'),[runname '_band_' istr '_useralign.fig']);
                close(get(ax(1),'parent'));
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
            [info(i).useralign.data,scale]=normalize(info(i).useralign.data);
            info(i).useralign.solution.amp=scale;
            info(i).useralign.solution.amperr=scale./snr;
            
            % ask if user is satisfied
            if(~user_satisfied)
                continue;
            end
            
            % ask to utilize for subsequent bands
            if(i<size(bank,1) && use_align)
                goodarr=info(i).useralign.solution.arr;
                goodidx=find(snridx);
                lags(goodidx,goodidx)=goodarr(:,ones(1,nn))-goodarr(:,ones(1,nn))';
                idx(goodidx,goodidx)=i+1;
            end
            i=i+1;
        end
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


function [filt]=bandpass_parameters(style,order,passes)

% set defaults
filt.style=style;
filt.order=order;
filt.passes=passes;

happy_user=false;
while(~happy_user)
    choice=menu('ADJUST FILTER OPTIONS?',...
        ['FILTER STYLE: ' upper(filt.style)],...
        ['FILTER ORDER: ' num2str(filt.order)],...
        ['NUMBER OF PASSES: ' num2str(filt.passes)],...
        'NO');
    switch choice
        case 1 % style
            choice=menu('CHOOSE A STYLE',...
                ['CURRENT (' upper(filt.style) ')'],...
                'BUTTERWORTH','CHEBYSHEV TYPE 1','ELLIPTIC');
            switch choice
                case 2
                    filt.style='butter';
                case 3
                    filt.style='cheby1';
                case 4
                    filt.style='ellip';
            end
        case 2 % order
            tmp=inputdlg(...
                ['Filter Order? [' num2str(filt.order) ']:'],...
                'Enter Filter Order',1,{num2str(filt.order)});
            if(~isempty(tmp))
                try
                    tmp=str2double(tmp{:});
                    if(~isnan(tmp) && ~isinf(tmp) ...
                            && tmp>0 && tmp==fix(tmp))
                        filt.order=tmp;
                    end
                catch
                    % do not change filt.order
                end
            end
        case 3 % passes
            choice=menu('CHOOSE NUMBER OF PASSES',...
                ['CURRENT (' num2str(filt.passes) ')'],...
                '1','2');
            switch choice
                case 2
                    filt.passes=1;
                case 3
                    filt.passes=2;
            end
        case 4 % done
            happy_user=true;
    end
end

end


function [snrcut]=get_snr_cutoff(snrcut)
happy_user=false;
while(~happy_user)
    choice=menu(['Adjust SNR CutOff (' num2str(snrcut) ') ?'],'YES','NO');
    switch choice
        case 1
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
        case 2
            happy_user=true;
    end
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


function [data,info,skip,skipall,fast,ax]=get_initial_window(data,info,fast)

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

% ask the user if they wish to continue
skip=false; skipall=false;
happy_user=false;
while(~happy_user)
    % plot records (so user can decide to skip or not)
    ax=plot0(cut(data,win(1),win(2)),...
        'normstyle','single','xlim',win,'title',...
        ['FILTER CORNERS: ' num2str(1/info.filter.corners(2)) 's  to  ' ...
        num2str(1/info.filter.corners(1)) 's']);

    choice=menu('What to do with this filter band?',...
        'Align The Waveforms!',...
        'Adjust Initial Window',...
        ['Change To ' modestr ' Mode'],...
        'Skip This Filter Band',...
        'Skip All Remaining Filter Bands');
    switch choice
        case 1 % process
            happy_user=true;
        case 2 % adjust window
            try
                close(get(ax,'parent'));
                [win,win,ax]=userwindow(data,win);
                close(get(ax,'parent'));
                win=win.limits;
            catch
            end
        case 3 % switch mode
            if(fast)
                fast=false;
                modestr='FAST';
            else
                fast=true;
                modestr='POWERUSER';
            end
        case 4 % skip band
            happy_user=true;
            skip=true;
        case 5 % skip all
            happy_user=true;
            skipall=true;
    end
end

% implement window
info.initwin=win;
if(choice<3); data=cut(data,win(1),win(2)); end

end
