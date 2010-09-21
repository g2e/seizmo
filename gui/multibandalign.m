function [info]=multibandalign(data,bank,varargin)
%MULTIBANDALIGN    Aligns signal at multiple frequency bands
%
%    Usage:    info=multibandalign(data,bank)
%              info=multibandalign(data,bank,'option',value,...)
%
%    Description:
%     INFO=MULTIBANDALIGN(DATA,BANK) presents the user with an interface to
%     align the data given in SEIZMO struct DATA at the frequency bands
%     specified in the filter bank BANK.  BANK should be formatted as
%     output from FILTER_BANK.  DATA should be records where the phase of
%     interest has already been isolated (not necessary but saves you from
%     having to do it for every frequency band).  You can do this using a
%     combination of USERMOVEOUT, USERWINDOW, and USERTAPER.  The output is
%     a struct containing the USERALIGN output from each aligned frequency
%     as well as SNR info and filter info.  All figures are automatically
%     saved to the current directory as .fig files.
%
%     INFO=MULTIBANDALIGN(DATA,BANK,'OPTION',VALUE,...) passes additional
%     options to USERALIGN.  See that function for details.
%
%    Notes:
%
%    Examples:
%     % 
%
%    See also: USERALIGN, USERSNR, FILTER_BANK, IIRFILTER

%     Version History:
%        Mar. 25, 2010 - initial version
%        Sep. 21, 2010 - working version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 21, 2010 at 12:25 GMT

% todo:

% check nargin
if(nargin<2)
    error('seizmo:multibandalign:notEnoughInputs',...
        'Not enough input arguments.');
elseif(nargin>4 && mod(nargin,2))
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
    error(lasterror)
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
    
    % preallocate output
    info([])=struct('info',[],'xc',[],'data',[],'filter',[],'snr',[]);
    
    % require o field is set to the same value (within +/- 2 millisec)
    o=getheader(data,'o');
    outc=cell2mat(getheader(data,'o utc'));
    if(any(abs(timediff(outc(1,:),outc))>0.002))
        error('seizmo:multibandalign:oFieldVaries',...
            'O field must correspond to one UTC time for all records!');
    end
    
    % get additional filter info (default to 4th order 1-pass butterworth)
    filt=bandpass_parameters('butter',4,1);
    
    % get snr cut (default to 3)
    snrcut=get_snr_cutoff(3);
    
    % build lag table (for pre-aligning each band using previous lags)
    lags=o(:,ones(1,nrecs))'-o(:,ones(1,nrecs));
    idx=ones(nrecs);
    goodarr=-o;
    goodidx=1:nrecs;
    
    % loop over each filter
    for i=1:size(bank,1)
        % get new relative arrivals & prealign
        % - lags are absolute
        % - use 10^idx as weight (higher weight to more recent freqs)
        % - use latest "good" arrivals to get absolute times
        tt=ttalign(lags,10.^(idx-i),goodarr,1,goodidx);
        data0=timeshift(data,-o-tt); % shift to origin then new times
        info(i).tt_start=tt;
        
        % filter the records
        data0=iirfilter(data0,'bp',filt.style,'c',bank(i,2:3),...
            'o',filt.order,'p',filt.passes);
        filt.corners=bank(i,2:3);
        info(i).filter=filt;
        
        % plot records (so user can decide to skip or not)
        ax=plot0(data0,'title',...
            ['FILTER CORNERS: ' num2str(1/bank(i,3)) 's  to  ' ...
            num2str(1/bank(i,2)) 's']);
        
        % ask the user if they wish to continue
        skip=false; skipall=false;
        happy_user=false;
        while(~happy_user)
            choice=menu('What to do with this filter band?',...
                'Process','Skip','Skip All Remaining');
            switch choice
                case 1
                    happy_user=true;
                case 2
                    happy_user=true;
                    skip=true;
                case 3
                    happy_user=true;
                    skipall=true;
            end
        end
        
        % handle skipping
        istr=num2str(i);
        saveas(get(ax,'parent'),['multibandalign_' istr '_preview.fig']);
        close(get(ax,'parent'));
        if(skip); continue; end
        if(skipall); break; end
        
        % user estimated snr
        [snr,s,ax]=usersnr(data0);
        info(i).snr=s;
        info(i).snr.snrcut=snrcut;
        info(i).snr.snr=snr;
        saveas(get(ax,'parent'),['multibandalign_' istr '_usersnr.fig']);
        close(get(ax,'parent'));
        
        % trimming dataset to snr significant
        snridx=snr>=snrcut;
        data0=data0(snridx);
        snr=snr(snridx);
        nn=numel(data0);
        
        % align
        [tmp,xc,data1]=useralign(data0,...
            'estarr',zeros(nn,1),'snr',snr,varargin{:});
        info(i).info=tmp;
        info(i).xc=xc;
        info(i).data=data1;
        clear tmp xc data1
        ax=info.handles(ishandle(info.handles));
        saveas(get(ax(1),'parent'),['multibandalign_' istr '_useralign.fig']);
        close(get(ax(1),'parent'));
        
        % amplitude analysis
        % - need peak2peak amplitude
        % - amplitude "standard error" is given my amp/(2*snr)
        
        % ask to utilize for subsequent bands
        if(use_align)
            goodarr=info(i).info.solution.arr;
            goodidx=find(snridx);
            lags(goodidx,goodidx)=...
                info(i).info.solution.arr(:,ones(1,nn))'...
                -info.solution.arr(:,ones(1,nn));
            idx(goodidx,goodidx)=i;
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
    error(lasterror)
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
    choice=menu('Utilize alignment for subsequent filters?','YES','NO');
    switch choice
        case 1
            lgc=true;
        case 2
            lgc=false;
    end
end
end


