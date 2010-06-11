function []=multibandalign(data,bank,projname,varargin)
%MULTIBANDALIGN    Aligns signal at multiple frequency bands
%
%    Usage:    []=multibandalign(data,bank,projname)
%              []=multibandalign(data,bank,projname,...,useralign_options,...)
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also: USERALIGN

%     Version History:
%        Mar. 25, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2010 at 12:25 GMT

% todo:
% - menu to utilize alignment for subsequent bands

% check nargin
error(nargchk(3,inf,nargin));

% check data (dep)
versioninfo(data,'dep');

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
    % check required inputs (bank, projname)
    if(size(bank,2)~=3 ...
            || any(bank(:)<=0 | isnan(bank(:)) | isinf(bank(:))) ...
            || any(bank(:,1)<=bank(:,2) | bank(:,3)<=bank(:,1)))
        error('seizmo:multibandalign:badInput',...
            'BANK must be in the format from FILTER_BANK!');
    end
    if(~ischar(projname) || size(projname,1)~=1)
        error('seizmo:multibandalign:badInput',...
            'PROJNAME must be a string!');
    end
    
    % fix projname
    [projpath,name,ext]=fileparts(projname);
    if(isempty(ext)); ext='.mat'; end
    projname=fullfile(projpath,name);
    
    % save inputs to a file
    save([projname ext],'data','bank','projname','varargin');
    
    % require o field is set to the same value (within +/- 2 millisec)
    o=getheader(data,'o');
    outc=cell2mat(getheader(data,'o utc'));
    if(any(abs(timediff(outc(1,:),outc))>0.002))
        error('seizmo:multibandalign:oFieldVaries',...
            'O field must correspond to one UTC time for all records!');
    end
    
    % get additional filter info
    filt=bandpass_parameters('butter',4,1);
    
    % get snr cut (default to 3)
    snrcut=get_snr_cutoff(3);
    
    % number of records
    nrecs=numel(data);
    
    % build align table (for pre-aligning each band using previous align)
    lags=o(:,ones(1,nrecs))'-o(:,ones(1,nrecs));
    idx=zeros(nrecs);
    goodarr=-o;
    goodidx=1:nrecs;
    
    % loop over each filter
    for i=1:size(bank,1)
        % get new relative arrivals & prealign
        % - lags are absolute
        % - use 10^idx as weight
        % - use latest "good" arrivals to get absolute times
        tt=ttalign(lags,10.^idx,goodarr,1,goodidx);
        data0=timeshift(data,-o-tt);
        
        % filter the records
        data0=iirfilter(data0,'bp',filt.style,'c',bank(i,2:3),...
            'o',filt.order,'p',filt.passes);
        filt.corners=bank(i,2:3);
        
        % plot records
        fh0=p0(data0,'title',['FILTER: ' num2str(bank(i,2)) 'Hz  to  ' ...
            num2str(bank(i,3)) 'Hz']);
        
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
        if(skip); continue; end
        if(skipall); break; end
        
        % user estimated snr
        [snr,s,fh]=usersnr(data0);
        
        % trimming dataset to snr significant
        snridx=snr>=snrcut;
        data0=data0(snridx);
        snr=snr(snridx);
        nn=numel(data0);
        
        % align
        [info,xc,data1]=useralign(data0,...
            'estarr',zeros(nn,1),'snr',snr,varargin{:});
        
        % add snr params to info
        info.usersnr=s;
        info.figurehandles(end+1)=fh;
        info.snrcut.value=snrcut;
        info.snrcut.snridx=snridx;
    
        % perform cluster analysis on results
        [info.usercluster,info.figurehandles(end+1)]=usercluster(data1,...
            xc.cg(:,:,1));
        
        % amplitude analysis
        % - need peak2peak amplitude
        % - amplitude "standard error" is given my amp/(2*snr)
        
        % ask to utilize for subsequent bands
        if(use_align)
            % biggest cluster only
            pop=histc(info.usercluster.T,1:max(info.usercluster.T));
            [idx,idx]=max(pop);
            biggestclusteridx=idx==info.usercluster.T;
            goodidx=find(snridx);
            goodidx=goodidx(biggestclusteridx);
            goodarr=info.solution.arr(biggestclusteridx);
            lags(goodidx,goodidx)=info.solution.arr(:,ones(1,nn))'...
                -info.solution.arr(:,ones(1,nn));
            idx(goodidx,goodidx)=i;
        end
        
        % save everything (add in filter info, align table too!)
        info.iirfilter=filt;
        save([projname '_' num2str(i,'%02d') ext],...
            'info','xc','data1',...
            'goodidx','goodarr','lags','idx');
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


