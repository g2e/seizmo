function [varargout]=stft(data,varargin)
%STFT    Short Time Fourier Transform (aka Sliding Spectra) of SEIZMO data
%
%    Usage:    stft(data)
%              stft(data,...,'width',seconds,...)
%              stft(data,...,'units',type,...)
%              stft(data,...,'overlap',percent,...)
%              stft(data,...,'pow2pad',power,...)
%              stft(data,...,'colormap',cmap,...)
%              stft(data,...,'freqrange',frange,...)
%              stft(data,...,'dbrange',dbrange,...)
%              stft(data,...,'colorbar',show,...)
%              stft(data,...,'axis',handle,...)
%              data=stft(...)
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also: FTAN

%     Version History:
%        Apr. 26, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 26, 2010 at 19:15 GMT

% todo:
% options
% - window width in seconds
% - overlap in %
% - pow2pad
% plot vs output
% - no plot if output
% - plot options (plot 1 style)
%   - colormap
%   - db range
%   - freq range
%   - hotrod

% check nargin
msg=nargchk(1,inf,nargin);
if(~isempty(msg)); error(msg); end
if(mod(nargin-1,2))
    error('seizmo:stft:badNumInputs',...
        'Unpaired Option/Value!');
end

% check data structure
versioninfo(data,'dep');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt spectrogram
try
    % check headers
    data=checkheader(data);
    
    % verbosity
    verbose=seizmoverbose;
    
    % number of records
    nrecs=numel(data);
    
    % require evenly spaced time series
    iftype=getenumid(data,'iftype');
    leven=getlgc(data,'leven');
    if(any(strcmpi(leven,'false')))
        error('seizmo:stft:badLEVEN',...
            ['Record(s):\n' sprintf('%d ',find(strcmpi(leven,'false'))) ...
            '\nInvalid operation on unevenly sampled record(s)!']);
    elseif(any(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy')))
        error('seizmo:stft:badIFTYPE',...
            ['Record(s):\n' sprintf('%d ',...
            find(~strcmpi(iftype,'itime') & ~strcmpi(iftype,'ixy'))) ...
            '\nDatatype of record(s) in DATA must be Timeseries or XY!']);
    end
    
    % valid string option values
    valid.units={'%' 'per' 'percent' ...
                 'n' 'samp' 'samples' ...
                 's' 'sec' 'seconds'};
    
    % option defaults
    varargin=[{'u' '%' 'w' 10 'o' 75 'p' 1 'c' 'fire' 'h' [2/3 1] ...
        'd' [] 'f' []} varargin];
    
    % check options
    if(~iscellstr(varargin(1:2:end)))
        error('seizmo:stft:badOption',...
            'All Options must be specified with a string!');
    end
    for i=1:2:numel(varargin)
        % skip empty
        skip=false;
        if(isempty(varargin{i+1})); skip=true; end
        
        % check option is available
        switch lower(varargin{i})
            case {'u' 'unit' 'units'}
                if(skip); continue; end
                if(ischar(varargin{i+1}))
                    varargin{i+1}=cellstr(varargin{i+1});
                end
                if(~iscellstr(varargin{i+1}) ...
                        || any(~ismember(varargin{i+1},valid.units)) ...
                        || ~any(numel(varargin{i+1})==[1 nrecs]))
                    error('seizmo:stft:badUnits',...
                        'UNITS is fucked');
                end
                units=varargin{i+1};
                if(isscalar(units)); units(1:nrecs,1)=units; end
            case {'w' 'window' 'width' 'length' 'len'}
                if(skip); continue; end
                if(~isreal(varargin{i+1}) || any(varargin{i+1}<=0))
                    error('seizmo:stft:badWidth',...
                        'WIDTH is fucked');
                end
                width=varargin{i+1};
                if(isscalar(width)); width(1:nrecs,1)=width; end
            case {'o' 'over' 'olap' 'overlap'}
                if(skip); continue; end
                %
                overlap=varargin{i+1};
                if(isscalar(overlap)); overlap(1:nrecs,1)=overlap; end
            case {'p' 'p2p' 'padpower' 'pow2pad'}
                if(skip); continue; end
                %
                pow2pad=varargin{i+1};
                if(isscalar(pow2pad)); pow2pad(1:nrecs,1)=pow2pad; end
            case {'c' 'cmap' 'color' 'colormap'}
                if(skip); continue; end
                if(ischar(varargin{i+1}))
                    varargin{i+1}=cellstr(varargin{i+1});
                end
                %
                cmap=varargin{i+1};
                if(isscalar(cmap)); cmap(1:nrecs,1)=cmap; end
            case {'z' 'zlim' 'h' 'hr' 'hotrod'}
                if(skip); continue; end
                %
                hotrod=varargin{i+1};
                if(size(hotrod,1)==1); hotrod=hotrod(ones(nrecs,1),:); end
            case {'d' 'db' 'dbr' 'dbrange'}
                if(skip); continue; end
                %
                dbrange=varargin{i+1};
                if(size(dbrange,1)==1); dbrange=dbrange(ones(nrecs,1),:); end
            case {'f' 'fr' 'freq' 'frange' 'freqrange'}
                if(skip); continue; end
                %
                frange=varargin{i+1};
                if(size(frange,1)==1); frange=frange(ones(nrecs,1),:); end
            otherwise
                error('seizmo:stft:badOption',...
                    'Unknown Option: %s',varargin{i});
        end
    end
    
    % identify width units
    pw=ismember(units,{'%' 'per' 'percent'});
    nw=ismember(units,{'n' 'samp' 'samples'});
    sw=ismember(units,{'s' 'sec' 'seconds'});
    
    % convert width/overlap/pow2pad to proper units (samples)
    [b,delta,npts,ncmp]=getheader(data,'b','delta','npts','ncmp');
    if(any(pw)); width(pw)=ceil(width(pw).*npts(pw)/100); end
    if(any(nw)); width(nw)=ceil(width(nw)); end
    if(any(sw)); width(sw)=ceil(width./delta); end
    overlap=ceil(overlap.*width/100);
    pow2pad=2.^(nextpow2n(width)+pow2pad);
    
    % detail message
    if(verbose)
        disp('Getting Short Time Fourier Transform of Record(s)');
        print_time_left(0,nrecs);
    end
    
    % loop over records
    e=nan(nrecs,1); nxsize=e; nysize=e; depmen=e; depmin=e; depmax=e;
    for i=1:nrecs
        % skip dataless
        if(~npts(i)); continue; end
        
        % return spectrograms
        if(nargout)
            % get spectrogram
            P=cell(ncmp(i),1);
            for j=1:ncmp(i)
                [P{j},F,T,P{j}]=spectrogram(double(data(i).dep(:,j)),...
                    width(i),overlap(i),pow2pad(i),1/delta(i));
                P{j}=P{j}(:); % make column vector
            end
            
            % assign power spectra to dep
            data(i).dep=cell2mat(P);
            
            % get fields
            % Differences from SAC (!!!):
            % - b/e/delta are timing related
            % - xminimum/xmaximum/yminimum/ymaximum are time/freq values
            b(i)=b(i)+T(1);
            e(i)=b(i)+T(end);
            delta(i)=(e(i)-b(i))/(numel(T)-1);
            npts(i)=numel(P{1});
            nxsize(i)=numel(T);
            nysize(i)=numel(F);
            depmen(i)=mean(data(i).dep(:));
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        else % plotting
            % how this is gonna look
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %        record name        %
            %    +---------------+-+    %
            % amp|  seismogram   |c|    %
            %    +---------------+b|    %
            %  f |               |a| db %
            %  r |  spectrogram  |r|    %
            %  e |               | |    %
            %  q +---------------+-+    %
            %        time (sec)         %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            figure;
            spectrogram(double(data(i).dep),...
                width(i),overlap(i),pow2pad(i),1/delta(i),'yaxis');
            colormap(cmap{i});
            z=zlim;
            zlim(z(1)+diff(z)*hotrod(i,:));
            set(gcf,'color','k');
            set(gca,'color','k','xcolor','w','ycolor','w');
            %set(get(gca,'XLabel'),'color','w');
            %set(get(gca,'YLabel'),'color','w');
            die
        end
        
        % detail message
        if(verbose); print_time_left(i,nrecs); end
    end
    
    % update headers if there is output args
    if(~nargout); return; end
    varargout{1}=changeheader(data,'b',b,'e',e,'delta',delta,...
        'npts',npts,'iftype','ixyz','nxsize',nxsize,'nysize',nysize,...
        'depmen',depmen,'depmin',depmin,'depmax',depmax,...
        'xminimum',b,'xmaximum',e,'yminimum',0,'ymaximum',1./(2*delta));
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
