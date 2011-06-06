function [tt]=ttsurf(model,wave,speed,period,evlalo,stlalo)
%TTSURF    Return surface wave travel time from event(s) to station(s)
%
%    Usage:    tt=ttsurf(model,wave,speed,period,evlalo,stlalo)
%
%    Description:
%     TT=TTSURF(MODEL,WAVE,SPEED,PERIOD,EVLALO,STLALO) computes surface
%     wave travel times between event(s) and station(s) at the given
%     period(s).  MODEL indicates the surface wave dispersion model
%     (default is CUB2), WAVE indicates the wave type (default is
%     'Rayleigh'), & SPEED indicates the wave speed type (default is
%     'group').  PERIOD is the period in seconds at which the travel times
%     are computed.  EVLALO & STLALO are the event and station locations
%     given as Nx2 matrices of [LAT LON].  The output TT is NxM where N is
%     the number of event to station paths and M is the number of periods
%     in PERIOD.
%
%    Notes:
%     - Currently this function is under-developed.  Only CUB2 Rayleigh
%       wave group travel times are implemented.
%
%    Examples:
%     % Compute and contour the CUB2 25s travel time map from 0 Lat, 0 Lon:
%     [lat,lon]=meshgrid(90:-5:-90,-180:5:180);
%     tt=ttsurf('cub2','rayleigh','group',25,[0 0],[lat(:) lon(:)]);
%     figure; contourf(lon',lat',reshape(tt,size(lat))',100:100:7000);
%
%    See also: CUB2, PREM_DISPERSION, AK135_DISPERSION

%     Version History:
%        Jan. 22, 2011 - initial version
%        June  5, 2011 - added prem/ak135 phase velocity, pac group
%                        velocity call
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  5, 2011 at 09:15 GMT

% todo
% - love waves, group and phase, more models ...

% check nargin
error(nargchk(6,6,nargin));

% defaults
if(nargin<1 || isempty(model)); model='cub2'; end
if(nargin<2 || isempty(wave)); wave='rayleigh'; end
if(nargin<3 || isempty(speed)); speed='group'; end

% valid options
validwave={'rayleigh' 'love'};
validspeed={'group' 'phase'};
validmodel={'cub2' 'prem' 'ak135' 'pac'};

% check options
if(~isstring(wave) || ~any(strcmpi(wave,validwave)))
    error('seizmo:ttsurf:badInput',...
        'WAVETYPE must be either ''Rayleigh'' or ''Love''!');
elseif(~isstring(speed) || ~any(strcmpi(speed,validspeed)))
    error('seizmo:ttsurf:badInput',...
        'VELOTYPE must be either ''Group'' or ''Phase''!');
elseif(~isstring(model) || ~any(strcmpi(model,validmodel)))
    error('seizmo:ttsurf:badInput',...
        ['MODEL must be one of the following:\n' ...
        sprintf('''%s'' ',validmodel{:})]);
end

% lowercase options
wave=lower(wave);
speed=lower(speed);
model=lower(model);

% check periods
if(~isreal(period) || any(period<=0))
    error('seizmo:ttsurf:badInput',...
        'PERIOD must be positive reals!');
end
period=period(:);
nper=numel(period);

% check lat/lon points
if(~isreal(stlalo) || size(stlalo,2)~=2)
    error('seizmo:ttsurf:badInput',...
        'STLALO must be a Nx2 array of [LAT LON]!');
elseif(~isreal(evlalo) || size(evlalo,2)~=2)
    error('seizmo:ttsurf:badInput',...
        'EVLALO must be a Nx2 array of [LAT LON]!');
elseif(size(stlalo,1)~=1 && size(evlalo,1)~=1 ...
        && size(stlalo,1)~=size(evlalo,1))
    error('seizmo:ttsurf:badInput',...
        'STLALO & EVLALO must be equal sized or scalar (position)!');
end

% expand lat/lon
if(size(stlalo,1)==1); stlalo=stlalo(ones(size(evlalo,1),1),:); end
if(size(evlalo,1)==1); evlalo=evlalo(ones(size(stlalo,1),1),:); end
npaths=size(stlalo,1);

% preallocate tt
tt=nan(npaths,nper);

% use global for model caching
global SEIZMO

% path step size
step=2; % in degrees

% split by model
switch model
    case 'cub2'
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.CUB2;
        catch
            % not there so load and cache
            model=load('CUB2');
            SEIZMO.MANTLEMODEL.CUB2=model;
        end
        
        % check period
        if(any(period<min(model.period) || period>max(model.period)))
            error('seizmo:ttsurf:badInput',...
                'PERIOD outside model range!');
        end
        
        % split by wave
        switch wave
            case 'rayleigh'
                % split by speed
                switch speed
                    case 'group'
                        % get great circle paths (for all)
                        gcarc=sphericalinv(evlalo(:,1),evlalo(:,2),...
                            stlalo(:,1),stlalo(:,2));
                        npts=max(2,ceil(gcarc/step+1));
                        plat=cell(npaths,1);
                        plon=cell(npaths,1);
                        for i=1:npaths
                            [plat{i},plon{i}]=gcarc2latlon(...
                                evlalo(i,1),evlalo(i,2),...
                                stlalo(i,1),stlalo(i,2),npts(i));
                            % force to 0-360
                            plon{i}(plon{i}<0)=plon{i}(plon{i}<0)+360;
                        end
                        
                        % loop over period
                        for i=1:nper
                            % get map for this period
                            % - get surrounding maps
                            % - simple linear interpolation
                            pi1=find(period(i)>=model.period,1,'last');
                            pi2=find(period(i)<=model.period,1,'first');
                            if(pi1==pi2)
                                map=model.rayleigh_group(:,:,pi1);
                            else
                                d1=(model.period(pi2)-period)/diff(model.period(pi1:pi2));
                                d2=(period-model.period(pi1))/diff(model.period(pi1:pi2));
                                map=d1*model.rayleigh_group(:,:,pi1)...
                                    +d2*model.rayleigh_group(:,:,pi2);
                            end
                            
                            % get travel time(s)
                            % - assuming no gradients
                            for j=1:npaths
                                % get traveltime
                                tt(j,i)=sum(6371*pi/180*gcarc(j)/npts(j)...
                                    ./map(181*round(plon{j})+91-round(plat{j})));
                            end
                        end
                    case 'phase'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                end
            case 'love'
                % split by speed
                switch speed
                    case 'group'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                    case 'phase'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                end
        end
    case 'prem'
        % split by wave
        switch wave
            case 'rayleigh'
                % split by speed
                switch speed
                    case 'group'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                    case 'phase'
                        % get great circle paths (for all)
                        gcarc=sphericalinv(evlalo(:,1),evlalo(:,2),...
                            stlalo(:,1),stlalo(:,2))*6371*pi/180;
                        
                        % get phase velocity at periods
                        phvel=prem_dispersion(1./period);
                        
                        % phase travel times
                        tt=gcarc(:,ones(1,nper))./phvel(:,ones(1,npaths))';
                end
            case 'love'
                % split by speed
                switch speed
                    case 'group'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                    case 'phase'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                end
        end
    case 'ak135'
        % split by wave
        switch wave
            case 'rayleigh'
                % split by speed
                switch speed
                    case 'group'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                    case 'phase'
                        % get great circle paths (for all)
                        gcarc=sphericalinv(evlalo(:,1),evlalo(:,2),...
                            stlalo(:,1),stlalo(:,2))*6371*pi/180;
                        
                        % get phase velocity at periods
                        phvel=prem_dispersion(1./period);
                        
                        % phase travel times
                        tt=gcarc(:,ones(1,nper))./phvel(:,ones(1,npaths))';
                end
            case 'love'
                % split by speed
                switch speed
                    case 'group'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                    case 'phase'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                end
        end
    case 'pac'
        % split by wave
        switch wave
            case 'rayleigh'
                % split by speed
                switch speed
                    case 'group'
                        % get great circle paths (for all)
                        gcarc=sphericalinv(evlalo(:,1),evlalo(:,2),...
                            stlalo(:,1),stlalo(:,2))*6371*pi/180;
                        
                        % get group velocity at periods
                        grpvel=pac_dispersion(1./period);
                        
                        % phase travel times
                        tt=gcarc(:,ones(1,nper))./grpvel(:,ones(1,npaths))';
                    case 'phase'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                end
            case 'love'
                % split by speed
                switch speed
                    case 'group'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                    case 'phase'
                        error('seizmo:ttsurf:notImplemented',...
                            ['Sorry, this combination of MODEL, ' ...
                            'WAVETYPE, & VELOTYPE are not ' ...
                            'implemented yet.']);
                end
        end
end

end
