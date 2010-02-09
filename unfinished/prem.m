function [mout]=prem(varargin)
%PREM    Returns the PREM model
%
%    Usage:    model=prem()
%              model=prem('refperiod',period)
%              model=prem('anisotropic',true|false)
%              model=prem('ocean',true|false)
%              model=prem('depths',depths)
%              model=prem('spvw',nsamples)
%              model=prem('maxdepth',maxdepth)
%              model=prem('ddptype',type)
%              model=prem('ddpthickness',thickness)
%              model=prem('ddpvp',poly)
%              model=prem('ddpvs',poly)
%              model=prem('ulvzthickness',thickness)
%              model=prem('ulvzvelocity',velocity)
%
%    Description:
%
%    Notes:
%     - MODEL: depth rho vpv vsv qk qu vph vsh eta
%              depth rho vp vs qk qu
%
%    Examples:
%
%    See also: iasp91, ak135

%     Version History:
%        Jan. 17, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 17, 2010 at 07:40 GMT

% todo:
% - mass polynomials
%
%               re
%              f
%     m(r)  =  | p(r)*4*pi*r^2*dr
%              j
%             0
%
% - gravity polynomials
%
%                           re
%                          f
%     g(r)  =  Gm(r)/r  =  | G*p(r)*4*pi*r*dr
%                          j
%                         0
%
% - pressure polynomials
%
%               0
%              f
%     P(r)  =  | p(r)*g(r)*dr
%              j
%            re
%
% - dv/dr (just take dirivative)

% id funtion
me=mfilename;
pkgme=['seizmo:' me ':'];

% all inputs must be 'option',value pairs
if(mod(nargin,2))
    error([pkgme 'unpairedOption'],...
        'Option missing a value!');
end

% get global
global SEIZMO

% define some constants
A=pi*4e12;
G=6.67428e-17;
Re=6371;
Re3=Re^3;

% valid values for string options
valid.DDPTYPE={'default'};

% option defaults
option.REFPERIOD=1; % reference period of 1sec
option.ANISOTROPY=true; % anisotropy
option.OCEAN=true; % ocean
option.DEPTHS=[]; % no predefined depths
option.SPVW=1; % minimum of 1 sample(s) per vertical wavelength
option.MAXDEPTH=[]; % no max depth
option.DDPTYPE='default'; % d" type ('default' is normal)
option.DDPTHICKNESS=150; % d" thickness (km)
option.DDPVP=0; % d" vp % gradient adjustment
option.DDPVS=0; % d" vs % gradient adjustment
option.ULVZTHICKNESS=0; % ulvz thickness (km)
option.ULVZVELOCITY=-10; % ulvz % velocity change

% get options from SEIZMO global
ME=upper(me);
try
    fields=fieldnames(SEIZMO.(ME));
    for i=1:numel(fields)
        if(~isempty(SEIZMO.(ME).(fields{i})))
            option.(fields{i})=SEIZMO.(ME).(fields{i});
        end
    end
catch
end

% get options from command line
for i=1:2:nargin-1
    if(~ischar(varargin{i}))
        error([pkgme 'badInput'],...
            'Options must be specified as a string!');
    end
    if(~isempty(varargin{i+1}))
        option.(upper(varargin{i}))=varargin{i+1};
    end
end

% check options


% model details
% is it anisotropic?
if(option.ANISOTROPY)
    % do we have an ocean?
    if(option.OCEAN)
        % the ocean (3km)
        model.order={'ocean'};
        model.ocean.radius=[6368 6371];
        model.ocean.depth=[0 3];
        model.ocean.rho=1.02;
        model.ocean.vpv=1.45;
        model.ocean.vsv=0;
        model.ocean.qk=57823;
        model.ocean.qu=99999.9;
        model.ocean.vph=1.45;
        model.ocean.vsh=0;
        model.ocean.eta=1;
        
        % the upper crust
        model.order=[model.order {'uppercrust'}];
        model.uppercrust.radius=[6356 6368];
        model.uppercrust.depth=[3 15];
        model.uppercrust.rho=2.6;
        model.uppercrust.vpv=5.8;
        model.uppercrust.vsv=3.2;
        model.uppercrust.qk=57823;
        model.uppercrust.qu=600;
        model.uppercrust.vph=5.8;
        model.uppercrust.vsh=3.2;
        model.uppercrust.eta=1;
    else % NO OCEAN
        % the upper crust (to the surface)
        model.order={'uppercrust'};
        model.uppercrust.radius=[6356 6371];
        model.uppercrust.depth=[0 15];
        model.uppercrust.rho=2.6;
        model.uppercrust.vpv=5.8;
        model.uppercrust.vsv=3.2;
        model.uppercrust.qk=57823;
        model.uppercrust.qu=600;
        model.uppercrust.vph=5.8;
        model.uppercrust.vsh=3.2;
        model.uppercrust.eta=1;
    end
    
    % the lower crust
    model.order=[model.order {'lowercrust'}];
    model.lowercrust.radius=[6346.6 6356];
    model.lowercrust.depth=[15 24.4];
    model.lowercrust.rho=2.9;
    model.lowercrust.vpv=6.8;
    model.lowercrust.vsv=3.9;
    model.lowercrust.qk=57823;
    model.lowercrust.qu=600;
    model.lowercrust.vph=6.8;
    model.lowercrust.vsh=3.9;
    model.lowercrust.eta=1;
    
    % the mantle lithosphere (lid)
    model.order=[model.order {'lithosphere'}];
    model.lithosphere.radius=[6291 6346.6];
    model.lithosphere.depth=[24.4 80];
    model.lithosphere.rho=[2.6910 0.6924];
    model.lithosphere.vpv=[0.8317 7.2180];
    model.lithosphere.vsv=[5.8582 -1.4678];
    model.lithosphere.qk=57823;
    model.lithosphere.qu=600;
    model.lithosphere.vph=[3.5908 4.6172];
    model.lithosphere.vsh=[-1.0839 5.7176];
    model.lithosphere.eta=[3.3687 -2.4778];
    
    % the mantle asthenosphere (lvz)
    model.order=[model.order {'asthenosphere'}];
    model.asthenosphere.radius=[6151 6291];
    model.asthenosphere.depth=[80 220];
    model.asthenosphere.rho=[2.6910 0.6924];
    model.asthenosphere.vpv=[0.8317 7.2180];
    model.asthenosphere.vsv=[5.8582 -1.4678];
    model.asthenosphere.qk=57823;
    model.asthenosphere.qu=80;
    model.asthenosphere.vph=[3.5908 4.6172];
    model.asthenosphere.vsh=[-1.0839 5.7176];
    model.asthenosphere.eta=[3.3687 -2.4778];
    
    % the upper mantle
    model.order=[model.order {'uppermantle'}];
    model.uppermantle.radius=[5971 6151];
    model.uppermantle.depth=[220 400];
    model.uppermantle.rho=[7.1089 -3.8045];
    model.uppermantle.vpv=[20.3926 -12.2569];
    model.uppermantle.vsv=[8.9496 -4.4597];
    model.uppermantle.qk=57823;
    model.uppermantle.qu=143;
    model.uppermantle.vph=[20.3926 -12.2569];
    model.uppermantle.vsh=[8.9496 -4.4597];
    model.uppermantle.eta=1;
    
    % upper transition zone
    model.order=[model.order {'uppertz'}];
    model.uppertz.radius=[5771 5971];
    model.uppertz.depth=[400 600];
    model.uppertz.rho=[11.2494 -8.0298];
    model.uppertz.vpv=[39.7027 -32.6166];
    model.uppertz.vsv=[22.3512 -18.5856];
    model.uppertz.qk=57823;
    model.uppertz.qu=143;
    model.uppertz.vph=[39.7027 -32.6166];
    model.uppertz.vsh=[22.3512 -18.5856];
    model.uppertz.eta=1;
    
    % middle transition zone
    model.order=[model.order {'middletz'}];
    model.middletz.radius=[5701 5771];
    model.middletz.depth=[600 670];
    model.middletz.rho=[5.3197 -1.4836];
    model.middletz.vpv=[19.0957 -9.8672];
    model.middletz.vsv=[9.9839 -4.9324];
    model.middletz.qk=57823;
    model.middletz.qu=143;
    model.middletz.vph=[19.0957 -9.8672];
    model.middletz.vsh=[9.9839 -4.9324];
    model.middletz.eta=1;
    
    % lower transition zone
    model.order=[model.order {'lowertz'}];
    model.lowertz.radius=[5600 5701];
    model.lowertz.depth=[670 771];
    model.lowertz.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.lowertz.vpv=[29.2766 -23.6027 5.5242 -2.5514];
    model.lowertz.vsv=[22.3459 -17.2473 -2.0834 0.9783];
    model.lowertz.qk=57823;
    model.lowertz.qu=312;
    model.lowertz.vph=[29.2766 -23.6027 5.5242 -2.5514];
    model.lowertz.vsh=[22.3459 -17.2473 -2.0834 0.9783];
    model.lowertz.eta=1;
    
    % lower mantle
    model.order=[model.order {'lowermantle'}];
    %model.lowermantle.radius=[3630 5600];
    %model.lowermantle.depth=[771 2741];
    model.lowermantle.radius=[3480+option.DDPTHICKNESS 5600];
    model.lowermantle.depth=[771 2891-option.DDPTHICKNESS];
    model.lowermantle.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.lowermantle.vpv=[24.9520 -40.4673 51.4832 -26.6419];
    model.lowermantle.vsv=[11.1671 -13.7818 17.4575 -9.2777];
    model.lowermantle.qk=57823;
    model.lowermantle.qu=312;
    model.lowermantle.vph=[24.9520 -40.4673 51.4832 -26.6419];
    model.lowermantle.vsh=[11.1671 -13.7818 17.4575 -9.2777];
    model.lowermantle.eta=1;
    
    % normal or custom d"?
    switch option.DDPTYPE
        case 'default'
            % d"
            model.order=[model.order {'ddp'}];
            %model.ddp.radius=[3480 3630];
            %model.ddp.depth=[2741 2891];
            model.ddp.radius=[3480 3480+option.DDPTHICKNESS];
            model.ddp.depth=[2891-option.DDPTHICKNESS 2891];
            model.ddp.rho=[7.9565 -6.4761 5.5283 -3.0807];
            model.ddp.vpv=[15.3891 -5.3181 5.5242 -2.5514];
            model.ddp.vsv=[6.9254 1.4672 -2.0834 0.9783];
            model.ddp.qk=57823;
            model.ddp.qu=312;
            model.ddp.vph=[15.3891 -5.3181 5.5242 -2.5514];
            model.ddp.vsh=[6.9254 1.4672 -2.0834 0.9783];
            model.ddp.eta=1;

            % adjust velocity gradient
            model.ddp.vpv(2:end)=model.ddp.vpv(2:end)*(100+option.DDPVP)/100;
            model.ddp.vsv(2:end)=model.ddp.vsv(2:end)*(100+option.DDPVS)/100;
            model.ddp.vph(2:end)=model.ddp.vph(2:end)*(100+option.DDPVP)/100;
            model.ddp.vsh(2:end)=model.ddp.vsh(2:end)*(100+option.DDPVS)/100;

            % correct for changes
            r=model.ddp.radius(2)/6371;
            model.ddp.vpv(1)=model.ddp.vpv(1) ...
                +polyval(fliplr(model.lowermantle.vpv),r) ...
                -polyval(fliplr(model.ddp.vpv),r);
            model.ddp.vsv(1)=model.ddp.vsv(1) ...
                +polyval(fliplr(model.lowermantle.vsv),r) ...
                -polyval(fliplr(model.ddp.vsv),r);
            model.ddp.vph(1)=model.ddp.vph(1) ...
                +polyval(fliplr(model.lowermantle.vph),r) ...
                -polyval(fliplr(model.ddp.vph),r);
            model.ddp.vsh(1)=model.ddp.vsh(1) ...
                +polyval(fliplr(model.lowermantle.vsh),r) ...
                -polyval(fliplr(model.ddp.vsh),r);
        case 'custom'
            % how do we approach this?
            % thickness field function value_top value_bottom option
            %           field function value_top value_bottom option
            %           field function value_top value_bottom option
            % function: 1 - power
            %           2 - error function
            %           3 - comp. error function
            % - how do supply parameters affecting the shape of erf?
            %   - max x of erf(x) (def to 2)
            % - how about something like what i did for adv geodynamics?
            %
            % ddp.layer(1).thickness=50;
            %             .vpv.function='power';
            %                 .top=15;
            %                 .bottom=16;
            %                 .option=1;
            %
            % fields that can be edited:
            %  vpv, vph, vsv, vsh, qu, qk, eta
            %
            % fields that cannot be edited:
            %  rho (need to do lots of integration stuff for erf/erfc)
            %
            % antiderivative of erf(z):
            %  z*erf(z) + e^(-z^2)/sqrt(pi) - 1/sqrt(pi)
            %
            % antiderivative of erfc(z):
            %  z - z*erf(z) - e^(-z^2)/sqrt(pi) + 1/sqrt(pi)
            %
            % ok what to do:
            % - change thickness to depths/radii
            % - functions to handle custom types
    end
    
    % outer core
    model.order=[model.order {'outercore'}];
    model.outercore.radius=[1221.5 3480];
    model.outercore.depth=[2891 5149.5];
    model.outercore.rho=[12.5815 -1.2638 -3.6526 -5.5281];
    model.outercore.vpv=[11.0487 -4.0362 4.8023 -13.5732];
    model.outercore.vsv=0;
    model.outercore.qk=57823;
    model.outercore.qu=99999.9;
    model.outercore.vph=[11.0487 -4.0362 4.8023 -13.5732];
    model.outercore.vsh=0;
    model.outercore.eta=1;
    
    % inner core
    model.order=[model.order {'innercore'}];
    model.innercore.radius=[0 1221.5];
    model.innercore.depth=[5149.5 6371];
    model.innercore.rho=[13.0885 0 -8.8381];
    model.innercore.vpv=[11.2622 0 -6.3640];
    model.innercore.vsv=[3.6678 0 -4.4475];
    model.innercore.qk=1327.7;
    model.innercore.qu=84.6;
    model.innercore.vph=[11.2622 0 -6.3640];
    model.innercore.vsh=[3.6678 0 -4.4475];
    model.innercore.eta=1;
    
    % get mass polynomials
    no=numel(model.order);
    for i=1:no
        nm=numel(model.(model.order{i}).rho);
        model.(model.order{i}).m=...
            A*Re3*[0 0 0 model.(model.order{i}).rho./(3:2+nm)];
    end
    
    % correct mass polynomials
    for i=no-1:-1:1
        model.(model.order{i}).m(1)=...
            polyval(fliplr(model.(model.order{i+1}).m),...
            model.(model.order{i+1}).radius(2)/Re)...
            -polyval(fliplr(model.(model.order{i}).m),...
            model.(model.order{i}).radius(1)/Re);
    end
    
    % gravity polynomials start at r^-2 (ie must be divided by r^2)
    for i=1:no
        model.(model.order{i}).g=G/Re^2*model.(model.order{i}).m;
    end
    
    % pressure polynomials start at r^-2 (ie must be divided by r^2)
    for i=1:no
        model.(model.order{i}).p=...
            conv(model.(model.order{i}).g,model.(model.order{i}).rho);
        model.(model.order{i}).p
    end
    
    % output structure
    fields={'depth' 'rho' 'vpv' 'vsv' 'qk' 'qu' ...
        'vph' 'vsh' 'eta' 'm' 'g'};
    nf=numel(fields);
    vpv=strcmp(fields,'vpv'); vsv=strcmp(fields,'vsv');
    vph=strcmp(fields,'vph'); vsh=strcmp(fields,'vsh');
    qk=strcmp(fields,'qk');   qu=strcmp(fields,'qu');
    g=strcmp(fields,'g');
    
    % specific depths
    if(~isempty(option.DEPTHS))
        % preallocate output
        mout=nan(numel(option.DEPTHS),nf);
        
        % loop over depths
        cnt=0;
        for i=1:numel(option.DEPTHS)
            % loop over regions
            for j=1:no
                % is depth within this region? (or at the boundary)
                if(option.DEPTHS(i)>=model.(model.order{j}).depth(1) ...
                        && option.DEPTHS(i)...
                                        <=model.(model.order{j}).depth(2))
                    % yes, so increment
                    cnt=cnt+1;
                    
                    % allow for custom
                    if(strncmpi(model.order{j},'ddp',3) ...
                            && strcmpi(option.DDPTYPE,'custom'))
                        % custom method
                        
                    else
                        % default polynomial method
                        r=(Re-option.DEPTHS(i))/Re;
                        mout(cnt,1)=option.DEPTHS(i);
                        for k=2:nf
                            mout(cnt,k)=polyval(fliplr(...
                                model.(model.order{j}).(fields{k})),r);
                        end
                        if(r>0); mout(cnt,g)=mout(cnt,g)/r^2; end
                    end
                    
                    % correct velocities to be at REFPERIOD
                    % L is based on isotropic velocities for Voigt solid
                    vsvoigt=sqrt((2*mout(cnt,vsv)^2+mout(cnt,vsh)^2)/3);
                    vpvoigt=sqrt((mout(cnt,vpv)^2+4*mout(cnt,vph)^2)/5);
                    L=4/3*(vsvoigt/vpvoigt)^2;
                    Ap=1-log(option.REFPERIOD)/pi*...
                        (L/mout(cnt,qu)+(1-L)/mout(cnt,qk));
                    As=1-log(option.REFPERIOD)/pi/mout(cnt,qu);
                    mout(cnt,vpv)=mout(cnt,vpv)*Ap;
                    mout(cnt,vph)=mout(cnt,vph)*Ap;
                    mout(cnt,vsv)=mout(cnt,vsv)*As;
                    mout(cnt,vsh)=mout(cnt,vsh)*As;
                    
                    % exit loop (avoids extra output at discontinuities)
                    %break;
                end
            end
        end
    else % SPVW at REFPERIOD
        % setup / preallocate output
        mout=nan(1100,nf); % 1100 nodes is a guess
        n=1;
        
        % loop over regions
        for i=1:no
            % start at top of region
            depth=model.(model.order{i}).depth(1);
            r=(Re-depth)/Re;
            mout(n,1)=depth;
            for k=2:nf
                mout(n,k)=polyval(...
                    fliplr(model.(model.order{i}).(fields{k})),r);
            end
            if(r>0); mout(n,g)=mout(n,g)/r^2; end
            
            % step to next depth
            if(mout(n,vsv)>0)
                vmin=mout(n,vsv);
            else
                vmin=mout(n,vpv);
            end
            wavelength=vmin*option.REFPERIOD;
            stepsize=wavelength/option.SPVW;
            depth=depth+stepsize;
            n=n+1;
            
            % finished if below max depth
            if(depth>option.MAXDEPTH); break; end
            
            % step until at/below bottom of region (within a meter)
            while(depth<(model.(model.order{i}).depth(2)-0.001))
                r=(Re-depth)/Re;
                mout(n,1)=depth;
                for k=2:nf
                    mout(n,k)=polyval(...
                        fliplr(model.(model.order{i}).(fields{k})),r);
                end
                if(r>0); mout(n,g)=mout(n,g)/r^2; end

                % step to next depth
                if(mout(n,vsv)>0)
                    vmin=mout(n,vsv);
                else
                    vmin=mout(n,vpv);
                end
                wavelength=vmin*option.REFPERIOD;
                stepsize=wavelength/option.SPVW;
                depth=depth+stepsize;
                n=n+1;
                
                % finished if below max depth
                if(depth>option.MAXDEPTH); break; end
            end
            
            % finish at the bottom of region
            depth=model.(model.order{i}).depth(2);
            if(depth>option.MAXDEPTH); break; end
            r=(Re-depth)/Re;
            mout(n,1)=depth;
            for k=2:nf
                mout(n,k)=polyval(...
                    fliplr(model.(model.order{i}).(fields{k})),r);
            end
            if(r>0); mout(n,g)=mout(n,g)/r^2; end
            n=n+1;
        end
        
        % truncate unused mout
        n=n-1;
        mout=mout(1:n,:);
    end
else % ISOTROPIC
    % do we have an ocean?
    if(option.OCEAN)
        % the ocean (3km)
        model.order={'ocean'};
        model.ocean.radius=[6368 6371];
        model.ocean.depth=[0 3];
        model.ocean.rho=[1.02 0 0 0 0];
        model.ocean.vp=[1.45 0 0 0 0];
        model.ocean.vs=[0 0 0 0 0];
        model.ocean.qk=[57823 0 0 0 0];
        model.ocean.qu=[99999.9 0 0 0 0];
        
        % the upper crust
        model.order=[model.order {'uppercrust'}];
        model.uppercrust.radius=[6356 6368];
        model.uppercrust.depth=[3 15];
        model.uppercrust.rho=[2.6 0 0 0 0];
        model.uppercrust.vp=[5.8 0 0 0 0];
        model.uppercrust.vs=[3.2 0 0 0 0];
        model.uppercrust.qk=[57823 0 0 0 0];
        model.uppercrust.qu=[600 0 0 0 0];
    else % NO OCEAN
        % the upper crust (to the surface)
        model.order={'uppercrust'};
        model.uppercrust.radius=[6356 6371];
        model.uppercrust.depth=[0 15];
        model.uppercrust.rho=[2.6 0 0 0 0];
        model.uppercrust.vp=[5.8 0 0 0 0];
        model.uppercrust.vs=[3.2 0 0 0 0];
        model.uppercrust.qk=[57823 0 0 0 0];
        model.uppercrust.qu=[600 0 0 0 0];
    end
    
    % the lower crust
    model.order=[model.order {'lowercrust'}];
    model.lowercrust.radius=[6346.6 6356];
    model.lowercrust.depth=[15 24.4];
    model.lowercrust.rho=[2.9 0 0 0 0];
    model.lowercrust.vp=[6.8 0 0 0 0];
    model.lowercrust.vs=[3.9 0 0 0 0];
    model.lowercrust.qk=[57823 0 0 0 0];
    model.lowercrust.qu=[600 0 0 0 0];
    
    % the mantle lithosphere (lid)
    model.order=[model.order {'lithosphere'}];
    model.lithosphere.radius=[6291 6346.6];
    model.lithosphere.depth=[24.4 80];
    model.lithosphere.rho=[2.6910 0.6924 0 0 0];
    model.lithosphere.vpv=[4.1875 3.9382 0 0 0];
    model.lithosphere.vsv=[2.1519 2.3481 0 0 0];
    model.lithosphere.qk=[57823 0 0 0 0];
    model.lithosphere.qu=[600 0 0 0 0];
    
    % the mantle asthenosphere (lvz)
    model.order=[model.order {'asthenosphere'}];
    model.asthenosphere.radius=[6151 6291];
    model.asthenosphere.depth=[80 220];
    model.asthenosphere.rho=[2.6910 0.6924 0 0 0];
    model.asthenosphere.vpv=[4.1875 3.9382 0 0 0];
    model.asthenosphere.vsv=[2.1519 2.3481 0 0 0];
    model.asthenosphere.qk=[57823 0 0 0 0];
    model.asthenosphere.qu=[80 0 0 0 0];
    
    % the upper mantle
    model.order=[model.order {'uppermantle'}];
    model.uppermantle.radius=[5971 6151];
    model.uppermantle.depth=[220 400];
    model.uppermantle.rho=[7.1089 -3.8045 0 0 0];
    model.uppermantle.vp=[20.3926 -12.2569 0 0 0];
    model.uppermantle.vs=[8.9496 -4.4597 0 0 0];
    model.uppermantle.qk=[57823 0 0 0 0];
    model.uppermantle.qu=[143 0 0 0 0];
    
    % upper transition zone
    model.order=[model.order {'uppertz'}];
    model.uppertz.radius=[5771 5971];
    model.uppertz.depth=[400 600];
    model.uppertz.rho=[11.2494 -8.0298 0 0 0];
    model.uppertz.vp=[39.7027 -32.6166 0 0 0];
    model.uppertz.vs=[22.3512 -18.5856 0 0 0];
    model.uppertz.qk=[57823 0 0 0 0];
    model.uppertz.qu=[143 0 0 0 0];
    
    % middle transition zone
    model.order=[model.order {'middletz'}];
    model.middletz.radius=[5701 5771];
    model.middletz.depth=[600 670];
    model.middletz.rho=[5.3197 -1.4836 0 0 0];
    model.middletz.vp=[19.0957 -9.8672 0 0 0];
    model.middletz.vs=[9.9839 -4.9324 0 0 0];
    model.middletz.qk=[57823 0 0 0 0];
    model.middletz.qu=[143 0 0 0 0];
    
    % lower transition zone
    model.order=[model.order {'lowertz'}];
    model.lowertz.radius=[5600 5701];
    model.lowertz.depth=[670 771];
    model.lowertz.rho=[7.9565 -6.4761 5.5283 -3.0807 0];
    model.lowertz.vp=[29.2766 -23.6027 5.5242 -2.5514 0];
    model.lowertz.vs=[22.3459 -17.2473 -2.0834 0.9783 0];
    model.lowertz.qk=[57823 0 0 0 0];
    model.lowertz.qu=[312 0 0 0 0];
    
    % lower mantle
    model.order=[model.order {'lowermantle'}];
    model.lowermantle.radius=[3630 5600];
    model.lowermantle.depth=[771 2741];
    model.lowermantle.rho=[7.9565 -6.4761 5.5283 -3.0807 0];
    model.lowermantle.vp=[24.9520 -40.4673 51.4832 -26.6419 0];
    model.lowermantle.vs=[11.1671 -13.7818 17.4575 -9.2777 0];
    model.lowermantle.qk=[57823 0 0 0 0];
    model.lowermantle.qu=[312 0 0 0 0];
    
    % d"
    model.order=[model.order {'ddp'}];
    model.ddp.radius=[3480 3630];
    model.ddp.depth=[2741 2891];
    model.ddp.rho=[7.9565 -6.4761 5.5283 -3.0807 0];
    model.ddp.vp=[15.3891 -5.3181 5.5242 -2.5514 0];
    model.ddp.vs=[6.9254 1.4672 -2.0834 0.9783 0];
    model.ddp.qk=[57823 0 0 0 0];
    model.ddp.qu=[312 0 0 0 0];
    
    % outer core
    model.order=[model.order {'outercore'}];
    model.outercore.radius=[1221.5 3480];
    model.outercore.depth=[2891 5149.5];
    model.outercore.rho=[12.5815 -1.2638 -3.6526 -5.5281 0];
    model.outercore.vp=[11.0487 -4.0362 4.8023 -13.5732 0];
    model.outercore.vs=[0 0 0 0 0];
    model.outercore.qk=[57823 0 0 0 0];
    model.outercore.qu=[99999.9 0 0 0 0];
    
    % inner core
    model.order=[model.order {'innercore'}];
    model.innercore.radius=[0 1221.5];
    model.innercore.depth=[5149.5 6371];
    model.innercore.rho=[13.0885 0 -8.8381 0 0];
    model.innercore.vp=[11.2622 0 -6.3640 0 0];
    model.innercore.vs=[3.6678 0 -4.4475 0 0];
    model.innercore.qk=[1327.7 0 0 0 0];
    model.innercore.qu=[84.6 0 0 0 0];
    
    
    % output structure
    fields={'depth' 'rho' 'vp' 'vs' 'qk' 'qu'};
    vp=strcmp(fields,'vp'); vs=strcmp(fields,'vs');
    qk=strcmp(fields,'qk'); qu=strcmp(fields,'qu');
    
    % specific depths
    if(~isempty(option.DEPTHS))
        % preallocate output
        mout=nan(numel(option.DEPTHS),6);
        
        % loop over depths
        for i=1:numel(option.DEPTHS)
            % loop over regions
            for j=1:numel(model.order)
                % is depth within this region? (or at the boundary)
                if(option.DEPTHS(i)>=model.(model.order{j}).depth(1) ...
                        && option.DEPTHS(i)...
                                        <=model.(model.order{j}).depth(2))
                    % yes, so now evaluate
                    r=(6371-option.DEPTHS(i))/6371;
                    mout(i,1)=option.DEPTHS(i);
                    for k=2:numel(fields)
                        mout(i,k)=...
                            poly2val(r,model.(model.order{j}).(fields{k}));
                    end
                    
                    % correct velocities to be at REFPERIOD
                    % L is based on isotropic velocities for Voigt solid
                    vsvoigt=mout(i,vs);
                    vpvoigt=mout(i,vp);
                    L=4/3*(vsvoigt/vpvoigt)^2;
                    Ap=1-log(option.REFPERIOD)/pi*...
                        (L/mout(i,qu)+(1-L)/mout(i,qk));
                    As=1-log(option.REFPERIOD)/pi/mout(i,qu);
                    mout(i,vp)=mout(i,vp)*Ap;
                    mout(i,vs)=mout(i,vs)*As;
                    
                    % exit loop (avoids extra output at discontinuities)
                    break;
                end
            end
        end
    else % SPVW at REFPERIOD
        % setup / preallocate output
        mout=nan(1100,6); % 1100 nodes is a guess
        n=1;
        
        % loop over regions
        for i=1:numel(model.order)
            % start at top of region
            depth=model.(model.order{i}).depth(1);
            r=(6371-depth)/6371;
            mout(n,1)=depth;
            for k=2:numel(fields)
                mout(n,k)=poly2val(r,model.(model.order{i}).(fields{k}));
            end
            
            % step to next depth
            if(mout(n,vs)>0)
                vmin=mout(n,vs);
            else
                vmin=mout(n,vp);
            end
            wavelength=vmin*option.REFPERIOD;
            stepsize=wavelength/option.SPVW;
            depth=depth+stepsize;
            n=n+1;
            
            % finished if below max depth
            if(depth>option.MAXDEPTH); break; end
            
            % step until at/below bottom of region (within a meter)
            while(depth<(model.(model.order{i}).depth(2)-0.001))
                r=(6371-depth)/6371;
                mout(n,1)=depth;
                for k=2:numel(fields)
                    mout(n,k)=...
                        poly2val(r,model.(model.order{i}).(fields{k}));
                end

                % step to next depth
                if(mout(n,vs)>0)
                    vmin=mout(n,vs);
                else
                    vmin=mout(n,vp);
                end
                wavelength=vmin*option.REFPERIOD;
                stepsize=wavelength/option.SPVW;
                depth=depth+stepsize;
                n=n+1;
                
                % finished if below max depth
                if(depth>option.MAXDEPTH); break; end
            end
            
            % finish at the bottom of region
            depth=model.(model.order{i}).depth(2);
            if(depth>option.MAXDEPTH); break; end
            r=(6371-depth)/6371;
            mout(n,1)=depth;
            for k=2:numel(fields)
                mout(n,k)=poly2val(r,model.(model.order{i}).(fields{k}));
            end
            n=n+1;
        end
        
        % truncate unused mout
        n=n-1;
        mout=mout(1:n,:);
    end
end

end
