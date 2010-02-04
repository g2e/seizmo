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
%   - m=SpdV
% - gravity polynomials
%   - g=Gm/r
% - pressure polynomials
%   - 
% - dv/dr

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
        model.ocean.rho=[1.02 0 0 0 0];
        model.ocean.vpv=[1.45 0 0 0 0];
        model.ocean.vsv=[0 0 0 0 0];
        model.ocean.qk=[57823 0 0 0 0];
        model.ocean.qu=[99999.9 0 0 0 0];
        model.ocean.vph=[1.45 0 0 0 0];
        model.ocean.vsh=[0 0 0 0 0];
        model.ocean.eta=[1 0 0 0 0];
        
        % the upper crust
        model.order=[model.order {'uppercrust'}];
        model.uppercrust.radius=[6356 6368];
        model.uppercrust.depth=[3 15];
        model.uppercrust.rho=[2.6 0 0 0 0];
        model.uppercrust.vpv=[5.8 0 0 0 0];
        model.uppercrust.vsv=[3.2 0 0 0 0];
        model.uppercrust.qk=[57823 0 0 0 0];
        model.uppercrust.qu=[600 0 0 0 0];
        model.uppercrust.vph=[5.8 0 0 0 0];
        model.uppercrust.vsh=[3.2 0 0 0 0];
        model.uppercrust.eta=[1 0 0 0 0];
    else % NO OCEAN
        % the upper crust (to the surface)
        model.order={'uppercrust'};
        model.uppercrust.radius=[6356 6371];
        model.uppercrust.depth=[0 15];
        model.uppercrust.rho=[2.6 0 0 0 0];
        model.uppercrust.vpv=[5.8 0 0 0 0];
        model.uppercrust.vsv=[3.2 0 0 0 0];
        model.uppercrust.qk=[57823 0 0 0 0];
        model.uppercrust.qu=[600 0 0 0 0];
        model.uppercrust.vph=[5.8 0 0 0 0];
        model.uppercrust.vsh=[3.2 0 0 0 0];
        model.uppercrust.eta=[1 0 0 0 0];
    end
    
    % the lower crust
    model.order=[model.order {'lowercrust'}];
    model.lowercrust.radius=[6346.6 6356];
    model.lowercrust.depth=[15 24.4];
    model.lowercrust.rho=[2.9 0 0 0 0];
    model.lowercrust.vpv=[6.8 0 0 0 0];
    model.lowercrust.vsv=[3.9 0 0 0 0];
    model.lowercrust.qk=[57823 0 0 0 0];
    model.lowercrust.qu=[600 0 0 0 0];
    model.lowercrust.vph=[6.8 0 0 0 0];
    model.lowercrust.vsh=[3.9 0 0 0 0];
    model.lowercrust.eta=[1 0 0 0 0];
    
    % the mantle lithosphere (lid)
    model.order=[model.order {'lithosphere'}];
    model.lithosphere.radius=[6291 6346.6];
    model.lithosphere.depth=[24.4 80];
    model.lithosphere.rho=[2.6910 0.6924 0 0 0];
    model.lithosphere.vpv=[0.8317 7.2180 0 0 0];
    model.lithosphere.vsv=[5.8582 -1.4678 0 0 0];
    model.lithosphere.qk=[57823 0 0 0 0];
    model.lithosphere.qu=[600 0 0 0 0];
    model.lithosphere.vph=[3.5908 4.6172 0 0 0];
    model.lithosphere.vsh=[-1.0839 5.7176 0 0 0];
    model.lithosphere.eta=[3.3687 -2.4778 0 0 0];
    
    % the mantle asthenosphere (lvz)
    model.order=[model.order {'asthenosphere'}];
    model.asthenosphere.radius=[6151 6291];
    model.asthenosphere.depth=[80 220];
    model.asthenosphere.rho=[2.6910 0.6924 0 0 0];
    model.asthenosphere.vpv=[0.8317 7.2180 0 0 0];
    model.asthenosphere.vsv=[5.8582 -1.4678 0 0 0];
    model.asthenosphere.qk=[57823 0 0 0 0];
    model.asthenosphere.qu=[80 0 0 0 0];
    model.asthenosphere.vph=[3.5908 4.6172 0 0 0];
    model.asthenosphere.vsh=[-1.0839 5.7176 0 0 0];
    model.asthenosphere.eta=[3.3687 -2.4778 0 0 0];
    
    % the upper mantle
    model.order=[model.order {'uppermantle'}];
    model.uppermantle.radius=[5971 6151];
    model.uppermantle.depth=[220 400];
    model.uppermantle.rho=[7.1089 -3.8045 0 0 0];
    model.uppermantle.vpv=[20.3926 -12.2569 0 0 0];
    model.uppermantle.vsv=[8.9496 -4.4597 0 0 0];
    model.uppermantle.qk=[57823 0 0 0 0];
    model.uppermantle.qu=[143 0 0 0 0];
    model.uppermantle.vph=[20.3926 -12.2569 0 0 0];
    model.uppermantle.vsh=[8.9496 -4.4597 0 0 0];
    model.uppermantle.eta=[1 0 0 0 0];
    
    % upper transition zone
    model.order=[model.order {'uppertz'}];
    model.uppertz.radius=[5771 5971];
    model.uppertz.depth=[400 600];
    model.uppertz.rho=[11.2494 -8.0298 0 0 0];
    model.uppertz.vpv=[39.7027 -32.6166 0 0 0];
    model.uppertz.vsv=[22.3512 -18.5856 0 0 0];
    model.uppertz.qk=[57823 0 0 0 0];
    model.uppertz.qu=[143 0 0 0 0];
    model.uppertz.vph=[39.7027 -32.6166 0 0 0];
    model.uppertz.vsh=[22.3512 -18.5856 0 0 0];
    model.uppertz.eta=[1 0 0 0 0];
    
    % middle transition zone
    model.order=[model.order {'middletz'}];
    model.middletz.radius=[5701 5771];
    model.middletz.depth=[600 670];
    model.middletz.rho=[5.3197 -1.4836 0 0 0];
    model.middletz.vpv=[19.0957 -9.8672 0 0 0];
    model.middletz.vsv=[9.9839 -4.9324 0 0 0];
    model.middletz.qk=[57823 0 0 0 0];
    model.middletz.qu=[143 0 0 0 0];
    model.middletz.vph=[19.0957 -9.8672 0 0 0];
    model.middletz.vsh=[9.9839 -4.9324 0 0 0];
    model.middletz.eta=[1 0 0 0 0];
    
    % lower transition zone
    model.order=[model.order {'lowertz'}];
    model.lowertz.radius=[5600 5701];
    model.lowertz.depth=[670 771];
    model.lowertz.rho=[7.9565 -6.4761 5.5283 -3.0807 0];
    model.lowertz.vpv=[29.2766 -23.6027 5.5242 -2.5514 0];
    model.lowertz.vsv=[22.3459 -17.2473 -2.0834 0.9783 0];
    model.lowertz.qk=[57823 0 0 0 0];
    model.lowertz.qu=[312 0 0 0 0];
    model.lowertz.vph=[29.2766 -23.6027 5.5242 -2.5514 0];
    model.lowertz.vsh=[22.3459 -17.2473 -2.0834 0.9783 0];
    model.lowertz.eta=[1 0 0 0 0];
    
    % lower mantle
    model.order=[model.order {'lowermantle'}];
    %model.lowermantle.radius=[3630 5600];
    %model.lowermantle.depth=[771 2741];
    model.lowermantle.radius=[3480+option.DDPTHICKNESS 5600];
    model.lowermantle.depth=[771 2891-option.DDPTHICKNESS];
    model.lowermantle.rho=[7.9565 -6.4761 5.5283 -3.0807 0];
    model.lowermantle.vpv=[24.9520 -40.4673 51.4832 -26.6419 0];
    model.lowermantle.vsv=[11.1671 -13.7818 17.4575 -9.2777 0];
    model.lowermantle.qk=[57823 0 0 0 0];
    model.lowermantle.qu=[312 0 0 0 0];
    model.lowermantle.vph=[24.9520 -40.4673 51.4832 -26.6419 0];
    model.lowermantle.vsh=[11.1671 -13.7818 17.4575 -9.2777 0];
    model.lowermantle.eta=[1 0 0 0 0];
    
    % normal or custom d"?
    switch option.DDPTYPE
        case 'default'
            % d"
            model.order=[model.order {'ddp'}];
            %model.ddp.radius=[3480 3630];
            %model.ddp.depth=[2741 2891];
            model.ddp.radius=[3480 3480+option.DDPTHICKNESS];
            model.ddp.depth=[2891-option.DDPTHICKNESS 2891];
            model.ddp.rho=[7.9565 -6.4761 5.5283 -3.0807 0];
            model.ddp.vpv=[15.3891 -5.3181 5.5242 -2.5514 0];
            model.ddp.vsv=[6.9254 1.4672 -2.0834 0.9783 0];
            model.ddp.qk=[57823 0 0 0 0];
            model.ddp.qu=[312 0 0 0 0];
            model.ddp.vph=[15.3891 -5.3181 5.5242 -2.5514 0];
            model.ddp.vsh=[6.9254 1.4672 -2.0834 0.9783 0];
            model.ddp.eta=[1 0 0 0 0];

            % adjust velocity gradient
            model.ddp.vpv(2:5)=model.ddp.vpv(2:5)*(100+option.DDPVP)/100;
            model.ddp.vsv(2:5)=model.ddp.vsv(2:5)*(100+option.DDPVS)/100;
            model.ddp.vph(2:5)=model.ddp.vph(2:5)*(100+option.DDPVP)/100;
            model.ddp.vsh(2:5)=model.ddp.vsh(2:5)*(100+option.DDPVS)/100;

            % correct for changes
            r=model.ddp.radius(2)/6371;
            model.ddp.vpv(1)=model.ddp.vpv(1) ...
                +poly2val(r,model.lowermantle.vpv) ...
                -poly2val(r,model.ddp.vpv);
            model.ddp.vsv(1)=model.ddp.vsv(1) ...
                +poly2val(r,model.lowermantle.vsv) ...
                -poly2val(r,model.ddp.vsv);
            model.ddp.vph(1)=model.ddp.vph(1) ...
                +poly2val(r,model.lowermantle.vph) ...
                -poly2val(r,model.ddp.vph);
            model.ddp.vsh(1)=model.ddp.vsh(1) ...
                +poly2val(r,model.lowermantle.vsh) ...
                -poly2val(r,model.ddp.vsh);
        case 'custom'
            % how do we approach this?
            % thickness vp_top vp_bottom vp_function
            %           vs_top vs_bottom vs_function
            %           qu_top qu_bottom qu_function
            % function: 1 - linear gradient
            %           2 - parabola
            %           3 - error function
            % - how do supply parameters affecting the shape of erf?
            % - what about a complimentary error function?
            %   - decay to a certain value like wong? yes
            %   - looks like the decay is like erf(0) to erf(2)
            %   - need some way of doing erf2val rather than poly2val
            %       - base on number of values (5=poly)
            % - how about something like what i did for adv geodynamics?
    end
    
    % outer core
    model.order=[model.order {'outercore'}];
    model.outercore.radius=[1221.5 3480];
    model.outercore.depth=[2891 5149.5];
    model.outercore.rho=[12.5815 -1.2638 -3.6526 -5.5281 0];
    model.outercore.vpv=[11.0487 -4.0362 4.8023 -13.5732 0];
    model.outercore.vsv=[0 0 0 0 0];
    model.outercore.qk=[57823 0 0 0 0];
    model.outercore.qu=[99999.9 0 0 0 0];
    model.outercore.vph=[11.0487 -4.0362 4.8023 -13.5732 0];
    model.outercore.vsh=[0 0 0 0 0];
    model.outercore.eta=[1 0 0 0 0];
    
    % inner core
    model.order=[model.order {'innercore'}];
    model.innercore.radius=[0 1221.5];
    model.innercore.depth=[5149.5 6371];
    model.innercore.rho=[13.0885 0 -8.8381 0 0];
    model.innercore.vpv=[11.2622 0 -6.3640 0 0];
    model.innercore.vsv=[3.6678 0 -4.4475 0 0];
    model.innercore.qk=[1327.7 0 0 0 0];
    model.innercore.qu=[84.6 0 0 0 0];
    model.innercore.vph=[11.2622 0 -6.3640 0 0];
    model.innercore.vsh=[3.6678 0 -4.4475 0 0];
    model.innercore.eta=[1 0 0 0 0];
    
    % output structure
    fields={'depth' 'rho' 'vpv' 'vsv' 'qk' 'qu' 'vph' 'vsh' 'eta'};
    vpv=strcmp(fields,'vpv'); vsv=strcmp(fields,'vsv');
    vph=strcmp(fields,'vph'); vsh=strcmp(fields,'vsh');
    qk=strcmp(fields,'qk');   qu=strcmp(fields,'qu');
    
    % specific depths
    if(~isempty(option.DEPTHS))
        % preallocate output
        mout=nan(numel(option.DEPTHS),9);
        
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
                    vsvoigt=sqrt((2*mout(i,vsv)^2+mout(i,vsh)^2)/3);
                    vpvoigt=sqrt((mout(i,vpv)^2+4*mout(i,vph)^2)/5);
                    L=4/3*(vsvoigt/vpvoigt)^2;
                    Ap=1-log(option.REFPERIOD)/pi*...
                        (L/mout(i,qu)+(1-L)/mout(i,qk));
                    As=1-log(option.REFPERIOD)/pi/mout(i,qu);
                    mout(i,vpv)=mout(i,vpv)*Ap;
                    mout(i,vph)=mout(i,vph)*Ap;
                    mout(i,vsv)=mout(i,vsv)*As;
                    mout(i,vsh)=mout(i,vsh)*As;
                    
                    % exit loop (avoids extra output at discontinuities)
                    break;
                end
            end
        end
    else % SPVW at REFPERIOD
        % setup / preallocate output
        mout=nan(1100,9); % 1100 nodes is a guess
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
                r=(6371-depth)/6371;
                mout(n,1)=depth;
                for k=2:numel(fields)
                    mout(n,k)=...
                        poly2val(r,model.(model.order{i}).(fields{k}));
                end

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

function [v]=poly2val(r,p)
v=p(1)+r*(p(2)+r*(p(3)+r*(p(4)+r*p(5))));
end
