function [mout]=prem(varargin)
%PREM    Returns the PREM model
%
%    Usage:    model=prem()
%              model=prem('refperiod',period)
%              model=prem('anisotropic',true|false)
%              model=prem('crust',true|false)
%              model=prem('ocean',true|false)
%              model=prem('depths',depths)
%              model=prem('spvw',nsamples)
%              model=prem('maxdepth',maxdepth)
%
%             Experimental Options:
%              model=prem('ddpthickness',thickness)
%              model=prem('ddpvp',poly)
%              model=prem('ddpvs',poly)
%
%             Not Implemented Yet:
%              model=prem('ddptype',type)
%              model=prem('ulvzthickness',thickness)
%              model=prem('ulvzvelocity',velocity)
%
%    Description:
%
%    Notes:
%     - MODEL: depth rho vpv vsv qk qu vph vsh eta m g
%              depth rho vp vs qk qu m g
%
%       to go: p k u y l dvp_dr dvs_dr d2vp_dr2 d2vs_dr2
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
% - pressure polynomials
% - bulk sound speed derivatives?
% - cmb customization
% - ulvz options

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
option.CRUST=true; % crust
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
    % do we have a crust?
    if(option.CRUST)
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
            model.ocean.drho_dr=0;
            model.ocean.dvpv_dr=0;
            model.ocean.dvsv_dr=0;
            model.ocean.dvph_dr=0;
            model.ocean.dvsh_dr=0;
            model.ocean.deta_dr=0;
            model.ocean.d2rho_dr2=0;
            model.ocean.d2vpv_dr2=0;
            model.ocean.d2vsv_dr2=0;
            model.ocean.d2vph_dr2=0;
            model.ocean.d2vsh_dr2=0;
            model.ocean.d2eta_dr2=0;

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
            model.uppercrust.drho_dr=0;
            model.uppercrust.dvpv_dr=0;
            model.uppercrust.dvsv_dr=0;
            model.uppercrust.dvph_dr=0;
            model.uppercrust.dvsh_dr=0;
            model.uppercrust.deta_dr=0;
            model.uppercrust.d2rho_dr2=0;
            model.uppercrust.d2vpv_dr2=0;
            model.uppercrust.d2vsv_dr2=0;
            model.uppercrust.d2vph_dr2=0;
            model.uppercrust.d2vsh_dr2=0;
            model.uppercrust.d2eta_dr2=0;
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
            model.uppercrust.drho_dr=0;
            model.uppercrust.dvpv_dr=0;
            model.uppercrust.dvsv_dr=0;
            model.uppercrust.dvph_dr=0;
            model.uppercrust.dvsh_dr=0;
            model.uppercrust.deta_dr=0;
            model.uppercrust.d2rho_dr2=0;
            model.uppercrust.d2vpv_dr2=0;
            model.uppercrust.d2vsv_dr2=0;
            model.uppercrust.d2vph_dr2=0;
            model.uppercrust.d2vsh_dr2=0;
            model.uppercrust.d2eta_dr2=0;
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
        model.lowercrust.drho_dr=0;
        model.lowercrust.dvpv_dr=0;
        model.lowercrust.dvsv_dr=0;
        model.lowercrust.dvph_dr=0;
        model.lowercrust.dvsh_dr=0;
        model.lowercrust.deta_dr=0;
        model.lowercrust.d2rho_dr2=0;
        model.lowercrust.d2vpv_dr2=0;
        model.lowercrust.d2vsv_dr2=0;
        model.lowercrust.d2vph_dr2=0;
        model.lowercrust.d2vsh_dr2=0;
        model.lowercrust.d2eta_dr2=0;
        
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
        model.lithosphere.drho_dr=0.6924;
        model.lithosphere.dvpv_dr=7.2180;
        model.lithosphere.dvsv_dr=-1.4678;
        model.lithosphere.dvph_dr=4.6172;
        model.lithosphere.dvsh_dr=5.7176;
        model.lithosphere.deta_dr=-2.4778;
        model.lithosphere.d2rho_dr2=0;
        model.lithosphere.d2vpv_dr2=0;
        model.lithosphere.d2vsv_dr2=0;
        model.lithosphere.d2vph_dr2=0;
        model.lithosphere.d2vsh_dr2=0;
        model.lithosphere.d2eta_dr2=0;
    else % NO CRUST
        % the mantle lithosphere (lid)
        model.order={'lithosphere'};
        model.lithosphere.radius=[6291 6371];
        model.lithosphere.depth=[0 80];
        model.lithosphere.rho=[2.6910 0.6924];
        model.lithosphere.vpv=[0.8317 7.2180];
        model.lithosphere.vsv=[5.8582 -1.4678];
        model.lithosphere.qk=57823;
        model.lithosphere.qu=600;
        model.lithosphere.vph=[3.5908 4.6172];
        model.lithosphere.vsh=[-1.0839 5.7176];
        model.lithosphere.eta=[3.3687 -2.4778];
        model.lithosphere.drho_dr=0.6924;
        model.lithosphere.dvpv_dr=7.2180;
        model.lithosphere.dvsv_dr=-1.4678;
        model.lithosphere.dvph_dr=4.6172;
        model.lithosphere.dvsh_dr=5.7176;
        model.lithosphere.deta_dr=-2.4778;
        model.lithosphere.d2rho_dr2=0;
        model.lithosphere.d2vpv_dr2=0;
        model.lithosphere.d2vsv_dr2=0;
        model.lithosphere.d2vph_dr2=0;
        model.lithosphere.d2vsh_dr2=0;
        model.lithosphere.d2eta_dr2=0;
    end
    
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
    model.asthenosphere.drho_dr=0.6924;
    model.asthenosphere.dvpv_dr=7.2180;
    model.asthenosphere.dvsv_dr=-1.4678;
    model.asthenosphere.dvph_dr=4.6172;
    model.asthenosphere.dvsh_dr=5.7176;
    model.asthenosphere.deta_dr=-2.4778;
    model.asthenosphere.d2rho_dr2=0;
    model.asthenosphere.d2vpv_dr2=0;
    model.asthenosphere.d2vsv_dr2=0;
    model.asthenosphere.d2vph_dr2=0;
    model.asthenosphere.d2vsh_dr2=0;
    model.asthenosphere.d2eta_dr2=0;
    
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
    model.uppermantle.drho_dr=-3.8045;
    model.uppermantle.dvpv_dr=-12.2569;
    model.uppermantle.dvsv_dr=-4.4597;
    model.uppermantle.dvph_dr=-12.2569;
    model.uppermantle.dvsh_dr=-4.4597;
    model.uppermantle.deta_dr=0;
    model.uppermantle.d2rho_dr2=0;
    model.uppermantle.d2vpv_dr2=0;
    model.uppermantle.d2vsv_dr2=0;
    model.uppermantle.d2vph_dr2=0;
    model.uppermantle.d2vsh_dr2=0;
    model.uppermantle.d2eta_dr2=0;
    
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
    model.uppertz.drho_dr=-8.0298;
    model.uppertz.dvpv_dr=-32.6166;
    model.uppertz.dvsv_dr=-18.5856;
    model.uppertz.dvph_dr=-32.6166;
    model.uppertz.dvsh_dr=-18.5856;
    model.uppertz.deta_dr=0;
    model.uppertz.d2rho_dr2=0;
    model.uppertz.d2vpv_dr2=0;
    model.uppertz.d2vsv_dr2=0;
    model.uppertz.d2vph_dr2=0;
    model.uppertz.d2vsh_dr2=0;
    model.uppertz.d2eta_dr2=0;
    
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
    model.middletz.drho_dr=-1.4836;
    model.middletz.dvpv_dr=-9.8672;
    model.middletz.dvsv_dr=-4.9324;
    model.middletz.dvph_dr=-9.8672;
    model.middletz.dvsh_dr=-4.9324;
    model.middletz.deta_dr=0;
    model.middletz.d2rho_dr2=0;
    model.middletz.d2vpv_dr2=0;
    model.middletz.d2vsv_dr2=0;
    model.middletz.d2vph_dr2=0;
    model.middletz.d2vsh_dr2=0;
    model.middletz.d2eta_dr2=0;
    
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
    model.lowertz.drho_dr=[-6.4761 2*5.5283 3*-3.0807];
    model.lowertz.dvpv_dr=[-23.6027 2*5.5242 3*-2.5514];
    model.lowertz.dvsv_dr=[-17.2473 2*-2.0834 3*0.9783];
    model.lowertz.dvph_dr=[-23.6027 2*5.5242 3*-2.5514];
    model.lowertz.dvsh_dr=[-17.2473 2*-2.0834 3*0.9783];
    model.lowertz.deta_dr=0;
    model.lowertz.d2rho_dr2=[2*5.5283 6*-3.0807];
    model.lowertz.d2vpv_dr2=[2*5.5242 6*-2.5514];
    model.lowertz.d2vsv_dr2=[2*-2.0834 6*0.9783];
    model.lowertz.d2vph_dr2=[2*5.5242 6*-2.5514];
    model.lowertz.d2vsh_dr2=[2*-2.0834 6*0.9783];
    model.lowertz.d2eta_dr2=0;
    
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
    model.lowermantle.drho_dr=[-6.4761 2*5.5283 3*-3.0807];
    model.lowermantle.dvpv_dr=[-40.4673 2*51.4832 3*-26.6419];
    model.lowermantle.dvsv_dr=[-13.7818 2*17.4575 3*-9.2777];
    model.lowermantle.dvph_dr=[-40.4673 2*51.4832 3*-26.6419];
    model.lowermantle.dvsh_dr=[-13.7818 2*17.4575 3*-9.2777];
    model.lowermantle.deta_dr=0;
    model.lowermantle.d2rho_dr2=[2*5.5283 6*-3.0807];
    model.lowermantle.d2vpv_dr2=[2*51.4832 6*-26.6419];
    model.lowermantle.d2vsv_dr2=[2*17.4575 6*-9.2777];
    model.lowermantle.d2vph_dr2=[2*51.4832 6*-26.6419];
    model.lowermantle.d2vsh_dr2=[2*17.4575 6*-9.2777];
    model.lowermantle.d2eta_dr2=0;
    
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
            model.ddp.drho_dr=[-6.4761 2*5.5283 3*-3.0807];
            model.ddp.dvpv_dr=[-5.3181 2*5.5242 3*-2.5514];
            model.ddp.dvsv_dr=[1.4672 2*-2.0834 3*0.9783];
            model.ddp.dvph_dr=[-5.3181 2*5.5242 3*-2.5514];
            model.ddp.dvsh_dr=[1.4672 2*-2.0834 3*0.9783];
            model.ddp.deta_dr=0;
            model.ddp.d2rho_dr2=[2*5.5283 6*-3.0807];
            model.ddp.d2vpv_dr2=[2*5.5242 6*-2.5514];
            model.ddp.d2vsv_dr2=[2*-2.0834 6*0.9783];
            model.ddp.d2vph_dr2=[2*5.5242 6*-2.5514];
            model.ddp.d2vsh_dr2=[2*-2.0834 6*0.9783];
            model.ddp.d2eta_dr2=0;

            % adjust velocity gradient
            model.ddp.vpv(2:end)=model.ddp.vpv(2:end)*(100+option.DDPVP)/100;
            model.ddp.vsv(2:end)=model.ddp.vsv(2:end)*(100+option.DDPVS)/100;
            model.ddp.vph(2:end)=model.ddp.vph(2:end)*(100+option.DDPVP)/100;
            model.ddp.vsh(2:end)=model.ddp.vsh(2:end)*(100+option.DDPVS)/100;
            model.ddp.dvpv_dr=model.ddp.dvpv_dr*(100+option.DDPVP)/100;
            model.ddp.dvsv_dr=model.ddp.dvsv_dr*(100+option.DDPVS)/100;
            model.ddp.dvph_dr=model.ddp.dvph_dr*(100+option.DDPVP)/100;
            model.ddp.dvsh_dr=model.ddp.dvsh_dr*(100+option.DDPVS)/100;
            model.ddp.d2vpv_dr2=model.ddp.d2vpv_dr2*(100+option.DDPVP)/100;
            model.ddp.d2vsv_dr2=model.ddp.d2vsv_dr2*(100+option.DDPVS)/100;
            model.ddp.d2vph_dr2=model.ddp.d2vph_dr2*(100+option.DDPVP)/100;
            model.ddp.d2vsh_dr2=model.ddp.d2vsh_dr2*(100+option.DDPVS)/100;

            % correct for changes (remove the introduced discontinuities)
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
    model.outercore.drho_dr=[-1.2638 2*-3.6526 3*-5.5281];
    model.outercore.dvpv_dr=[-4.0362 2*4.8023 3*-13.5732];
    model.outercore.dvsv_dr=0;
    model.outercore.dvph_dr=[-4.0362 2*4.8023 3*-13.5732];
    model.outercore.dvsh_dr=0;
    model.outercore.deta_dr=0;
    model.outercore.d2rho_dr2=[2*-3.6526 6*-5.5281];
    model.outercore.d2vpv_dr2=[2*4.8023 6*-13.5732];
    model.outercore.d2vsv_dr2=0;
    model.outercore.d2vph_dr2=[2*4.8023 6*-13.5732];
    model.outercore.d2vsh_dr2=0;
    model.outercore.d2eta_dr2=0;
    
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
    model.innercore.drho_dr=[0 2*-8.8381];
    model.innercore.dvpv_dr=[0 2*-6.3640];
    model.innercore.dvsv_dr=[0 2*-4.4475];
    model.innercore.dvph_dr=[0 2*-6.3640];
    model.innercore.dvsh_dr=[0 2*-4.4475];
    model.innercore.deta_dr=0;
    model.innercore.d2rho_dr2=2*-8.8381;
    model.innercore.d2vpv_dr2=2*-6.3640;
    model.innercore.d2vsv_dr2=2*-4.4475;
    model.innercore.d2vph_dr2=2*-6.3640;
    model.innercore.d2vsh_dr2=2*-4.4475;
    model.innercore.d2eta_dr2=0;
    
    % get mass polynomials
    %
    %               re
    %              f
    %     m(r)  =  | rho(r)*4*pi*r^2*dr
    %              j
    %             0
    %
    no=numel(model.order);
    for i=1:no
        nm=numel(model.(model.order{i}).rho);
        model.(model.order{i}).m=...
            A*Re3*[0 0 0 model.(model.order{i}).rho./(3:2+nm)];
    end
    
    % correct mass polynomials (remove introduced discontinuities)
    for i=no-1:-1:1
        model.(model.order{i}).m(1)=...
            polyval(fliplr(model.(model.order{i+1}).m),...
            model.(model.order{i+1}).radius(2)/Re)...
            -polyval(fliplr(model.(model.order{i}).m),...
            model.(model.order{i}).radius(1)/Re);
    end
    
    % gravity polynomials start at r^-2 (ie must be divided by r^2)
    %
    %                                     re
    %                                    f
    %     g(r)  =  Gm(r)/r^2  =  G/r^2 * | rho(r)*4*pi*r^2*dr
    %                                    j
    %                                   0
    %
    for i=1:no
        model.(model.order{i}).g=G/Re^2*model.(model.order{i}).m;
    end
    
    % pressure polynomials start at r^-2 (ie must be divided by r^2)
    %
    %                re
    %               f
    %     P(r)  = - | rho(r)*g(r)*dr
    %               j
    %              0
    %
    %tmp=[-2 -1 1 1:10];
    %for i=1:no
        %model.(model.order{i}).rho'
        %model.(model.order{i}).m'
        %model.(model.order{i}).g'
        %model.(model.order{i}).p=...
        %    conv(model.(model.order{i}).g,model.(model.order{i}).rho);
        %nm=numel(model.(model.order{i}).p);
        %model.(model.order{i}).p=1e3*model.(model.order{i}).p./tmp(1:nm);
        %if(nm>2)
        %    model.(model.order{i}).p(3)=log(model.(model.order{i}).p(3));
        %end
    %end
    
    % output structure
    fields=fieldnames(model.(model.order{1}));
    nf=numel(fields);
    
    % specific depths
    if(~isempty(option.DEPTHS))
        % preallocate output (assumes no discon & < MaxDepth)
        for i=1:nf
            mout.(fields{i})=nan(numel(option.DEPTHS),1);
        end
        
        % start counter
        cnt=0;
        
        % loop over depths
        for i=1:numel(option.DEPTHS)
            % skip if below max depth
            if(option.DEPTHS(i)>option.MAXDEPTH); continue; end
            
            % loop over regions
            for j=1:no
                % is depth within this region? (or at the boundary)
                if(option.DEPTHS(i)>=model.(model.order{j}).depth(1) ...
                && option.DEPTHS(i)<=model.(model.order{j}).depth(2))
                    % yes, so increment
                    cnt=cnt+1;
                    
                    % allow for custom d"
                    if(strncmpi(model.order{j},'ddp',3) ...
                            && strcmpi(option.DDPTYPE,'custom'))
                        % custom method
                        
                    else
                        % default polynomial method
                        mout.depth(cnt)=option.DEPTHS(i);
                        mout.radius(cnt)=Re-mout.depth(cnt);
                        r=mout.radius(cnt)/Re;
                        for k=3:nf
                            mout.(fields{k})(cnt)=polyval(fliplr(...
                                model.(model.order{j}).(fields{k})),r);
                        end
                        if(r>0)
                            mout.g(cnt)=mout.g(cnt)/r^2;
                            mout.p(cnt)=mout.p(cnt)/r^2;
                        end
                    end
                    
                    % correct velocities to be at REFPERIOD
                    Lv=4/3*(mout.vsv(cnt)/mout.vpv(cnt))^2;
                    Lh=4/3*(mout.vsh(cnt)/mout.vph(cnt))^2;
                    Apv=1-log(option.REFPERIOD)/pi*...
                        (Lv/mout.qu(cnt)+(1-Lv)/mout.qk(cnt));
                    Aph=1-log(option.REFPERIOD)/pi*...
                        (Lh/mout.qu(cnt)+(1-Lh)/mout.qk(cnt));
                    As=1-log(option.REFPERIOD)/pi/mout.qu(cnt);
                    mout.vpv(cnt)=mout.vpv(cnt)*Apv;
                    mout.vph(cnt)=mout.vph(cnt)*Aph;
                    mout.vsv(cnt)=mout.vsv(cnt)*As;
                    mout.vsh(cnt)=mout.vsh(cnt)*As;
                    mout.dvpv_dr(cnt)=mout.dvpv_dr(cnt)*Apv;
                    mout.dvph_dr(cnt)=mout.dvph_dr(cnt)*Aph;
                    mout.dvsv_dr(cnt)=mout.dvsv_dr(cnt)*As;
                    mout.dvsh_dr(cnt)=mout.dvsh_dr(cnt)*As;
                    mout.d2vpv_dr2(cnt)=mout.d2vpv_dr2(cnt)*Apv;
                    mout.d2vph_dr2(cnt)=mout.d2vph_dr2(cnt)*Aph;
                    mout.d2vsv_dr2(cnt)=mout.d2vsv_dr2(cnt)*As;
                    mout.d2vsh_dr2(cnt)=mout.d2vsh_dr2(cnt)*As;
                    
                    % L is based on isotropic velocities for Voigt solid
                    %vsvoigt=sqrt((2*mout.vsv(cnt)^2+mout.vsh(cnt)^2)/3);
                    %vpvoigt=sqrt((mout.vpv(cnt)^2+4*mout.vph(cnt)^2)/5);
                    %L=4/3*(vsvoigt/vpvoigt)^2;
                    %Ap=1-log(option.REFPERIOD)/pi*...
                    %    (L/mout.qu(cnt)+(1-L)/mout.qk(cnt));
                    %As=1-log(option.REFPERIOD)/pi/mout.qu(cnt);
                    %mout.vpv(cnt)=mout.vpv(cnt)*Ap;
                    %mout.vph(cnt)=mout.vph(cnt)*Ap;
                    %mout.vsv(cnt)=mout.vsv(cnt)*As;
                    %mout.vsh(cnt)=mout.vsh(cnt)*As;
                    %mout.dvpv_dr(cnt)=mout.dvpv_dr(cnt)*Ap;
                    %mout.dvph_dr(cnt)=mout.dvph_dr(cnt)*Ap;
                    %mout.dvsv_dr(cnt)=mout.dvsv_dr(cnt)*As;
                    %mout.dvsh_dr(cnt)=mout.dvsh_dr(cnt)*As;
                    %mout.d2vpv_dr2(cnt)=mout.d2vpv_dr2(cnt)*Ap;
                    %mout.d2vph_dr2(cnt)=mout.d2vph_dr2(cnt)*Ap;
                    %mout.d2vsv_dr2(cnt)=mout.d2vsv_dr2(cnt)*As;
                    %mout.d2vsh_dr2(cnt)=mout.d2vsh_dr2(cnt)*As;
                end
            end
        end
        
        % truncate unused mout
        for i=1:nf
            mout.(fields{i})=mout.(fields{i})(1:cnt);
        end
    else % SPVW at REFPERIOD
        % setup / preallocate output
        for i=1:nf
            mout.(fields{i})=nan(1100,1);
        end
        
        % start counter
        n=1;
        
        % loop over regions
        for i=1:no
            % start at top of region
            depth=model.(model.order{i}).depth(1);
            if(strncmpi(model.order{i},'ddp',3) ...
                    && strcmpi(option.DDPTYPE,'custom'))
                % custom method
                
            else
                % default polynomial method
                mout.depth(n)=depth;
                mout.radius(n)=Re-mout.depth(n);
                r=mout.radius(n)/Re;
                for k=3:nf
                    mout.(fields{k})(n)=polyval(fliplr(...
                        model.(model.order{i}).(fields{k})),r);
                end
                if(r>0); mout.g(n)=mout.g(n)/r^2; end
            end
            
            % correct velocities to be at REFPERIOD
            Lv=4/3*(mout.vsv(n)/mout.vpv(n))^2;
            Lh=4/3*(mout.vsh(n)/mout.vph(n))^2;
            Apv=1-log(option.REFPERIOD)/pi*...
                (Lv/mout.qu(n)+(1-Lv)/mout.qk(n));
            Aph=1-log(option.REFPERIOD)/pi*...
                (Lh/mout.qu(n)+(1-Lh)/mout.qk(n));
            As=1-log(option.REFPERIOD)/pi/mout.qu(n);
            mout.vpv(n)=mout.vpv(n)*Apv;
            mout.vph(n)=mout.vph(n)*Aph;
            mout.vsv(n)=mout.vsv(n)*As;
            mout.vsh(n)=mout.vsh(n)*As;
            mout.dvpv_dr(n)=mout.dvpv_dr(n)*Apv;
            mout.dvph_dr(n)=mout.dvph_dr(n)*Aph;
            mout.dvsv_dr(n)=mout.dvsv_dr(n)*As;
            mout.dvsh_dr(n)=mout.dvsh_dr(n)*As;
            mout.d2vpv_dr2(n)=mout.d2vpv_dr2(n)*Apv;
            mout.d2vph_dr2(n)=mout.d2vph_dr2(n)*Aph;
            mout.d2vsv_dr2(n)=mout.d2vsv_dr2(n)*As;
            mout.d2vsh_dr2(n)=mout.d2vsh_dr2(n)*As;
            
            % L is based on isotropic velocities for Voigt solid
            %vsvoigt=sqrt((2*mout.vsv(n)^2+mout.vsh(n)^2)/3);
            %vpvoigt=sqrt((mout.vpv(n)^2+4*mout.vph(n)^2)/5);
            %L=4/3*(vsvoigt/vpvoigt)^2;
            %Ap=1-log(option.REFPERIOD)/pi*...
            %    (L/mout.qu(n)+(1-L)/mout.qk(n));
            %As=1-log(option.REFPERIOD)/pi/mout.qu(n);
            %mout.vpv(n)=mout.vpv(n)*Ap;
            %mout.vph(n)=mout.vph(n)*Ap;
            %mout.vsv(n)=mout.vsv(n)*As;
            %mout.vsh(n)=mout.vsh(n)*As;
            %mout.dvpv_dr(n)=mout.dvpv_dr(n)*Ap;
            %mout.dvph_dr(n)=mout.dvph_dr(n)*Ap;
            %mout.dvsv_dr(n)=mout.dvsv_dr(n)*As;
            %mout.dvsh_dr(n)=mout.dvsh_dr(n)*As;
            %mout.d2vpv_dr2(n)=mout.d2vpv_dr2(n)*Ap;
            %mout.d2vph_dr2(n)=mout.d2vph_dr2(n)*Ap;
            %mout.d2vsv_dr2(n)=mout.d2vsv_dr2(n)*As;
            %mout.d2vsh_dr2(n)=mout.d2vsh_dr2(n)*As;
            
            % step to next depth using lowest vertical velocity
            if(mout.vsv(n)>0)
                vmin=mout.vsv(n);
            else
                vmin=mout.vpv(n);
            end
            wavelength=vmin*option.REFPERIOD;
            stepsize=wavelength/option.SPVW;
            depth=depth+stepsize;
            n=n+1;
            
            % finished if below max depth
            if(depth>option.MAXDEPTH); break; end
            
            % step until at/below bottom of region (within a meter)
            while(depth<(model.(model.order{i}).depth(2)-0.001))
                if(strncmpi(model.order{i},'ddp',3) ...
                        && strcmpi(option.DDPTYPE,'custom'))
                    % custom method
                    
                else
                    % default polynomial method
                    mout.depth(n)=depth;
                    mout.radius(n)=Re-mout.depth(n);
                    r=mout.radius(n)/Re;
                    for k=3:nf
                        mout.(fields{k})(n)=polyval(fliplr(...
                            model.(model.order{i}).(fields{k})),r);
                    end
                    if(r>0); mout.g(n)=mout.g(n)/r^2; end
                end

                % correct velocities to be at REFPERIOD
                Lv=4/3*(mout.vsv(n)/mout.vpv(n))^2;
                Lh=4/3*(mout.vsh(n)/mout.vph(n))^2;
                Apv=1-log(option.REFPERIOD)/pi*...
                    (Lv/mout.qu(n)+(1-Lv)/mout.qk(n));
                Aph=1-log(option.REFPERIOD)/pi*...
                    (Lh/mout.qu(n)+(1-Lh)/mout.qk(n));
                As=1-log(option.REFPERIOD)/pi/mout.qu(n);
                mout.vpv(n)=mout.vpv(n)*Apv;
                mout.vph(n)=mout.vph(n)*Aph;
                mout.vsv(n)=mout.vsv(n)*As;
                mout.vsh(n)=mout.vsh(n)*As;
                mout.dvpv_dr(n)=mout.dvpv_dr(n)*Apv;
                mout.dvph_dr(n)=mout.dvph_dr(n)*Aph;
                mout.dvsv_dr(n)=mout.dvsv_dr(n)*As;
                mout.dvsh_dr(n)=mout.dvsh_dr(n)*As;
                mout.d2vpv_dr2(n)=mout.d2vpv_dr2(n)*Apv;
                mout.d2vph_dr2(n)=mout.d2vph_dr2(n)*Aph;
                mout.d2vsv_dr2(n)=mout.d2vsv_dr2(n)*As;
                mout.d2vsh_dr2(n)=mout.d2vsh_dr2(n)*As;
                
                % step to next depth using lowest vertical velocity
                if(mout.vsv(n)>0)
                    vmin=mout.vsv(n);
                else
                    vmin=mout.vpv(n);
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
            if(strncmpi(model.order{i},'ddp',3) ...
                    && strcmpi(option.DDPTYPE,'custom'))
                % custom method
                
            else
                % default polynomial method
                mout.depth(n)=depth;
                mout.radius(n)=Re-mout.depth(n);
                r=mout.radius(n)/Re;
                for k=3:nf
                    mout.(fields{k})(n)=polyval(fliplr(...
                        model.(model.order{i}).(fields{k})),r);
                end
                if(r>0); mout.g(n)=mout.g(n)/r^2; end
            end
            
            % correct velocities to be at REFPERIOD
            Lv=4/3*(mout.vsv(n)/mout.vpv(n))^2;
            Lh=4/3*(mout.vsh(n)/mout.vph(n))^2;
            Apv=1-log(option.REFPERIOD)/pi*...
                (Lv/mout.qu(n)+(1-Lv)/mout.qk(n));
            Aph=1-log(option.REFPERIOD)/pi*...
                (Lh/mout.qu(n)+(1-Lh)/mout.qk(n));
            As=1-log(option.REFPERIOD)/pi/mout.qu(n);
            mout.vpv(n)=mout.vpv(n)*Apv;
            mout.vph(n)=mout.vph(n)*Aph;
            mout.vsv(n)=mout.vsv(n)*As;
            mout.vsh(n)=mout.vsh(n)*As;
            mout.dvpv_dr(n)=mout.dvpv_dr(n)*Apv;
            mout.dvph_dr(n)=mout.dvph_dr(n)*Aph;
            mout.dvsv_dr(n)=mout.dvsv_dr(n)*As;
            mout.dvsh_dr(n)=mout.dvsh_dr(n)*As;
            mout.d2vpv_dr2(n)=mout.d2vpv_dr2(n)*Apv;
            mout.d2vph_dr2(n)=mout.d2vph_dr2(n)*Aph;
            mout.d2vsv_dr2(n)=mout.d2vsv_dr2(n)*As;
            mout.d2vsh_dr2(n)=mout.d2vsh_dr2(n)*As;
            
            % increment
            n=n+1;
        end
        
        % truncate unused mout
        n=n-1;
        for i=1:nf
            mout.(fields{i})=mout.(fields{i})(1:n);
        end
    end
    
    % bulk sound speed
    mout.vbv=sqrt(mout.vpv.^2-4/3*mout.vsv.^2);
    mout.vbh=sqrt(mout.vph.^2-4/3*mout.vsh.^2);
    
    % now get moduli
    mout.poissonv=0.5*(mout.vpv.^2-2*mout.vsv.^2)./(mout.vpv.^2-mout.vsv.^2);
    mout.poissonh=0.5*(mout.vph.^2-2*mout.vsh.^2)./(mout.vph.^2-mout.vsh.^2);
    mout.shearv=mout.vsv.^2.*mout.rho*1e9;
    mout.shearh=mout.vsh.^2.*mout.rho*1e9;
    mout.bulkv=mout.vpv.^2.*mout.rho*1e9-4/3*mout.shearv;
    mout.bulkh=mout.vph.^2.*mout.rho*1e9-4/3*mout.shearh;
    mout.youngsv=2*mout.shearv.*(1+mout.poissonv);
    mout.youngsh=2*mout.shearh.*(1+mout.poissonh);
    mout.lambdav=mout.bulkv-2/3*mout.shearv;
    mout.lambdah=mout.bulkh-2/3*mout.shearh;
else % ISOTROPIC
    % do we have a crust?
    if(option.CRUST)
        % do we have an ocean?
        if(option.OCEAN) % YES
            % the ocean (3km)
            model.order={'ocean'};
            model.ocean.radius=[6368 6371];
            model.ocean.depth=[0 3];
            model.ocean.rho=1.02;
            model.ocean.vp=1.45;
            model.ocean.vs=0;
            model.ocean.qk=57823;
            model.ocean.qu=99999.9;
            model.ocean.drho_dr=0;
            model.ocean.dvp_dr=0;
            model.ocean.dvs_dr=0;
            model.ocean.d2rho_dr2=0;
            model.ocean.d2vp_dr2=0;
            model.ocean.d2vs_dr2=0;

            % the upper crust
            model.order=[model.order {'uppercrust'}];
            model.uppercrust.radius=[6356 6368];
            model.uppercrust.depth=[3 15];
            model.uppercrust.rho=2.6;
            model.uppercrust.vp=5.8;
            model.uppercrust.vs=3.2;
            model.uppercrust.qk=57823;
            model.uppercrust.qu=600;
            model.uppercrust.drho_dr=0;
            model.uppercrust.dvp_dr=0;
            model.uppercrust.dvs_dr=0;
            model.uppercrust.d2rho_dr2=0;
            model.uppercrust.d2vp_dr2=0;
            model.uppercrust.d2vs_dr2=0;
        else % NO OCEAN
            % the upper crust (to the surface)
            model.order={'uppercrust'};
            model.uppercrust.radius=[6356 6371];
            model.uppercrust.depth=[0 15];
            model.uppercrust.rho=2.6;
            model.uppercrust.vp=5.8;
            model.uppercrust.vs=3.2;
            model.uppercrust.qk=57823;
            model.uppercrust.qu=600;
            model.uppercrust.drho_dr=0;
            model.uppercrust.dvp_dr=0;
            model.uppercrust.dvs_dr=0;
            model.uppercrust.d2rho_dr2=0;
            model.uppercrust.d2vp_dr2=0;
            model.uppercrust.d2vs_dr2=0;
        end

        % the lower crust
        model.order=[model.order {'lowercrust'}];
        model.lowercrust.radius=[6346.6 6356];
        model.lowercrust.depth=[15 24.4];
        model.lowercrust.rho=2.9;
        model.lowercrust.vp=6.8;
        model.lowercrust.vs=3.9;
        model.lowercrust.qk=57823;
        model.lowercrust.qu=600;
        model.lowercrust.drho_dr=0;
        model.lowercrust.dvp_dr=0;
        model.lowercrust.dvs_dr=0;
        model.lowercrust.d2rho_dr2=0;
        model.lowercrust.d2vp_dr2=0;
        model.lowercrust.d2vs_dr2=0;

        % the mantle lithosphere (lid)
        model.order=[model.order {'lithosphere'}];
        model.lithosphere.radius=[6291 6346.6];
        model.lithosphere.depth=[24.4 80];
        model.lithosphere.rho=[2.6910 0.6924];
        model.lithosphere.vp=[4.1875 3.9382];
        model.lithosphere.vs=[2.1519 2.3481];
        model.lithosphere.qk=57823;
        model.lithosphere.qu=600;
        model.lithosphere.drho_dr=0.6924;
        model.lithosphere.dvp_dr=3.9382;
        model.lithosphere.dvs_dr=2.3481;
        model.lithosphere.d2rho_dr2=0;
        model.lithosphere.d2vp_dr2=0;
        model.lithosphere.d2vs_dr2=0;
    else
        % the mantle lithosphere (lid)
        model.order={'lithosphere'};
        model.lithosphere.radius=[6291 6371];
        model.lithosphere.depth=[0 80];
        model.lithosphere.rho=[2.6910 0.6924];
        model.lithosphere.vp=[4.1875 3.9382];
        model.lithosphere.vs=[2.1519 2.3481];
        model.lithosphere.qk=57823;
        model.lithosphere.qu=600;
        model.lithosphere.drho_dr=0.6924;
        model.lithosphere.dvp_dr=3.9382;
        model.lithosphere.dvs_dr=2.3481;
        model.lithosphere.d2rho_dr2=0;
        model.lithosphere.d2vp_dr2=0;
        model.lithosphere.d2vs_dr2=0;
    end
    
    % the mantle asthenosphere (lvz)
    model.order=[model.order {'asthenosphere'}];
    model.asthenosphere.radius=[6151 6291];
    model.asthenosphere.depth=[80 220];
    model.asthenosphere.rho=[2.6910 0.6924];
    model.asthenosphere.vp=[4.1875 3.9382];
    model.asthenosphere.vs=[2.1519 2.3481];
    model.asthenosphere.qk=57823;
    model.asthenosphere.qu=80;
    model.asthenosphere.drho_dr=0.6924;
    model.asthenosphere.dvp_dr=3.9382;
    model.asthenosphere.dvs_dr=2.3481;
    model.asthenosphere.d2rho_dr2=0;
    model.asthenosphere.d2vp_dr2=0;
    model.asthenosphere.d2vs_dr2=0;
    
    % the upper mantle
    model.order=[model.order {'uppermantle'}];
    model.uppermantle.radius=[5971 6151];
    model.uppermantle.depth=[220 400];
    model.uppermantle.rho=[7.1089 -3.8045];
    model.uppermantle.vp=[20.3926 -12.2569];
    model.uppermantle.vs=[8.9496 -4.4597];
    model.uppermantle.qk=57823;
    model.uppermantle.qu=143;
    model.uppermantle.drho_dr=-3.8045;
    model.uppermantle.dvp_dr=-12.2569;
    model.uppermantle.dvs_dr=-4.4597;
    model.uppermantle.d2rho_dr2=0;
    model.uppermantle.d2vp_dr2=0;
    model.uppermantle.d2vs_dr2=0;
    
    % upper transition zone
    model.order=[model.order {'uppertz'}];
    model.uppertz.radius=[5771 5971];
    model.uppertz.depth=[400 600];
    model.uppertz.rho=[11.2494 -8.0298];
    model.uppertz.vp=[39.7027 -32.6166];
    model.uppertz.vs=[22.3512 -18.5856];
    model.uppertz.qk=57823;
    model.uppertz.qu=143;
    model.uppertz.drho_dr=-8.0298;
    model.uppertz.dvp_dr=-32.6166;
    model.uppertz.dvs_dr=-18.5856;
    model.uppertz.d2rho_dr2=0;
    model.uppertz.d2vp_dr2=0;
    model.uppertz.d2vs_dr2=0;
    
    % middle transition zone
    model.order=[model.order {'middletz'}];
    model.middletz.radius=[5701 5771];
    model.middletz.depth=[600 670];
    model.middletz.rho=[5.3197 -1.4836];
    model.middletz.vp=[19.0957 -9.8672];
    model.middletz.vs=[9.9839 -4.9324];
    model.middletz.qk=57823;
    model.middletz.qu=143;
    model.middletz.drho_dr=-1.4836;
    model.middletz.dvp_dr=-9.8672;
    model.middletz.dvs_dr=-4.9324;
    model.middletz.d2rho_dr2=0;
    model.middletz.d2vp_dr2=0;
    model.middletz.d2vs_dr2=0;
    
    % lower transition zone
    model.order=[model.order {'lowertz'}];
    model.lowertz.radius=[5600 5701];
    model.lowertz.depth=[670 771];
    model.lowertz.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.lowertz.vp=[29.2766 -23.6027 5.5242 -2.5514];
    model.lowertz.vs=[22.3459 -17.2473 -2.0834 0.9783];
    model.lowertz.qk=57823;
    model.lowertz.qu=312;
    model.lowertz.drho_dr=[-6.4761 2*5.5283 3*-3.0807];
    model.lowertz.dvp_dr=[-23.6027 2*5.5242 3*-2.5514];
    model.lowertz.dvs_dr=[-17.2473 2*-2.0834 3*0.9783];
    model.lowertz.d2rho_dr2=[2*5.5283 6*-3.0807];
    model.lowertz.d2vp_dr2=[2*5.5242 6*-2.5514];
    model.lowertz.d2vs_dr2=[2*-2.0834 6*0.9783];
    
    % lower mantle
    model.order=[model.order {'lowermantle'}];
    %model.lowermantle.radius=[3630 5600];
    %model.lowermantle.depth=[771 2741];
    model.lowermantle.radius=[3480+option.DDPTHICKNESS 5600];
    model.lowermantle.depth=[771 2891-option.DDPTHICKNESS];
    model.lowermantle.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.lowermantle.vp=[24.9520 -40.4673 51.4832 -26.6419];
    model.lowermantle.vs=[11.1671 -13.7818 17.4575 -9.2777];
    model.lowermantle.qk=57823;
    model.lowermantle.qu=312;
    model.lowermantle.drho_dr=[-6.4761 2*5.5283 3*-3.0807];
    model.lowermantle.dvp_dr=[-40.4673 2*51.4832 3*-26.6419];
    model.lowermantle.dvs_dr=[-13.7818 2*17.4575 3*-9.2777];
    model.lowermantle.d2rho_dr2=[2*5.5283 6*-3.0807];
    model.lowermantle.d2vp_dr2=[2*51.4832 6*-26.6419];
    model.lowermantle.d2vs_dr2=[2*17.4575 6*-9.2777];
    
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
            model.ddp.vp=[15.3891 -5.3181 5.5242 -2.5514];
            model.ddp.vs=[6.9254 1.4672 -2.0834 0.9783];
            model.ddp.qk=57823;
            model.ddp.qu=312;
            model.ddp.drho_dr=[-6.4761 2*5.5283 3*-3.0807];
            model.ddp.dvp_dr=[-5.3181 2*5.5242 3*-2.5514];
            model.ddp.dvs_dr=[1.4672 2*-2.0834 3*0.9783];
            model.ddp.d2rho_dr2=[2*5.5283 6*-3.0807];
            model.ddp.d2vp_dr2=[2*5.5242 6*-2.5514];
            model.ddp.d2vs_dr2=[2*-2.0834 6*0.9783];

            % adjust velocity gradient
            model.ddp.vp(2:end)=model.ddp.vp(2:end)*(100+option.DDPVP)/100;
            model.ddp.vs(2:end)=model.ddp.vs(2:end)*(100+option.DDPVS)/100;
            model.ddp.dvp_dr=model.ddp.dvp_dr*(100+option.DDPVP)/100;
            model.ddp.dvs_dr=model.ddp.dvs_dr*(100+option.DDPVS)/100;
            model.ddp.d2vp_dr2=model.ddp.d2vp_dr2*(100+option.DDPVP)/100;
            model.ddp.d2vs_dr2=model.ddp.d2vs_dr2*(100+option.DDPVS)/100;

            % correct for changes (remove the introduced discontinuities)
            r=model.ddp.radius(2)/6371;
            model.ddp.vp(1)=model.ddp.vp(1) ...
                +polyval(fliplr(model.lowermantle.vp),r) ...
                -polyval(fliplr(model.ddp.vp),r);
            model.ddp.vs(1)=model.ddp.vs(1) ...
                +polyval(fliplr(model.lowermantle.vs),r) ...
                -polyval(fliplr(model.ddp.vs),r);
        case 'custom'
            % see anisotropy section
    end
    
    % outer core
    model.order=[model.order {'outercore'}];
    model.outercore.radius=[1221.5 3480];
    model.outercore.depth=[2891 5149.5];
    model.outercore.rho=[12.5815 -1.2638 -3.6526 -5.5281];
    model.outercore.vp=[11.0487 -4.0362 4.8023 -13.5732];
    model.outercore.vs=0;
    model.outercore.qk=57823;
    model.outercore.qu=99999.9;
    model.outercore.drho_dr=[-1.2638 2*-3.6526 3*-5.5281];
    model.outercore.dvp_dr=[-4.0362 2*4.8023 3*-13.5732];
    model.outercore.dvs_dr=0;
    model.outercore.d2rho_dr2=[2*-3.6526 6*-5.5281];
    model.outercore.d2vp_dr2=[2*4.8023 6*-13.5732];
    model.outercore.d2vs_dr2=0;
    
    % inner core
    model.order=[model.order {'innercore'}];
    model.innercore.radius=[0 1221.5];
    model.innercore.depth=[5149.5 6371];
    model.innercore.rho=[13.0885 0 -8.8381];
    model.innercore.vp=[11.2622 0 -6.3640];
    model.innercore.vs=[3.6678 0 -4.4475];
    model.innercore.qk=1327.7;
    model.innercore.qu=84.6;
    model.innercore.drho_dr=[0 2*-8.8381];
    model.innercore.dvp_dr=[0 2*-6.3640];
    model.innercore.dvs_dr=[0 2*-4.4475];
    model.innercore.d2rho_dr2=2*-8.8381;
    model.innercore.d2vp_dr2=2*-6.3640;
    model.innercore.d2vs_dr2=2*-4.4475;
    
    % get mass polynomials
    %
    %               re
    %              f
    %     m(r)  =  | p(r)*4*pi*r^2*dr
    %              j
    %             0
    %
    no=numel(model.order);
    for i=1:no
        % we assume rho is a vector of polynomial constants
        % applied to r^0 thru r^(nm-1)
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
    %
    %                           re
    %                          f
    %     g(r)  =  Gm(r)/r  =  | G*p(r)*4*pi*r*dr
    %                          j
    %                         0
    %
    for i=1:no
        model.(model.order{i}).g=G/Re^2*model.(model.order{i}).m;
    end
    
    % output structure
    fields=fieldnames(model.(model.order{1}));
    nf=numel(fields);
    
    % specific depths
    if(~isempty(option.DEPTHS))
        % preallocate output (assumes no discon & < MaxDepth)
        for i=1:nf
            mout.(fields{i})=nan(numel(option.DEPTHS),1);
        end
        
        % start counter
        cnt=0;
        
        % loop over depths
        for i=1:numel(option.DEPTHS)
            % skip if below max depth
            if(option.DEPTHS(i)>option.MAXDEPTH); continue; end
            
            % loop over regions
            for j=1:no
                % is depth within this region? (or at the boundary)
                if(option.DEPTHS(i)>=model.(model.order{j}).depth(1) ...
                && option.DEPTHS(i)<=model.(model.order{j}).depth(2))
                    % yes, so increment
                    cnt=cnt+1;
                    
                    % allow for custom d"
                    if(strncmpi(model.order{j},'ddp',3) ...
                            && strcmpi(option.DDPTYPE,'custom'))
                        % custom method
                        
                    else
                        % default polynomial method
                        mout.depth(cnt)=option.DEPTHS(i);
                        mout.radius(cnt)=Re-mout.depth(cnt);
                        r=mout.radius(cnt)/Re;
                        for k=3:nf
                            mout.(fields{k})(cnt)=polyval(fliplr(...
                                model.(model.order{j}).(fields{k})),r);
                        end
                        if(r>0); mout.g(cnt)=mout.g(cnt)/r^2; end
                    end
                    
                    % correct velocities to be at REFPERIOD
                    L=4/3*(mout.vs(cnt)/mout.vp(cnt))^2;
                    Ap=1-log(option.REFPERIOD)/pi*...
                        (L/mout.qu(cnt)+(1-L)/mout.qk(cnt));
                    As=1-log(option.REFPERIOD)/pi/mout.qu(cnt);
                    mout.vp(cnt)=mout.vp(cnt)*Ap;
                    mout.vs(cnt)=mout.vs(cnt)*As;
                    mout.dvp_dr(cnt)=mout.dvp_dr(cnt)*Ap;
                    mout.dvs_dr(cnt)=mout.dvs_dr(cnt)*As;
                    mout.d2vp_dr2(cnt)=mout.d2vp_dr2(cnt)*Ap;
                    mout.d2vs_dr2(cnt)=mout.d2vs_dr2(cnt)*As;
                end
            end
        end
        
        % truncate unused mout
        for i=1:nf
            mout.(fields{i})=mout.(fields{i})(1:cnt);
        end
    else % SPVW at REFPERIOD
        % setup / preallocate output
        for i=1:nf
            mout.(fields{i})=nan(1100,1);
        end
        
        % start counter
        n=1;
        
        % loop over regions
        for i=1:no
            % start at top of region
            depth=model.(model.order{i}).depth(1);
            if(strncmpi(model.order{i},'ddp',3) ...
                    && strcmpi(option.DDPTYPE,'custom'))
                % custom method
                
            else
                % default polynomial method
                mout.depth(n)=depth;
                mout.radius(n)=Re-mout.depth(n);
                r=mout.radius(n)/Re;
                for k=3:nf
                    mout.(fields{k})(n)=polyval(fliplr(...
                        model.(model.order{i}).(fields{k})),r);
                end
                if(r>0); mout.g(n)=mout.g(n)/r^2; end
            end
            
            % correct velocities to be at REFPERIOD
            L=4/3*(mout.vs(n)/mout.vp(n))^2;
            Ap=1-log(option.REFPERIOD)/pi*...
                (L/mout.qu(n)+(1-L)/mout.qk(n));
            As=1-log(option.REFPERIOD)/pi/mout.qu(n);
            mout.vp(n)=mout.vp(n)*Ap;
            mout.vs(n)=mout.vs(n)*As;
            mout.dvp_dr(n)=mout.dvp_dr(n)*Ap;
            mout.dvs_dr(n)=mout.dvs_dr(n)*As;
            mout.d2vp_dr2(n)=mout.d2vp_dr2(n)*Ap;
            mout.d2vs_dr2(n)=mout.d2vs_dr2(n)*As;
            
            % step to next depth using lowest vertical velocity
            if(mout.vs(n)>0)
                vmin=mout.vs(n);
            else
                vmin=mout.vp(n);
            end
            wavelength=vmin*option.REFPERIOD;
            stepsize=wavelength/option.SPVW;
            depth=depth+stepsize;
            n=n+1;
            
            % finished if below max depth
            if(depth>option.MAXDEPTH); break; end
            
            % step until at/below bottom of region (within a meter)
            while(depth<(model.(model.order{i}).depth(2)-0.001))
                if(strncmpi(model.order{i},'ddp',3) ...
                        && strcmpi(option.DDPTYPE,'custom'))
                    % custom method
                    
                else
                    % default polynomial method
                    mout.depth(n)=depth;
                    mout.radius(n)=Re-mout.depth(n);
                    r=mout.radius(n)/Re;
                    for k=3:nf
                        mout.(fields{k})(n)=polyval(fliplr(...
                            model.(model.order{i}).(fields{k})),r);
                    end
                    if(r>0); mout.g(n)=mout.g(n)/r^2; end
                end

                % correct velocities to be at REFPERIOD
                L=4/3*(mout.vs(n)/mout.vp(n))^2;
                Ap=1-log(option.REFPERIOD)/pi*...
                    (L/mout.qu(n)+(1-L)/mout.qk(n));
                As=1-log(option.REFPERIOD)/pi/mout.qu(n);
                mout.vp(n)=mout.vp(n)*Ap;
                mout.vs(n)=mout.vs(n)*As;
                mout.dvp_dr(n)=mout.dvp_dr(n)*Ap;
                mout.dvs_dr(n)=mout.dvs_dr(n)*As;
                mout.d2vp_dr2(n)=mout.d2vp_dr2(n)*Ap;
                mout.d2vs_dr2(n)=mout.d2vs_dr2(n)*As;
                
                % step to next depth using lowest vertical velocity
                if(mout.vs(n)>0)
                    vmin=mout.vs(n);
                else
                    vmin=mout.vp(n);
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
            if(strncmpi(model.order{i},'ddp',3) ...
                    && strcmpi(option.DDPTYPE,'custom'))
                % custom method
                
            else
                % default polynomial method
                mout.depth(n)=depth;
                mout.radius(n)=Re-mout.depth(n);
                r=mout.radius(n)/Re;
                for k=3:nf
                    mout.(fields{k})(n)=polyval(fliplr(...
                        model.(model.order{i}).(fields{k})),r);
                end
                if(r>0); mout.g(n)=mout.g(n)/r^2; end
            end
            
            % correct velocities to be at REFPERIOD
            L=4/3*(mout.vs(n)/mout.vp(n))^2;
            Ap=1-log(option.REFPERIOD)/pi*...
                (L/mout.qu(n)+(1-L)/mout.qk(n));
            As=1-log(option.REFPERIOD)/pi/mout.qu(n);
            mout.vp(n)=mout.vp(n)*Ap;
            mout.vs(n)=mout.vs(n)*As;
            mout.dvp_dr(n)=mout.dvp_dr(n)*Ap;
            mout.dvs_dr(n)=mout.dvs_dr(n)*As;
            mout.d2vp_dr2(n)=mout.d2vp_dr2(n)*Ap;
            mout.d2vs_dr2(n)=mout.d2vs_dr2(n)*As;
            
            % increment
            n=n+1;
        end
        
        % truncate unused mout
        n=n-1;
        for i=1:nf
            mout.(fields{i})=mout.(fields{i})(1:n);
        end
    end
    
    % bulk sound speed
    mout.vb=sqrt(mout.vp.^2-4/3*mout.vs.^2);
    
    % now get moduli
    mout.poisson=0.5*(mout.vp.^2-2*mout.vs.^2)./(mout.vp.^2-mout.vs.^2);
    mout.shear=mout.vs.^2.*mout.rho*1e9;
    mout.bulk=mout.vp.^2.*mout.rho*1e9-4/3*mout.shear;
    mout.youngs=2*mout.shear.*(1+mout.poisson);
    mout.lambda=mout.bulk-2/3*mout.shear;
end

end
