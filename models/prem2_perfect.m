function [mout]=prem2_perfect(varargin)
%PREM2_PERFECT    Returns the PREM2 model
%
%    Usage:    model=prem2_perfect()
%              model=prem2_perfect('refperiod',period)
%              model=prem2_perfect('anisotropic',true|false)
%              model=prem2_perfect('crust',true|false)
%              model=prem2_perfect('ocean',true|false)
%              model=prem2_perfect('depths',depths)
%              model=prem2_perfect('spvw',nsamples)
%              model=prem2_perfect('maxdepth',maxdepth)
%
%    Description:
%     MODEL=PREM2_PERFECT() returns the PREM2 model of Song & Helmberger,
%     which is a slight modification of the p-wave speeds (vp & vb) near
%     the core-mantle boundary and the inner-outer core boundary.  The
%     model is at a reference period of 1Hz & is sampled once per vertical
%     wavelength for the entire Earth.  The model is anisotropic and
%     includes the oceanic and crustal layers.  The output includes
%     velocities, moduli, quality factors, mass, density and gravity.  Note
%     that all discontinuities are included and are sampled on both sides.
%
%     MODEL=PREM2_PERFECT('REFPERIOD',PERIOD) changes the reference period
%     of the model.  This will alter the velocites returned as they are
%     compensated for physical dispersion using Azimi's law.  You may use
%     the 'REFFREQUENCY' option if you prefer.  The default is 1s.
%
%     MODEL=PREM2_PERFECT('ANISOTROPIC',TRUE|FALSE) sets if the anisotropic
%     or isotropic version is returned.  You may use the 'ISOTROPIC' option
%     if you prefer.  The default is TRUE.
%
%     MODEL=PREM2_PERFECT('CRUST',TRUE|FALSE) sets if PREM2 includes
%     the crust layer.  If FALSE, then the layer below (the mantle) is
%     extended to the surface (no ocean too).  The default is TRUE.
%
%     MODEL=PREM2_PERFECT('OCEAN',TRUE|FALSE) sets if PREM2 includes
%     the ocean layer.  If FALSE, then the layer below (the crust) is
%     extended to the surface.  Please note that to include the ocean both
%     the OCEAN and CRUST options must be TRUE.  The default is TRUE.
%
%     MODEL=PREM2_PERFECT('DEPTHS',DEPTHS) indicates specific depths (in
%     km) to return the PREM model parameters at.  The default is [] (empty
%     matrix) and will instead use the SPVW option to sample the PREM2
%     model.  Note that depths specified at discontinuities will return
%     values for both sides.
%
%     MODEL=PREM2_PERFECT('SPVW',NSAMPLES) number of samples (or knots if
%     you prefer) per vertical wavelength to sample the PREM2 model at.
%     The default is 1.  This option is ignored if the DEPTHS option is
%     non-empty.  The vertical wavelength is determined by the reference
%     period (set with the REFPERIOD option) and the vertical s-wave
%     velocity (except when that quantity is 0km/s in the liquid layers, in
%     which case the vertical p-wave velocity is used).  This is useful for
%     sampling the model appropriately for numerical calculations.
%
%     MODEL=PREM2_PERFECT('MAXDEPTH',MAXDEPTH) is the maximum depth to
%     sample the model to.  The default is [] (empty matrix) which will
%     sample the entire model.
%
%    Notes:
%     - PREM reference:
%        Dziewonski & Anderson 1981, Preliminary reference earth model,
%        Phys. Earth planet. Inter. 25, pp. 297-356
%     - PREM2 reference:
%        Song & Helmberger 1995, A P-wave model of the Earth's core, J.
%        geophys. Res. 100, pp. 1917-1930
%     - Missing pressure values.
%     - It is slow, I know.  The implementation isn't optimal for speed.
%
%    Examples:
%     % compare PREM & PREM2
%     model(1)=prem_perfect;
%     model(2)=prem2_perfect;
%     plot1dmodel(model);
%
%    See also: PREM_PERFECT, PREM, IASP91, AK135

%     Version History:
%        Aug. 17, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2010 at 07:40 GMT

% todo:

% id funtion
me=mfilename;
pkgme=['seizmo:' me ':'];

% all inputs must be 'option',value pairs
if(mod(nargin,2))
    error([pkgme 'unpairedOption'],...
        'Option missing a value!');
end

% define some constants
A=pi*4e12;
G=6.67428e-17;
Re=6371;
Re3=Re^3;

% OPTION DEFAULTS:
% reference period of 1s
% anisotropy
% ocean
% crust
% no predefined depths
% minimum of 1 sample(s) per vertical wavelength
% no max depth
varargin=[{'rp' 1 'ani' true 'oc' true 'cru' true 'dep' [] 'spvw' 1 ...
    'maxdep' []} varargin];

% all options must be specified with a string
if(~iscellstr(varargin(1:2:end)))
    error([pkgme 'invalidOption'],...
        'Options must be specified with a string!');
end

% check options
for i=1:2:numel(varargin)
    val=varargin{i+1};
    switch lower(varargin{i})
        case {'p' 'rp' 'refper' 'period' 'per' 'refperiod'}
            if(isempty(val)); continue; end
            if(~isreal(val) || ~isscalar(val) || val<=0)
                error([pkgme 'badInput'],...
                    'REFPERIOD must be a positive real scalar!');
            end
            option.REFPERIOD=val;
        case {'f' 'rf' 'reffreq' 'frequency' 'freq' 'reffrequency'}
            if(isempty(val)); continue; end
            if(~isreal(val) || ~isscalar(val) || val<=0)
                error([pkgme 'badInput'],...
                    'REFFREQUENCY must be a positive real scalar!');
            end
            option.REFPERIOD=1/val;
        case {'a' 'ani' 'anisotropy' 'anisotropic'}
            if(isempty(val)); continue; end
            if(~isscalar(val) || (~isreal(val) && ~islogical(val)))
                error([pkgme 'badInput'],...
                    'ANISOTROPY must be TRUE or FALSE!');
            end
            option.ANISOTROPY=val;
        case {'i' 'iso' 'isotropy' 'isotropic'}
            if(isempty(val)); continue; end
            if(~isscalar(val) || (~isreal(val) && ~islogical(val)))
                error([pkgme 'badInput'],...
                    'ISOTROPY must be TRUE or FALSE!');
            end
            option.ANISOTROPY=~val;
        case {'o' 'oc' 'ocean'}
            if(isempty(val)); continue; end
            if(~isscalar(val) || (~isreal(val) && ~islogical(val)))
                error([pkgme 'badInput'],...
                    'OCEAN must be TRUE or FALSE!');
            end
            option.OCEAN=val;
        case {'no' 'noo' 'nooc' 'noocean'}
            if(isempty(val)); continue; end
            if(~isscalar(val) || (~isreal(val) && ~islogical(val)))
                error([pkgme 'badInput'],...
                    'NOOCEAN must be TRUE or FALSE!');
            end
            option.OCEAN=~val;
        case {'c' 'cr' 'cru' 'crust'}
            if(isempty(val)); continue; end
            if(~isscalar(val) || (~isreal(val) && ~islogical(val)))
                error([pkgme 'badInput'],...
                    'CRUST must be TRUE or FALSE!');
            end
            option.CRUST=val;
        case {'nc' 'noc' 'nocr' 'nocru' 'nocrust'}
            if(isempty(val)); continue; end
            if(~isscalar(val) || (~isreal(val) && ~islogical(val)))
                error([pkgme 'badInput'],...
                    'NOCRUST must be TRUE or FALSE!');
            end
            option.CRUST=~val;
        case {'d' 'dep' 'depth' 'depths'}
            if(~isempty(val) && (~isreal(val) || any(val<0 | val>Re)))
                error([pkgme 'badInput'],...
                    'DEPTHS must be between 0-6371km!');
            end
            option.DEPTHS=val;
        case {'s' 'spvw' 'kpvw' 'k'}
            if(isempty(val)); continue; end
            if(~isreal(val) || ~isscalar(val) || val<=0)
                error([pkgme 'badInput'],...
                    'SPVW must be a positive real scalar!');
            end
            option.SPVW=val;
        case {'m' 'md' 'max' 'maxd' 'maxdep' 'maxdepth'}
            if(~isempty(val) && (~isreal(val) || ~isscalar(val)))
                error([pkgme 'badInput'],...
                    'MAXDEPTH must be a real-valued scalar in km!');
            end
            option.MAXDEPTH=val;
        otherwise
            error([pkgme 'unknownOption'],...
                'Unknown Option: %s !',varargin{i});
    end
end

% model name
mout.name='PREM2';

% model details
% is it anisotropic?
if(option.ANISOTROPY)
    % do we have a crust?
    if(option.CRUST)
        % do we have an ocean?
        if(option.OCEAN)
            mout.ocean=true;
            mout.crust=true;
            
            % the ocean (3km)
            model.order={'ocean'};
            model.ocean.radius=[6368 6371];
            model.ocean.depth=[0 3];
            model.ocean.rho=1.02;
            model.ocean.vpv=1.45;
            model.ocean.vsv=0;
            model.ocean.qk=57823;
            model.ocean.qu=inf;
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
            mout.ocean=false;
            mout.crust=true;
            
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
    else % NO CRUST
        mout.ocean=false;
        mout.crust=false;
        
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
    model.lowermantle.radius=[3840 5600];
    model.lowermantle.depth=[771 2531];
    model.lowermantle.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.lowermantle.vpv=[24.9520 -40.4673 51.4832 -26.6419];
    model.lowermantle.vsv=[11.1671 -13.7818 17.4575 -9.2777];
    model.lowermantle.qk=57823;
    model.lowermantle.qu=312;
    model.lowermantle.vph=[24.9520 -40.4673 51.4832 -26.6419];
    model.lowermantle.vsh=[11.1671 -13.7818 17.4575 -9.2777];
    model.lowermantle.eta=1;
    
    % upper d"
    model.order=[model.order {'upperddp'}];
    model.upperddp.radius=[3630 3840];
    model.upperddp.depth=[2531 2741];
    model.upperddp.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.upperddp.vpv=[14.2743 -1.3998];
    model.upperddp.vsv=[11.1671 -13.7818 17.4575 -9.2777];
    model.upperddp.qk=57823;
    model.upperddp.qu=312;
    model.upperddp.vph=[14.2743 -1.3998];
    model.upperddp.vsh=[11.1671 -13.7818 17.4575 -9.2777];
    model.upperddp.eta=1;
    
    % lower d"
    model.order=[model.order {'lowerddp'}];
    model.lowerddp.radius=[3480 3630];
    model.lowerddp.depth=[2741 2891];
    model.lowerddp.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.lowerddp.vpv=[14.2743 -1.3998];
    model.lowerddp.vsv=[6.9254 1.4672 -2.0834 0.9783];
    model.lowerddp.qk=57823;
    model.lowerddp.qu=312;
    model.lowerddp.vph=[14.2743 -1.3998];
    model.lowerddp.vsh=[6.9254 1.4672 -2.0834 0.9783];
    model.lowerddp.eta=1;
    
    % outer core
    model.order=[model.order {'outercore'}];
    model.outercore.radius=[1621.5 3480];
    model.outercore.depth=[2891 4749.5];
    model.outercore.rho=[12.5815 -1.2638 -3.6526 -5.5281];
    model.outercore.vpv=[11.0487 -4.0362 4.8023 -13.5732];
    model.outercore.vsv=0;
    model.outercore.qk=57823;
    model.outercore.qu=inf;
    model.outercore.vph=[11.0487 -4.0362 4.8023 -13.5732];
    model.outercore.vsh=0;
    model.outercore.eta=1;
    
    % lower outer core
    model.order=[model.order {'loweroutercore'}];
    model.loweroutercore.radius=[1221.5 1621.5];
    model.loweroutercore.depth=[4749.5 5149.5];
    model.loweroutercore.rho=[12.5815 -1.2638 -3.6526 -5.5281];
    model.loweroutercore.vpv=[4.0354 82.0080 -347.7690 468.7860];
    model.loweroutercore.vsv=0;
    model.loweroutercore.qk=57823;
    model.loweroutercore.qu=inf;
    model.loweroutercore.vph=[4.0354 82.0080 -347.7690 468.7860];
    model.loweroutercore.vsh=0;
    model.loweroutercore.eta=1;
    
    % upper inner core
    model.order=[model.order {'upperinnercore'}];
    model.upperinnercore.radius=[1010 1221.5];
    model.upperinnercore.depth=[5149.5 5361];
    model.upperinnercore.rho=[13.0885 0 -8.8381];
    model.upperinnercore.vpv=[11.3041 -1.2730];
    model.upperinnercore.vsv=[3.6678 0 -4.4475];
    model.upperinnercore.qk=1327.7;
    model.upperinnercore.qu=84.6;
    model.upperinnercore.vph=[11.3041 -1.2730];
    model.upperinnercore.vsh=[3.6678 0 -4.4475];
    model.upperinnercore.eta=1;
    
    % inner core
    model.order=[model.order {'innercore'}];
    model.innercore.radius=[0 1010];
    model.innercore.depth=[5361 6371];
    model.innercore.rho=[13.0885 0 -8.8381];
    model.innercore.vpv=[11.2622 0 -6.3640];
    model.innercore.vsv=[3.6678 0 -4.4475];
    model.innercore.qk=1327.7;
    model.innercore.qu=84.6;
    model.innercore.vph=[11.2622 0 -6.3640];
    model.innercore.vsh=[3.6678 0 -4.4475];
    model.innercore.eta=1;
    
    % number of layers
    no=numel(model.order);
    
    % velocity and density derivatives
    deriv={'rho' 'vpv' 'vsv' 'vph' 'vsh' 'eta'};
    for i=1:no
        for j=1:numel(deriv)
            nd=numel(model.(model.order{i}).(deriv{j}));
            if(nd>1)
                model.(model.order{i}).(['d' deriv{j} '_dr'])=...
                    model.(model.order{i}).(deriv{j})(2:end).*(1:nd-1);
            else
                model.(model.order{i}).(['d' deriv{j} '_dr'])=0;
            end
            if(nd>2)
                model.(model.order{i}).(['d2' deriv{j} '_dr2'])=...
                    model.(model.order{i}).(deriv{j})(3:end)...
                    .*(2:nd-1).*(1:nd-2);
            else
                model.(model.order{i}).(['d2' deriv{j} '_dr2'])=0;
            end
        end
    end
    
    % mass polynomials
    %
    %               re
    %              f
    %     m(r)  =  | rho(r)*4*pi*r^2*dr
    %              j
    %             0
    %
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
    
    % gravity polynomials
    % - start at r^-2 (ie must be divided by r^2)
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
    
    % pressure polynomials
    % - start at r^-2 (ie must be divided by r^2)
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
    
    % remaining info
    mout.isotropic=false;
    mout.refperiod=option.REFPERIOD;
    mout.flattened=false;
    
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
                        %mout.p(cnt)=mout.p(cnt)/r^2;
                    end
                end
            end
        end
        
        % truncate unused mout
        for i=1:nf
            mout.(fields{i})=mout.(fields{i})(1:cnt);
        end
        
        % correct velocities to be at REFPERIOD
        mout=phys_disp_ani(mout,option);
    else % SPVW at REFPERIOD
        % setup / preallocate output
        % 1100 is a guesstimate
        for i=1:nf
            mout.(fields{i})=nan(1100,1);
        end
        
        % start counter
        n=0;
        
        % loop over regions
        for i=1:no
            % start at top of region
            depth=model.(model.order{i}).depth(1);
            
            % finished if below max depth
            if(depth>option.MAXDEPTH); break; end
            n=n+1;
            
            % default polynomial method
            mout.depth(n)=depth;
            mout.radius(n)=Re-mout.depth(n);
            r=mout.radius(n)/Re;
            for k=3:nf
                mout.(fields{k})(n)=polyval(fliplr(...
                    model.(model.order{i}).(fields{k})),r);
            end
            if(r>0); mout.g(n)=mout.g(n)/r^2; end
            
            % correct velocities to be at REFPERIOD
            mout=phys_disp_ani(mout,option,n);
            
            % step to next depth using lowest vertical velocity
            if(mout.vsv(n)>0)
                vmin=mout.vsv(n);
            else
                vmin=mout.vpv(n);
            end
            wavelength=vmin*option.REFPERIOD;
            stepsize=wavelength/option.SPVW;
            depth=depth+stepsize;
            
            % step until at/below bottom of region (within a meter)
            while(depth<(model.(model.order{i}).depth(2)-0.001))
                % finished if below max depth
                if(depth>option.MAXDEPTH); break; end
                n=n+1;
                
                % default polynomial method
                mout.depth(n)=depth;
                mout.radius(n)=Re-mout.depth(n);
                r=mout.radius(n)/Re;
                for k=3:nf
                    mout.(fields{k})(n)=polyval(fliplr(...
                        model.(model.order{i}).(fields{k})),r);
                end
                if(r>0); mout.g(n)=mout.g(n)/r^2; end

                % correct velocities to be at REFPERIOD
                mout=phys_disp_ani(mout,option,n);
                
                % step to next depth using lowest vertical velocity
                if(mout.vsv(n)>0)
                    vmin=mout.vsv(n);
                else
                    vmin=mout.vpv(n);
                end
                wavelength=vmin*option.REFPERIOD;
                stepsize=wavelength/option.SPVW;
                depth=depth+stepsize;
            end
            
            % finish at the bottom of region
            depth=model.(model.order{i}).depth(2);
            
            % finished if below max depth
            if(depth>option.MAXDEPTH); break; end
            n=n+1;
            
            % default polynomial method
            mout.depth(n)=depth;
            mout.radius(n)=Re-mout.depth(n);
            r=mout.radius(n)/Re;
            for k=3:nf
                mout.(fields{k})(n)=polyval(fliplr(...
                    model.(model.order{i}).(fields{k})),r);
            end
            if(r>0); mout.g(n)=mout.g(n)/r^2; end
            
            % correct velocities to be at REFPERIOD
            mout=phys_disp_ani(mout,option,n);
        end
        
        % truncate unused mout
        for i=1:nf
            mout.(fields{i})=mout.(fields{i})(1:n);
        end
    end
    
    % bulk sound speed
    mout.vbv=sqrt(mout.vpv.^2-4/3*mout.vsv.^2);
    mout.vbh=sqrt(mout.vph.^2-4/3*mout.vsh.^2);
    mout.dvbv_dr=(mout.vpv.*mout.dvpv_dr ...
        -4/3*mout.vsv.*mout.dvsv_dr)./mout.vbv;
    mout.dvbh_dr=(mout.vph.*mout.dvph_dr ...
        -4/3*mout.vsh.*mout.dvsh_dr)./mout.vbh;
    mout.d2vbv_dr2=(mout.dvpv_dr.*mout.d2vpv_dr2 ...
        -4/3*mout.dvsv_dr.*mout.d2vsv_dr2)./mout.dvbv_dr;
    mout.d2vbh_dr2=(mout.dvph_dr.*mout.d2vph_dr2 ...
        -4/3*mout.dvsh_dr.*mout.d2vsh_dr2)./mout.dvbh_dr;
    
    % now get moduli
    mout.poissonv=0.5*...
        (mout.vpv.^2-2*mout.vsv.^2)./(mout.vpv.^2-mout.vsv.^2);
    mout.poissonh=0.5*...
        (mout.vph.^2-2*mout.vsh.^2)./(mout.vph.^2-mout.vsh.^2);
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
            mout.ocean=true;
            mout.crust=true;
            
            % the ocean (3km)
            model.order={'ocean'};
            model.ocean.radius=[6368 6371];
            model.ocean.depth=[0 3];
            model.ocean.rho=1.02;
            model.ocean.vp=1.45;
            model.ocean.vs=0;
            model.ocean.qk=57823;
            model.ocean.qu=inf;

            % the upper crust
            model.order=[model.order {'uppercrust'}];
            model.uppercrust.radius=[6356 6368];
            model.uppercrust.depth=[3 15];
            model.uppercrust.rho=2.6;
            model.uppercrust.vp=5.8;
            model.uppercrust.vs=3.2;
            model.uppercrust.qk=57823;
            model.uppercrust.qu=600;
        else % NO OCEAN
            mout.ocean=false;
            mout.crust=true;
            
            % the upper crust (to the surface)
            model.order={'uppercrust'};
            model.uppercrust.radius=[6356 6371];
            model.uppercrust.depth=[0 15];
            model.uppercrust.rho=2.6;
            model.uppercrust.vp=5.8;
            model.uppercrust.vs=3.2;
            model.uppercrust.qk=57823;
            model.uppercrust.qu=600;
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

        % the mantle lithosphere (lid)
        model.order=[model.order {'lithosphere'}];
        model.lithosphere.radius=[6291 6346.6];
        model.lithosphere.depth=[24.4 80];
        model.lithosphere.rho=[2.6910 0.6924];
        model.lithosphere.vp=[4.1875 3.9382];
        model.lithosphere.vs=[2.1519 2.3481];
        model.lithosphere.qk=57823;
        model.lithosphere.qu=600;
    else
        mout.ocean=false;
        mout.crust=false;
        
        % the mantle lithosphere (lid)
        model.order={'lithosphere'};
        model.lithosphere.radius=[6291 6371];
        model.lithosphere.depth=[0 80];
        model.lithosphere.rho=[2.6910 0.6924];
        model.lithosphere.vp=[4.1875 3.9382];
        model.lithosphere.vs=[2.1519 2.3481];
        model.lithosphere.qk=57823;
        model.lithosphere.qu=600;
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
    
    % the upper mantle
    model.order=[model.order {'uppermantle'}];
    model.uppermantle.radius=[5971 6151];
    model.uppermantle.depth=[220 400];
    model.uppermantle.rho=[7.1089 -3.8045];
    model.uppermantle.vp=[20.3926 -12.2569];
    model.uppermantle.vs=[8.9496 -4.4597];
    model.uppermantle.qk=57823;
    model.uppermantle.qu=143;
    
    % upper transition zone
    model.order=[model.order {'uppertz'}];
    model.uppertz.radius=[5771 5971];
    model.uppertz.depth=[400 600];
    model.uppertz.rho=[11.2494 -8.0298];
    model.uppertz.vp=[39.7027 -32.6166];
    model.uppertz.vs=[22.3512 -18.5856];
    model.uppertz.qk=57823;
    model.uppertz.qu=143;
    
    % middle transition zone
    model.order=[model.order {'middletz'}];
    model.middletz.radius=[5701 5771];
    model.middletz.depth=[600 670];
    model.middletz.rho=[5.3197 -1.4836];
    model.middletz.vp=[19.0957 -9.8672];
    model.middletz.vs=[9.9839 -4.9324];
    model.middletz.qk=57823;
    model.middletz.qu=143;
    
    % lower transition zone
    model.order=[model.order {'lowertz'}];
    model.lowertz.radius=[5600 5701];
    model.lowertz.depth=[670 771];
    model.lowertz.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.lowertz.vp=[29.2766 -23.6027 5.5242 -2.5514];
    model.lowertz.vs=[22.3459 -17.2473 -2.0834 0.9783];
    model.lowertz.qk=57823;
    model.lowertz.qu=312;
    
    % lower mantle
    model.order=[model.order {'lowermantle'}];
    model.lowermantle.radius=[3840 5600];
    model.lowermantle.depth=[771 2531];
    model.lowermantle.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.lowermantle.vp=[24.9520 -40.4673 51.4832 -26.6419];
    model.lowermantle.vs=[11.1671 -13.7818 17.4575 -9.2777];
    model.lowermantle.qk=57823;
    model.lowermantle.qu=312;
    
    % upper d"
    model.order=[model.order {'upperddp'}];
    model.upperddp.radius=[3630 3840];
    model.upperddp.depth=[2531 2741];
    model.upperddp.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.upperddp.vp=[14.2743 -1.3998];
    model.upperddp.vs=[11.1671 -13.7818 17.4575 -9.2777];
    model.upperddp.qk=57823;
    model.upperddp.qu=312;
    
    % lower d"
    model.order=[model.order {'lowerddp'}];
    model.lowerddp.radius=[3480 3630];
    model.lowerddp.depth=[2741 2891];
    model.lowerddp.rho=[7.9565 -6.4761 5.5283 -3.0807];
    model.lowerddp.vp=[14.2743 -1.3998];
    model.lowerddp.vs=[6.9254 1.4672 -2.0834 0.9783];
    model.lowerddp.qk=57823;
    model.lowerddp.qu=312;
    
    % outer core
    model.order=[model.order {'outercore'}];
    model.outercore.radius=[1621.5 3480];
    model.outercore.depth=[2891 4749.5];
    model.outercore.rho=[12.5815 -1.2638 -3.6526 -5.5281];
    model.outercore.vp=[11.0487 -4.0362 4.8023 -13.5732];
    model.outercore.vs=0;
    model.outercore.qk=57823;
    model.outercore.qu=inf;
    
    % lower outer core
    model.order=[model.order {'loweroutercore'}];
    model.loweroutercore.radius=[1221.5 1621.5];
    model.loweroutercore.depth=[4749.5 5149.5];
    model.loweroutercore.rho=[12.5815 -1.2638 -3.6526 -5.5281];
    model.loweroutercore.vp=[4.0354 82.0080 -347.7690 468.7860];
    model.loweroutercore.vs=0;
    model.loweroutercore.qk=57823;
    model.loweroutercore.qu=inf;
    
    % upper inner core
    model.order=[model.order {'upperinnercore'}];
    model.upperinnercore.radius=[1010 1221.5];
    model.upperinnercore.depth=[5149.5 5361];
    model.upperinnercore.rho=[13.0885 0 -8.8381];
    model.upperinnercore.vp=[11.3041 -1.2730];
    model.upperinnercore.vs=[3.6678 0 -4.4475];
    model.upperinnercore.qk=1327.7;
    model.upperinnercore.qu=84.6;
    
    % inner core
    model.order=[model.order {'innercore'}];
    model.innercore.radius=[0 1010];
    model.innercore.depth=[5361 6371];
    model.innercore.rho=[13.0885 0 -8.8381];
    model.innercore.vp=[11.2622 0 -6.3640];
    model.innercore.vs=[3.6678 0 -4.4475];
    model.innercore.qk=1327.7;
    model.innercore.qu=84.6;
    
    % number of layers
    no=numel(model.order);
    
    % velocity and density derivatives
    deriv={'rho' 'vp' 'vs'};
    for i=1:no
        for j=1:numel(deriv)
            nd=numel(model.(model.order{i}).(deriv{j}));
            if(nd>1)
                model.(model.order{i}).(['d' deriv{j} '_dr'])=...
                    model.(model.order{i}).(deriv{j})(2:end).*(1:nd-1);
            else
                model.(model.order{i}).(['d' deriv{j} '_dr'])=0;
            end
            if(nd>2)
                model.(model.order{i}).(['d2' deriv{j} '_dr2'])=...
                    model.(model.order{i}).(deriv{j})(3:end)...
                    .*(2:nd-1).*(1:nd-2);
            else
                model.(model.order{i}).(['d2' deriv{j} '_dr2'])=0;
            end
        end
    end
    
    % mass polynomials
    %
    %               re
    %              f
    %     m(r)  =  | p(r)*4*pi*r^2*dr
    %              j
    %             0
    %
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
    
    % gravity polynomials 
    % - start at r^-2 (ie must be divided by r^2)
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
    
    % remaining info
    mout.isotropic=true;
    mout.refperiod=option.REFPERIOD;
    mout.flattened=false;
    
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
            end
        end
        
        % truncate unused mout
        for i=1:nf
            mout.(fields{i})=mout.(fields{i})(1:cnt);
        end
        
        % correct velocities to be at REFPERIOD
        mout=phys_disp_iso(mout,option);
    else % SPVW at REFPERIOD
        % setup / preallocate output
        % 1100 is a guesstimate
        for i=1:nf
            mout.(fields{i})=nan(1100,1);
        end
        
        % start counter
        n=0;
        
        % loop over regions
        for i=1:no
            % start at top of region
            depth=model.(model.order{i}).depth(1);
            
            % finished if below max depth
            if(depth>option.MAXDEPTH); break; end
            n=n+1;
            
            % default polynomial method
            mout.depth(n)=depth;
            mout.radius(n)=Re-mout.depth(n);
            r=mout.radius(n)/Re;
            for k=3:nf
                mout.(fields{k})(n)=polyval(fliplr(...
                    model.(model.order{i}).(fields{k})),r);
            end
            if(r>0); mout.g(n)=mout.g(n)/r^2; end
            
            % correct velocities to be at REFPERIOD
            mout=phys_disp_iso(mout,option,n);
            
            % step to next depth using lowest vertical velocity
            if(mout.vs(n)>0)
                vmin=mout.vs(n);
            else
                vmin=mout.vp(n);
            end
            wavelength=vmin*option.REFPERIOD;
            stepsize=wavelength/option.SPVW;
            depth=depth+stepsize;
            
            % step until at/below bottom of region (within a meter)
            while(depth<(model.(model.order{i}).depth(2)-0.001))
                % finished if below max depth
                if(depth>option.MAXDEPTH); break; end
                n=n+1;
                
                % default polynomial method
                mout.depth(n)=depth;
                mout.radius(n)=Re-mout.depth(n);
                r=mout.radius(n)/Re;
                for k=3:nf
                    mout.(fields{k})(n)=polyval(fliplr(...
                        model.(model.order{i}).(fields{k})),r);
                end
                if(r>0); mout.g(n)=mout.g(n)/r^2; end

                % correct velocities to be at REFPERIOD
                mout=phys_disp_iso(mout,option,n);
                
                % step to next depth using lowest vertical velocity
                if(mout.vs(n)>0)
                    vmin=mout.vs(n);
                else
                    vmin=mout.vp(n);
                end
                wavelength=vmin*option.REFPERIOD;
                stepsize=wavelength/option.SPVW;
                depth=depth+stepsize;
            end
            
            % finish at the bottom of region
            depth=model.(model.order{i}).depth(2);
            
            % finished if below max depth
            if(depth>option.MAXDEPTH); break; end
            n=n+1;
            
            % default polynomial method
            mout.depth(n)=depth;
            mout.radius(n)=Re-mout.depth(n);
            r=mout.radius(n)/Re;
            for k=3:nf
                mout.(fields{k})(n)=polyval(fliplr(...
                    model.(model.order{i}).(fields{k})),r);
            end
            if(r>0); mout.g(n)=mout.g(n)/r^2; end
            
            % correct velocities to be at REFPERIOD
            mout=phys_disp_iso(mout,option,n);
        end
        
        % truncate unused mout
        for i=1:nf
            mout.(fields{i})=mout.(fields{i})(1:n);
        end
    end
    
    % bulk sound speed
    mout.vb=sqrt(mout.vp.^2-4/3*mout.vs.^2);
    mout.dvb_dr=(mout.vp.*mout.dvp_dr ...
        -4/3*mout.vs.*mout.dvs_dr)./mout.vb;
    mout.d2vb_dr2=(mout.dvp_dr.*mout.d2vp_dr2 ...
        -4/3*mout.dvs_dr.*mout.d2vs_dr2)./mout.dvb_dr;
    
    % now get moduli
    mout.poisson=0.5*(mout.vp.^2-2*mout.vs.^2)./(mout.vp.^2-mout.vs.^2);
    mout.shear=mout.vs.^2.*mout.rho*1e9;
    mout.bulk=mout.vp.^2.*mout.rho*1e9-4/3*mout.shear;
    mout.youngs=2*mout.shear.*(1+mout.poisson);
    mout.lambda=mout.bulk-2/3*mout.shear;
end

end


function [mout]=phys_disp_voigt(mout,option,idx)
%PHYS_DISP_VOIGT    Corrects for physical dispersion in aniso velocity
if(nargin<3 || isempty(idx)); idx=':'; end
vsvoigt=sqrt((2*mout.vsv(idx).^2+mout.vsh(idx).^2)./3);
vpvoigt=sqrt((mout.vpv(idx).^2+4*mout.vph(idx).^2)./5);
L=4/3*(vsvoigt./vpvoigt).^2;
Ap=1-log(option.REFPERIOD)/pi*...
    (L./mout.qu(idx)+(1-L)./mout.qk(idx));
As=1-log(option.REFPERIOD)/pi./mout.qu(idx);
mout.vpv(idx)=mout.vpv(idx).*Ap;
mout.vph(idx)=mout.vph(idx).*Ap;
mout.vsv(idx)=mout.vsv(idx).*As;
mout.vsh(idx)=mout.vsh(idx).*As;
mout.dvpv_dr(idx)=mout.dvpv_dr(idx).*Ap;
mout.dvph_dr(idx)=mout.dvph_dr(idx).*Ap;
mout.dvsv_dr(idx)=mout.dvsv_dr(idx).*As;
mout.dvsh_dr(idx)=mout.dvsh_dr(idx).*As;
mout.d2vpv_dr2(idx)=mout.d2vpv_dr2(idx).*Ap;
mout.d2vph_dr2(idx)=mout.d2vph_dr2(idx).*Ap;
mout.d2vsv_dr2(idx)=mout.d2vsv_dr2(idx).*As;
mout.d2vsh_dr2(idx)=mout.d2vsh_dr2(idx).*As;
end


function [mout]=phys_disp_ani(mout,option,idx)
%PHYS_DISP_ANI    Corrects for physical dispersion in anisotropic velocity
if(nargin<3 || isempty(idx)); idx=':'; end
Lv=4/3*(mout.vsv(idx)./mout.vpv(idx)).^2;
Lh=4/3*(mout.vsh(idx)./mout.vph(idx)).^2;
Apv=1-log(option.REFPERIOD)/pi*...
    (Lv./mout.qu(idx)+(1-Lv)./mout.qk(idx));
Aph=1-log(option.REFPERIOD)/pi*...
    (Lh./mout.qu(idx)+(1-Lh)./mout.qk(idx));
As=1-log(option.REFPERIOD)/pi./mout.qu(idx);
mout.vpv(idx)=mout.vpv(idx).*Apv;
mout.vph(idx)=mout.vph(idx).*Aph;
mout.vsv(idx)=mout.vsv(idx).*As;
mout.vsh(idx)=mout.vsh(idx).*As;
mout.dvpv_dr(idx)=mout.dvpv_dr(idx).*Apv;
mout.dvph_dr(idx)=mout.dvph_dr(idx).*Aph;
mout.dvsv_dr(idx)=mout.dvsv_dr(idx).*As;
mout.dvsh_dr(idx)=mout.dvsh_dr(idx).*As;
mout.d2vpv_dr2(idx)=mout.d2vpv_dr2(idx).*Apv;
mout.d2vph_dr2(idx)=mout.d2vph_dr2(idx).*Aph;
mout.d2vsv_dr2(idx)=mout.d2vsv_dr2(idx).*As;
mout.d2vsh_dr2(idx)=mout.d2vsh_dr2(idx).*As;
end


function [mout]=phys_disp_iso(mout,option,idx)
%PHYS_DISP_ISO    Corrects for physical dispersion in isotropic velocity
if(nargin<3 || isempty(idx)); idx=':'; end
L=4/3*(mout.vs(idx)./mout.vp(idx)).^2;
Ap=1-log(option.REFPERIOD)/pi*...
    (L./mout.qu(idx)+(1-L)./mout.qk(idx));
As=1-log(option.REFPERIOD)/pi./mout.qu(idx);
mout.vp(idx)=mout.vp(idx).*Ap;
mout.vs(idx)=mout.vs(idx).*As;
mout.dvp_dr(idx)=mout.dvp_dr(idx).*Ap;
mout.dvs_dr(idx)=mout.dvs_dr(idx).*As;
mout.d2vp_dr2(idx)=mout.d2vp_dr2(idx).*Ap;
mout.d2vs_dr2(idx)=mout.d2vs_dr2(idx).*As;
end

