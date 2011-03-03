function [varargout]=cmb_1dmodel_library(varargin)
%CMB_1DMODEL_LIBRARY    A library of CMB 1D Earth models
%
%    Usage:    models=cmb_1dmodel_library('modname1','modname2',...)
%              models=cmb_1dmodel_library(model,'modname1','modname2',...)
%              [name,pv]=cmb_1dmodel_library(...)
%
%    Description:
%     MODELS=CMB_1DMODEL_LIBRARY('MODNAME1','MODNAME2',...) returns a
%     struct array of 1D models corresponding to the model names given by
%     MODNAMEs.  These models are all PREM perturbations of the lower
%     mantle.  Some model names correspond to a set of models.  The list of
%     currently supported models is in the Notes section below.
%
%     MODELS=CMB_1DMODEL_LIBRARY(MODEL,'MODNAME1','MODNAME2',...) uses the
%     1D model struct MODEL to create the models given my 'MODNAME1' etc.
%     The core-mantle boundary for these models is assumed to correspond to
%     that of PREM, us user beware!
%
%     [NAMES,PV]=CMB_1DMODEL_LIBRARY(...) returns outputs
%     for creating the models given the by MODNAME inputs. These outputs
%     (NAMES & PV) can be passed to PERTURB_1DMODEL with a standard 1D
%     model (such as PREM) to create the model.  See the examples below to
%     learn how this can be done.  NAMES is a cell array of the model names
%     and PV is a cell array of cell arrays that contain the inputs
%     property/value pairs to be passed to PERTURB_1DMODEL.  These models
%     are all PREM perturbations of the lower mantle, so using another
%     model will likely produce different results -- look out for
%     differences in core-mantle boundary depth!!!
%
%    Notes:
%     Published models:
%      M1, SYLO, SYLOoc, SYL1, S175, SLHO, SLHE, SLHA, SWDK, SGLE, SGLE2,
%      SGHP, SGRD, LAM+S, LAM-S, SKNA1, SKNA2, SPAC, SPAC2, PWDK, P1, PPAC,
%      LAM+P, LAM-P, LAM+S, LAM-S
%
%     Published model sets:
%      L2000, L2400, L2600, L2700, L2750, LL2000, LL2400, LL2600,
%      LL2700, LL2750, LL2000', SYLO., SYLO_, SLHO_, SGHE, ALRG
%
%     Custom models:
%      M1ps, M1ps2, M1pq, M1/3p, PSKNA1, PKNA1, PSLHE, PLHE, PSYLO, PYLO,
%      ULVZ, CONV
%
%     Custom model sets:
%      GPS, GQu, GPSQu
%
%     Set of models from Raul Wong's PhD thesis (Chp. 3):
%      RAUL
%
%    Examples:
%     % Reproduce Figure 4 of Wysession et al 1998:
%     plot1dmodel([prem cmb_1dmodel_library('syl1','slhe','slha',...
%         'sylo','sghe','swdk','sgle','sghp')],'vs');
%     xlim([6.9 7.5]);
%     ylim([2200 2900]);
%     xlabel('Vs (km/sec)');
%     ylabel('Depth (km)');
%
%    See also: PREM, PERTURB_1DMODEL, PLOT1DMODEL

%     Version History:
%        May  30, 2010 - initial version
%        Aug.  8, 2010 - fixed case for GQu & GPSQu
%        Aug. 22, 2010 - added models from Raul's thesis
%        Sep. 19, 2010 - 1 or 2 outputs
%        Jan. 25, 2011 - allow startmodel argument, fixing some of the
%                        models to allow easy transfer to other parameters,
%                        use thk2n now, fixed a few cases of oversampling,
%                        added ulvz, PYLO
%        Jan. 30, 2011 - added GP & GS gradient series
%        Feb. 28, 2011 - added CONV model
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 28, 2011 at 23:00 GMT

% todo:
% - models should have named discontinuities field
%   ('moho' '410' '660' 'cmb' 'icb')
%   - this is tough to keep consistent with the model if we update it
% - what about named regions?
%   ('crust' 'lithosphere' 'asthenosphere' 'mantle' 'core' ...)
%
% - more p-wave models!
% - convection model series
% - ppv model series
% - ulvz series
% - two gradient series
% - gradient + discon + gradient series
% - lamellae series
% - 3 block series

% how many string inputs
nmods=0;
for i=1:nargin; if(isstring(varargin{i})); nmods=nmods+1; end; end

% get starting model (PREM!)
startmod=prem;
models(1:nmods)=startmod;

% get knot density
kpkm=numel(startmod.depth)/6371;

% loop over inputs
cnt=0; % number of output models
for i=1:nargin
    % check for start model change
    if(isstruct(varargin{i}))
        if(isempty(chk1dmodel(varargin{i})))
            % is a valid model
            startmod=varargin{i};
            continue;
        else
            % something is wrong
            error(chk1dmodel(varargin{i}));
        end
    end
    
    % require strings
    if(~isstring(varargin{i}))
        error('seizmo:cmb_1dmodel_library:badInput',...
            'All MODNAMEs must be a string!');
    end
    
    % get perturbation
    pv=cell(1,1);
    switch upper(varargin{i})
        % from Ritsema, Garnero, & Lay 1997
        % JGR Vol 102, No B9, Pg 20395-20411
        % "A strongly negative shear velocity gradient and lateral
        %  variability in the lowermost mantle beneath the Pacific"
        case 'M1'
            % their preferred model (aka LL2700_-3% & LL2700'_-3%)
            name={'M1'};
            pv{1}={'vs' [2000   0 4 0 thk2n(700,kpkm);
                         2700 -.5 4 0 thk2n(191,kpkm);
                         2891  -3 4 0 0]};
        case 'L2000'
            name={'L2000_-1%' 'L2000_-2%' 'L2000_-3%'};
            pv{1}={'vs' [2000   0 4 0 thk2n(891,kpkm); 2891  -1 4 0 0]};
            pv{2}={'vs' [2000   0 4 0 thk2n(891,kpkm); 2891  -2 4 0 0]};
            pv{3}={'vs' [2000   0 4 0 thk2n(891,kpkm); 2891  -3 4 0 0]};
        case 'L2400'
            name={'L2400_-1%' 'L2400_-2%' 'L2400_-3%'};
            pv{1}={'vs' [2400   0 4 0 thk2n(491,kpkm); 2891  -1 4 0 0]};
            pv{2}={'vs' [2400   0 4 0 thk2n(491,kpkm); 2891  -2 4 0 0]};
            pv{3}={'vs' [2400   0 4 0 thk2n(491,kpkm); 2891  -3 4 0 0]};
        case 'L2600'
            name={'L2600_-1%' 'L2600_-2%' 'L2600_-3%'};
            pv{1}={'vs' [2600   0 4 0 thk2n(291,kpkm); 2891  -1 4 0 0]};
            pv{2}={'vs' [2600   0 4 0 thk2n(291,kpkm); 2891  -2 4 0 0]};
            pv{3}={'vs' [2600   0 4 0 thk2n(291,kpkm); 2891  -3 4 0 0]};
        case 'L2700'
            name={'L2700_-1%' 'L2700_-2%' 'L2700_-3%'};
            pv{1}={'vs' [2700   0 4 0 thk2n(191,kpkm); 2891  -1 4 0 0]};
            pv{2}={'vs' [2700   0 4 0 thk2n(191,kpkm); 2891  -2 4 0 0]};
            pv{3}={'vs' [2700   0 4 0 thk2n(191,kpkm); 2891  -3 4 0 0]};
        case 'L2750'
            name={'L2750_-1%' 'L2750_-2%' 'L2750_-3%'};
            pv{1}={'vs' [2750   0 4 0 thk2n(141,kpkm); 2891  -1 4 0 0]};
            pv{2}={'vs' [2750   0 4 0 thk2n(141,kpkm); 2891  -2 4 0 0]};
            pv{3}={'vs' [2750   0 4 0 thk2n(141,kpkm); 2891  -3 4 0 0]};
        case 'LL2600'
            name={'LL2600_-1%' 'LL2600_-2%' 'LL2600_-3%'};
            pv{1}={'vs' [2000    0 4 0 thk2n(600,kpkm);
                         2600 -0.5 4 0 thk2n(291,kpkm);
                         2891   -3 4 0 0]};
            pv{2}={'vs' [2000    0 4 0 thk2n(600,kpkm);
                         2600 -1.0 4 0 thk2n(291,kpkm);
                         2891   -3 4 0 0]};
            pv{3}={'vs' [2000    0 4 0 thk2n(600,kpkm);
                         2600 -1.5 4 0 thk2n(291,kpkm);
                         2891   -3 4 0 0]};
        case 'LL2650'
            name={'LL2650_-1%' 'LL2650_-2%' 'LL2650_-3%'};
            pv{1}={'vs' [2000    0 4 0 thk2n(650,kpkm);
                         2650 -0.5 4 0 thk2n(241,kpkm);
                         2891   -3 4 0 0]};
            pv{2}={'vs' [2000    0 4 0 thk2n(650,kpkm);
                         2650 -1.0 4 0 thk2n(241,kpkm);
                         2891   -3 4 0 0]};
            pv{3}={'vs' [2000    0 4 0 thk2n(650,kpkm);
                         2650 -1.5 4 0 thk2n(241,kpkm);
                         2891   -3 4 0 0]};
        case 'LL2700'
            % note that LL2700_-3% is M1
            name={'LL2700_-1%' 'LL2700_-2%' 'LL2700_-3%'};
            pv{1}={'vs' [2000    0 4 0 thk2n(700,kpkm);
                         2700 -0.5 4 0 thk2n(191,kpkm);
                         2891   -3 4 0 0]};
            pv{2}={'vs' [2000    0 4 0 thk2n(700,kpkm);
                         2700 -1.0 4 0 thk2n(191,kpkm);
                         2891   -3 4 0 0]};
            pv{3}={'vs' [2000    0 4 0 thk2n(700,kpkm);
                         2700 -1.5 4 0 thk2n(191,kpkm);
                         2891   -3 4 0 0]};
        case 'LL2750'
            name={'LL2750_-1%' 'LL2750_-2%' 'LL2750_-3%'};
            pv{1}={'vs' [2000    0 4 0 thk2n(750,kpkm);
                         2750 -0.5 4 0 thk2n(141,kpkm);
                         2891   -3 4 0 0]};
            pv{2}={'vs' [2000    0 4 0 thk2n(750,kpkm);
                         2750 -1.0 4 0 thk2n(141,kpkm);
                         2891   -3 4 0 0]};
            pv{3}={'vs' [2000    0 4 0 thk2n(750,kpkm);
                         2750 -1.5 4 0 thk2n(141,kpkm);
                         2891   -3 4 0 0]};
        case 'LL2700'''
            % note that LL2700'_-3% is M1
            name={'LL2700''_-2%' 'LL2700''_-3%' 'LL2700''_-4%'};
            pv{1}={'vs' [2000   0 4 0 thk2n(700,kpkm);
                         2700 -.5 4 0 thk2n(191,kpkm);
                         2891  -2 4 0 0]};
            pv{2}={'vs' [2000   0 4 0 thk2n(700,kpkm);
                         2700 -.5 4 0 thk2n(191,kpkm);
                         2891  -3 4 0 0]};
            pv{3}={'vs' [2000   0 4 0 thk2n(700,kpkm);
                         2700 -.5 4 0 thk2n(191,kpkm);
                         2891  -4 4 0 0]};
        % from Lay & Young 1989
        % GRL Vol 16, No 7, Pg 605-608
        % "Waveform complexity in teleseismic broadband SH displacements:
        %  slab diffractions or deep mantle reflections?"
        %
        % Note that not all SYLO details have been published to
        % my knowledge.  We assume the drop from PREM above the
        % discontinuity is 1/4 the jump in velocity (this is based
        % on looking at Figure 9 from Garnero & Lay 1997).  The
        % constant velocity region just below the discontinuity
        % is assumed to be 43km thick because this gives a 200km
        % gradient layer and the velocity at the CMB is assumed to
        % be -0.085km/s below that at the top of D" (based on
        % matching SGLE's difference with PREM via Figure 3 of
        % Gaherty & Lay 1992)
        %
        % this needs work :(
        % - it appears all plots in papers indicate a jump of ~2.6%
        %   so we use that here and for all derivative models
        case 'SYLO'
            name={'SYLO'};
            sylo=100*((1/(1+.026/4))-1);
            pv{1}={'vs' [2350    0 4 0 thk2n(298,kpkm);
                         2648 sylo 4 0 0;
                         2648  2.6 2 0 thk2n(43,kpkm);
                         2691    0 2 0 thk2n(200,kpkm);
                         2891 -1.2 2 0 0]};
        case 'SYLO_'
            % SYLO but with 30-90 km width to the jump
            name={'SYLO_30' 'SYLO_60' 'SYLO_90'};
            dv30=100*((1/(1+.026/4))-1)*283/298;
            dv60=100*((1/(1+.026/4))-1)*268/298;
            dv90=100*((1/(1+.026/4))-1)*253/298;
            % these were found using
            %  sylo=cmb_1dmodel_library('sylo');
            %  vs=interpdc1(sylo.depth,sylo.vs,2663:15:2693);
            %vs30=7.3715359419987; % from 2.75%
            vs30=7.3635162692843; % 2.6% (consistent w/ SYLO figures)
            vs60=vs30;
            vs90=vs30-0.09/100;
            pv{1}={'vs' [2350    0 4 0 thk2n(283,kpkm);
                         2633 dv30 4 0 thk2n(30,kpkm);
                         2663 vs30 0 0 thk2n(28,kpkm);
                         2691    0 2 0 thk2n(200,kpkm);
                         2891 -1.2 2 0 0]};
            pv{2}={'vs' [2350    0 4 0 thk2n(268,kpkm);
                         2618 dv60 4 0 thk2n(60,kpkm);
                         2678 vs60 0 0 thk2n(13,kpkm);
                         2691    0 2 0 thk2n(200,kpkm);
                         2891 -1.2 2 0 0]};
            pv{3}={'vs' [2350    0 4 0 thk2n(253,kpkm);
                         2603 dv90 4 0 thk2n(90,kpkm);
                         2693 vs90 0 0 thk2n(198,kpkm);
                         2891 -1.2 2 0 0]};
        % from Young & Lay 1990
        % JGR Vol 95, Pg 17385-17402
        % "Multiple phase analysis of the shear velocity
        %  structure in the D" region beneath Alaska"
        %
        % This is just SYLO with a adjustment to the top of the outercore.
        case 'SYLOOC'
            name={'SYLOoc'};
            sylo=100*((1/(1+.026/4))-1);
            pv{1}={'vs' [2350    0 4 0 thk2n(298,kpkm);
                         2648 sylo 4 0 0;
                         2648  2.6 2 0 thk2n(43,kpkm);
                         2691    0 2 0 thk2n(200,kpkm);
                         2891 -1.2 2 0 0],...
                   'vp' [2891  7.9 0 0 thk2n(100,kpkm);
                         2991    0 4 0 0]};
        % from Garnero & Lay 1997
        % JGR Vol 102, No B4, Pg 8121-8135
        % "Lateral variations in lowermost mantle shear wave
        %  anisotropy beneath the north Pacific and Alaska"
        case 'SYLO.'
            % SYLO but with a decreased D" jump
            sylo=100*((1/(1+.026/4))-1);
            name={'SYLO.25' 'SYLO.50' 'SYLO.75'};
            pv{1}={'vs' [2350       0 4 0 thk2n(298,kpkm);
                         2648    sylo 4 0 0;
                         2648 1/4*2.6 2 0 thk2n(43,kpkm);
                         2691       0 2 0 thk2n(200,kpkm);
                         2891    -1.2 2 0 0]};
            pv{2}={'vs' [2350       0 4 0 thk2n(298,kpkm);
                         2648    sylo 4 0 0;
                         2648 1/2*2.6 2 0 thk2n(43,kpkm);
                         2691       0 2 0 thk2n(200,kpkm);
                         2891    -1.2 2 0 0]};
            pv{3}={'vs' [2350       0 4 0 thk2n(298,kpkm);
                         2648    sylo 4 0 0;
                         2648 3/4*2.6 2 0 thk2n(43,kpkm);
                         2691       0 2 0 thk2n(200,kpkm);
                         2891    -1.2 2 0 0]};
        case 'S175'
            name={'S175'};
            dv2716=100*((1/(1+.026/4))-1)*366/298;
            % value that the discontinuity jumps to was found
            % using:
            %  sylo=cmb_1dmodel_library('sylo');
            %  vs=interpdc1(sylo.depth,sylo.vs,2716);
            pv{1}={'vs' [2350        0 4 0 thk2n(366,kpkm);
                         2716   dv2716 4 0 0;
                         2716     2.16 2 0 thk2n(175,kpkm);
                         2891 -.09*7/8 1 0 0]};
        % from Young & Lay 1987
        % PEPI Vol 49, Pg 37-53
        % "Evidence for a shear velocity discontinuity in the
        % lower mantle beneath India and the Indian Ocean"
        %
        % the tradition of not divulging the model continues
        case 'SYL1'
            name={'SYL1'};
            pv{1}={'vs' [2300     0 4 0 thk2n(309,kpkm);
                         2609 7.125 0 0 0;
                         2609   3.1 2 0 thk2n(282,kpkm);
                         2891 7.195 0 0 0]};
        % from Lay & Helmberger 1983
        % GJRAS Vol 75, Pg 799-837
        % "A lower mantle S-wave triplication
        %  and the velocity structure of D""
        %
        % These are all based on the figures
        % & info from Wysession et al 1998
        %
        % note that these had a core radius of 3485km
        % prem radius 3480km
        %   jb radius 3473km
        % so the discontinuities are 5km deeper with PREM!
        % also in Wysession et al 1998 the discontinuities are 2km higher
        %
        % slho 278 2.75
        % slhe 318 2.75
        % slha 251 2.75
        case 'SLHO'
            name={'SLHO'};
            pv{1}={'vs' [2300              0 4 0 thk2n(311,kpkm);
                         2611   311/(6371*3) 1 0 0;
                         2611           2.75 2 0 thk2n(280,kpkm);
                         2891 280/(6371*1.1) 1 0 0]};
        case 'SLHE'
            name={'SLHE'};
            pv{1}={'vs' [2300     0 4 0 thk2n(271,kpkm);
                         2571   0.2 2 0 0;
                         2571  2.75 2 0 thk2n(320,kpkm);
                         2891 0.575 2 0 0]};
        case 'SLHA'
            name={'SLHA'};
            pv{1}={'vs' [2350              0 4 0 thk2n(288,kpkm);
                         2638   288/(6371*3) 1 0 0;
                         2638           2.75 2 0 thk2n(253,kpkm);
                         2891 253/(6371*1.9) 1 0 0]};
        case 'SLHO_'
            % 18km, 30km, 45km wide discon
            % values below discon found using:
            %  slho=cmb_1dmodel_library('slho');
            %  vs=interpdc1(slho.depth,slho.vs,[2620 2626 2633.5]);
            name={'SLHO_18' 'SLHO_30' 'SLHO_45'};
            pv{1}={'vs' [2300              0 4 0 thk2n(302,kpkm);
                         2602   302/(6371*3) 1 0 thk2n(18,kpkm);
                         2620         7.2805 0 0 thk2n(271,kpkm);
                         2891 271/(6371*1.1) 1 0 0]};
            pv{2}={'vs' [2300              0 4 0 thk2n(296,kpkm);
                         2596   296/(6371*3) 1 0 thk2n(30,kpkm);
                         2626         7.2814 0 0 thk2n(265,kpkm);
                         2891 265/(6371*1.1) 1 0 0]};
            pv{3}={'vs' [2300                  0 4 0 thk2n(288.5,kpkm);
                         2588.5   288.5/(6371*3) 1 0 thk2n(45,kpkm);
                         2633.5           7.2825 0 0 thk2n(257.5,kpkm);
                         2891   257.5/(6371*1.1) 1 0 0]};
        % from Garnero et al 1988
        % GRL Vol 15, No 6, Pg 609-612
        % "Lateral variations near the core-mantle boundary"
        %
        % - note that this is set relative to PREM
        % - used 253km as SGHE appears to be the same as SLHA
        case 'SGHE'
            name={'SGHE_A' 'SGHE_B'};
            pv{1}={'vs' [2350            0 4 0 thk2n(288,kpkm);
                         2638 288/(6371*3) 1 0 thk2n(253,kpkm);
                         2891          7.2 0 0 0]};
            pv{2}={'vs' [2350            0 4 0 thk2n(288,kpkm);
                         2638 288/(6371*3) 1 0 0;
                         2638         2.75 2 0 thk2n(253,kpkm);
                         2891          7.2 0 0 0]};
        % from Weber and Davis 1990
        % GJI Vol 102, Pg 231-255
        % "Evidence of a laterally variable lower
        %  mantle structure from p- and s-waves"
        case 'SWDK'
            name={'SWDK'};
            pv{1}={'vs' [2300              0 4 0 thk2n(305,kpkm);
                         2605   305/(6371*3) 1 0 0;
                         2605              3 2 0 thk2n(286,kpkm);
                         2891 286/(6371*1.3) 1 0 0]};
        case 'PWDK'
            name={'PWDK'};
            pv{1}={'vp' [2300             0 4 0 thk2n(305,kpkm);
                         2605 305/(6371*.5) 1 0 0;
                         2605             3 2 0 thk2n(286,kpkm);
                         2891          13.7 0 0 0]};
        % from Gaherty and Lay 1992
        % JGR Vol 97, No B1, Pg 417-435
        % "Investigation of laterally heterogeneous shear
        %  velocity structure in D" beneath Eurasia"
        %
        % SYLO-like model: starts at 2300 rather than 2350, but with same
        %                  gradient, 2.75% jump at 285km above CMB, flat
        %                  until 200km above CMB, gradient to same as SYLO
        %                  at CMB
        case 'SGLE'
            name={'SGLE'};
            sylo=100*((1/(1+.026/4))-1)*306/298;
            pv{1}={'vs' [2300      0 4 0 thk2n(306,kpkm);
                         2606   sylo 4 0 0;
                         2606   2.75 2 0 thk2n(85,kpkm);
                         2691      0 1 0 thk2n(200,kpkm);
                         2891 7.2771 0 0 0]};
        case 'SGLE2'
            name={'SGLE2'};
            pv{1}={'vs' [2300    0 4 0 thk2n(291,kpkm)
                         2591  7.1 0 0 0;
                         2591  2.6 2 0 thk2n(140,kpkm);
                         2731 7.23 0 0 0;
                         2731  2.6 2 0 thk2n(160,kpkm);
                         2891 7.34 0 0 0]};
        % from Garnero et al 1993
        % PEPI Vol 79, Pg 335-347
        % "Preliminary evidence for a lower mantle shear wave
        %  velocity discontinuity beneath the central Pacific"
        %
        % SGHP based on plots.  Discon 176km above CMB, 2.55% jump,
        % curved gradient under discon, back near prem at CMB.
        %
        % SGRD based on plots.  Not great but meh.
        case 'SGHP'
            name={'SGHP'};
            pv{1}={'vs' [2300             0 4 0 thk2n(415,kpkm);
                         2715 415/(6371*.7) 1 0 0;
                         2715          2.55 2 2 thk2n(65,kpkm);
                         2780         -.015 1 0 thk2n(111,kpkm);
                         2891        7.2666 0 0 0]};
        case 'SGRD'
            name={'SGRD'};
            pv{1}={'vs' [2300             0 4    0 thk2n(320,kpkm);
                         2620 320/(6371*.7) 1 1.25 thk2n(40,kpkm);
                         2660        7.1666 0    0 thk2n(80,kpkm);
                         2740           7.3 0   -2 thk2n(25,kpkm);
                         2765         7.325 0    0 thk2n(126,kpkm);
                         2891        7.2666 0    0 0]};
        % from Weber 1994
        % GRL Vol 21, No 23, Pg 2531-2534
        % "Lamellae in D"? An alternative model for lower mantle anomalies"
        %
        % Single lamellae centered on 2700km depth, 3% jump
        % 
        case 'LAM+P'
            name={'LAM+P'};
            lamth=27;
            lamdv=3;
            pv{1}={'vp' [2700-lamth/2     0 4 0 0;
                         2700-lamth/2 lamdv 2 0 thk2n(lamth,kpkm);
                         2700+lamth/2     0 2 0 0;
                         2700+lamth/2     0 4 0 0]};
        case 'LAM-P'
            name={'LAM-P'};
            lamth=27;
            lamdv=-3;
            pv{1}={'vp' [2700-lamth/2     0 4 0 0;
                         2700-lamth/2 lamdv 2 0 thk2n(lamth,kpkm);
                         2700+lamth/2     0 2 0 0;
                         2700+lamth/2     0 4 0 0]};
        case 'LAM+S'
            name={'LAM+S'};
            lamth=27;
            lamdv=3;
            pv{1}={'vs' [2700-lamth/2     0 4 0 0;
                         2700-lamth/2 lamdv 2 0 thk2n(lamth,kpkm);
                         2700+lamth/2     0 2 0 0;
                         2700+lamth/2     0 4 0 0]};
        case 'LAM-S'
            name={'LAM-S'};
            lamth=27;
            lamdv=-3;
            pv{1}={'vs' [2700-lamth/2     0 4 0 0;
                         2700-lamth/2 lamdv 2 0 thk2n(lamth,kpkm);
                         2700+lamth/2     0 2 0 0;
                         2700+lamth/2     0 4 0 0]};
        % from Kendall and Nangini 1996
        % GRL Vol 23, No 4, Pg 399-402
        % "Lateral variations in D" below the Caribbean"
        %
        % 1 : 250km (2641), 2.75, 10km transition
        % 2 : 290km (2601), 2.45, ''
        %     the gradients are the same above
        %
        % From their figure it looks like skna2 is 300km above
        % and skna1 goes to 7.4km/s (~2.90 jump).  This is what is below.
        case 'SKNA1'
            name={'SKNA1'};
            pv{1}={'vs' [2200    0 4 0 thk2n(436,kpkm);
                         2636 7.19 0 0 2;
                         2646 2.90 2 0 thk2n(245,kpkm);
                         2891 7.28 0 0 0]};
        case 'SKNA2'
            name={'SKNA2'};
            pv{1}={'vs' [2200     0 4 0 thk2n(386,kpkm);
                         2586 7.174 0 0 2;
                         2596  7.35 0 0 thk2n(295,kpkm);
                         2891  7.35 0 0 0]};
        % from Russell et al 2004
        % GRL Vol 28, No 11, Pg 2281-2284
        % "Coexisting shear- and compressional-wave seismic
        %  velocity discontinuities beneath the central Pacific"
        case 'SPAC'
            name={'SPAC'};
            pv{1}={'vs' [2200    0 4 0 thk2n(461,kpkm);
                         2661 7.15 0 0 0;
                         2661  1.7 2 0 thk2n(220,kpkm);
                         2881 7.15 0 0 0;
                         2881  -12 2 0 5;
                         2891    0 2 0 0]};
        case 'PPAC'
            name={'PPAC'};
            pv{1}={'vp' [2100      0 4 0 thk2n(561,kpkm);
                         2661 13.483 0 0 0;
                         2661   0.75 2 0 thk2n(220,kpkm);
                         2881   13.8 0 0 0;
                         2881     -4 2 0 5;
                         2891      0 2 0 0]};
        case 'P1'
            name={'P1'};
            pv{1}={'vp' [2000     0 4 0 thk2n(700,kpkm);
                         2700 13.57 0 0 thk2n(191,kpkm);
                         2891  -.01 3 0 0]};
        % from Avants et al 2004
        % JGR Vol 111, B05305, doi:10.1029/2004JB003270
        % "Shear velocity variation within the D"
        %  region beneath the central Pacific"
        case 'ALRG'
            name={'BIN_A' 'BIN_B' 'BIN_C' 'BIN_D' 'BIN_E' ...
                  'BIN_F' 'BIN_G' 'BIN_H' 'BIN_I'};
            pv{1}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2702 7.164 0 0 0;
                         2702 7.286 0 0 thk2n(89,kpkm);
                         2791     0 2 0 0;
                         2791 7.199 0 0 thk2n(71,kpkm);
                         2862     0 2 0 0;
                         2862 6.875 0 0 thk2n(29,kpkm);
                         2891     0 2 0 0]};
            pv{2}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2568 7.127 0 0 0;
                         2568 7.234 0 0 thk2n(292,kpkm);
                         2860     0 2 0 0;
                         2860 7.096 0 0 thk2n(31,kpkm);
                         2891     0 2 0 0]};
            pv{3}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2719 7.169 0 0 0;
                         2719 7.226 0 0 thk2n(92,kpkm);
                         2811     0 2 0 0;
                         2811 7.104 0 0 thk2n(41,kpkm);
                         2852     0 2 0 0;
                         2852 6.962 0 0 thk2n(39,kpkm);
                         2891     0 2 0 0]};
            pv{4}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2552 7.122 0 0 0;
                         2552 7.236 0 0 thk2n(278,kpkm);
                         2830     0 2 0 0;
                         2830 7.215 0 0 thk2n(44,kpkm);
                         2874     0 2 0 0;
                         2874 6.940 0 0 thk2n(17,kpkm);
                         2891     0 2 0 0]};
            pv{5}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2652 7.150 0 0 0;
                         2652 7.200 0 0 thk2n(151,kpkm);
                         2803     0 2 0 0;
                         2803 7.107 0 0 thk2n(59,kpkm);
                         2862     0 2 0 0;
                         2862 6.872 0 0 thk2n(29,kpkm);
                         2891     0 2 0 0]};
            pv{6}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2735 7.174 0 0 0;
                         2735 7.238 0 0 thk2n(81,kpkm);
                         2816     0 2 0 0;
                         2816 7.072 0 0 thk2n(75,kpkm);
                         2891     0 2 0 0]};
            pv{7}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2490 7.105 0 0 0;
                         2490 7.268 0 0 thk2n(307,kpkm);
                         2797     0 2 0 0;
                         2797 7.217 0 0 thk2n(66,kpkm);
                         2863     0 2 0 0;
                         2863 7.116 0 0 thk2n(28,kpkm);
                         2891     0 2 0 0]};
            pv{8}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2692 7.162 0 0 0;
                         2692 7.197 0 0 thk2n(116,kpkm);
                         2808     0 2 0 0;
                         2808 7.089 0 0 thk2n(59,kpkm);
                         2867     0 2 0 0;
                         2867 6.898 0 0 thk2n(24,kpkm);
                         2891     0 2 0 0]};
            pv{9}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2710 7.167 0 0 0;
                         2710 7.202 0 0 thk2n(134,kpkm);
                         2844     0 2 0 0;
                         2844 6.684 0 0 thk2n(36,kpkm);
                         2880     0 2 0 0;
                         2880 5.013 0 0 thk2n(11,kpkm);
                         2891     0 2 0 0]};
        case 'SPAC2'
            name={'SPAC2'};
            pv{1}={'vs' [2200     0 4 0 thk2n(432,kpkm);
                         2632 7.145 0 0 0;
                         2632 7.202 0 0 thk2n(180,kpkm);
                         2812     0 2 0 0;
                         2812 7.123 0 0 thk2n(45,kpkm);
                         2857     0 2 0 0;
                         2857 6.937 0 0 thk2n(34,kpkm);
                         2891     0 2 0 0]};
        % Michael's Custom Models
        % M1p == M1 but with P too
        case 'M1PS'
            name='M1ps';
            pv{1}={'vs' [2000   0 4 0 thk2n(700,kpkm);
                         2700 -.5 4 0 thk2n(191,kpkm);
                         2891  -3 4 0 0],...
                   'vp' [2000   0 4 0 thk2n(700,kpkm);
                         2700 -.5 4 0 thk2n(191,kpkm);
                         2891  -3 4 0 0]};
        % M1p2 == M1p but with a x^2 like gradient
        case 'M1PS2'
            name='M1ps2';
            pv{1}={'vs' [2000   0 4   0 thk2n(700,kpkm);
                         2700 -.5 4 1.7 thk2n(191,kpkm);
                         2891  -3 4   0 0],...
                   'vp' [2000   0 4   0 thk2n(700,kpkm);
                         2700 -.5 4 1.7 thk2n(191,kpkm);
                         2891  -3 4   0 0]};
        % M1pq == M1p but with a strong negative gradient in Qk & Qu too
        case 'M1PSQ'
            name='M1psq';
            pv{1}={'vs' [2000   0 4 0 thk2n(700,kpkm);
                         2700 -.5 4 0 thk2n(191,kpkm);
                         2891  -3 4 0 0],...
                   'vp' [2000   0 4 0 thk2n(700,kpkm);
                         2700 -.5 4 0 thk2n(191,kpkm);
                         2891  -3 4 0 0],...
                   'qk' [2700   0 4 0 thk2n(191,kpkm);
                         2891 100 0 0 0],...
                   'qu' [2700   0 4 0 thk2n(191,kpkm);
                         2891  10 0 0 0]};
        % 1/3 deviations of M1p
        case 'M1/3PS'
            name='M1/3ps';
            pv{1}={'vs' [2000     0 4 0 thk2n(700,kpkm);
                         2700 -.5/3 4 0 thk2n(191,kpkm);
                         2891    -1 4 0 0],...
                   'vp' [2000     0 4 0 thk2n(700,kpkm);
                         2700 -.5/3 4 0 thk2n(191,kpkm);
                         2891    -1 4 0 0]};
        % skna1 but with similar perturbations to p
        case 'PSKNA1'
            name='PSKNA1';
            pv{1}={'vs' [2200    0 4 0 thk2n(436,kpkm);
                         2636 -.39 4 0 2;
                         2646 2.90 2 0 thk2n(245,kpkm);
                         2891 .212 4 0 0],...
                   'vp' [2200    0 4 0 thk2n(436,kpkm);
                         2636 -.39 4 0 2;
                         2646 2.90 2 0 thk2n(245,kpkm);
                         2891 .212 4 0 0]};
        % skna1 perturbations but only in p
        case 'PKNA1'
            name='PKNA1';
            pv{1}={'vp' [2200    0 4 0 thk2n(436,kpkm);
                         2636 -.39 4 0 2;
                         2646 2.90 2 0 thk2n(245,kpkm);
                         2891 .212 4 0 0]};
        % slhe but with similar perturbations to p
        case 'PSLHE'
            name={'PSLHE'};
            pv{1}={'vp' [2300     0 4 0 thk2n(271,kpkm);
                         2571   0.2 2 0 0;
                         2571  2.75 2 0 thk2n(320,kpkm);
                         2891 0.575 2 0 0],...
                   'vs' [2300     0 4 0 thk2n(271,kpkm);
                         2571   0.2 2 0 0;
                         2571  2.75 2 0 thk2n(320,kpkm);
                         2891 0.575 2 0 0]};
        % slhe but with similar perturbations to p
        % - triple knots in perturbed region
        case 'PSLHE3'
            name={'PSLHE'};
            pv{1}={'vp' [2300     0 4 0 thk2n(271,3*kpkm);
                         2571   0.2 2 0 0;
                         2571  2.75 2 0 thk2n(320,3*kpkm);
                         2891 0.575 2 0 0],...
                   'vs' [2300     0 4 0 thk2n(271,3*kpkm);
                         2571   0.2 2 0 0;
                         2571  2.75 2 0 thk2n(320,3*kpkm);
                         2891 0.575 2 0 0]};        
        % slhe perturbations but only in p
        case 'PLHE'
            name={'PLHE'};
            pv{1}={'vp' [2300     0 4 0 thk2n(271,kpkm);
                         2571   0.2 2 0 0;
                         2571  2.75 2 0 thk2n(320,kpkm);
                         2891 0.575 2 0 0]};
        % ulvz (dvp=-10,dvs=-30,thk=20km)
        case 'ULVZ'
            name='ULVZ';
            pv{1}={'vp' [2871   0 4 0 0;
                         2871 -10 4 0 thk2n(20,kpkm);
                         2891 -10 4 0 0],...
                   'vs' [2871   0 4 0 0;
                         2871 -30 4 0 thk2n(20,kpkm);
                         2891 -30 4 0 0]};
        % CONV (d" convection)
        % - based on Solomatov & Moresi 2002
        case 'CONV'
            name='CONV';
            pv{1}={'vp' [2300     0 4   -2 thk2n(341,kpkm*2);
                         2641    -2 4    0 0;
                         2641     1 4 1.25 thk2n(250,kpkm*2);
                         2891  1.25 4    0 0],...
                   'vs' [2300     0 4   -2 thk2n(341,kpkm*2);
                         2641    -2 4    0 0;
                         2641     1 4 1.25 thk2n(250,kpkm*2);
                         2891  1.25 4    0 0]};
        % SYLO with equal % P perturbations
        case 'PSYLO'
            name={'PSYLO'};
            pylo=100*((1/(1+.026/4))-1); % 1/4 of the 2.6% jump
            pv{1}={'vs' [2350    0 4 0 thk2n(298,kpkm);
                         2648 pylo 4 0 0;
                         2648  2.6 2 0 thk2n(43,kpkm);
                         2691    0 2 0 thk2n(200,kpkm);
                         2891 -1.2 2 0 0],...
                   'vp' [2350    0 4 0 thk2n(298,kpkm);
                         2648 pylo 4 0 0;
                         2648  2.6 2 0 thk2n(43,kpkm);
                         2691    0 2 0 thk2n(200,kpkm);
                         2891 -1.2 2 0 0]};
        % P equivalent of SYLO
        case 'PYLO'
            name={'PYLO'};
            pylo=100*((1/(1+.026/4))-1); % 1/4 of the 2.6% jump
            pv{1}={'vs' [2350    0 4 0 thk2n(298,kpkm);
                         2648 pylo 4 0 0;
                         2648  2.6 2 0 thk2n(43,kpkm);
                         2691    0 2 0 thk2n(200,kpkm);
                         2891 -1.2 2 0 0]};
        % gradient series (the L-series of Ritsema et al 1997)
        case 'GS'
            hng=[2000 2400 2600 2650 2700 2750];
            dv=-4:4;
            nhng=numel(hng);
            ndv=numel(dv);
            name=cell(nhng*ndv,1);
            pv=cell(nhng*ndv,1);
            for a=1:nhng;
                for b=1:ndv;
                    name{(a-1)*ndv+b}=['GS' num2str(hng(a)) '_' ...
                        num2str(dv(b)) '%'];
                    pv{(a-1)*ndv+b}=...
                        {'vs' [hng(a) 0 4 0 thk2n(2891-hng(a),kpkm);
                               2891 dv(b) 4 0 0]};
                end
            end
        % gradient series (extension of Ritsema et al 1997)
        case 'GP'
            hng=[2000 2400 2600 2650 2700 2750];
            dv=-4:4;
            nhng=numel(hng);
            ndv=numel(dv);
            name=cell(nhng*ndv,1);
            pv=cell(nhng*ndv,1);
            for a=1:nhng;
                for b=1:ndv;
                    name{(a-1)*ndv+b}=['GP' num2str(hng(a)) '_' ...
                        num2str(dv(b)) '%'];
                    pv{(a-1)*ndv+b}=...
                        {'vp' [hng(a) 0 4 0 thk2n(2891-hng(a),kpkm);
                               2891 dv(b) 4 0 0]};
                end
            end
        % gradient series (extension of Ritsema et al 1997)
        case 'GPS'
            hng=[2000 2400 2600 2650 2700 2750];
            dv=-4:4;
            nhng=numel(hng);
            ndv=numel(dv);
            name=cell(nhng*ndv,1);
            pv=cell(nhng*ndv,1);
            for a=1:nhng;
                for b=1:ndv;
                    name{(a-1)*ndv+b}=['GPS' num2str(hng(a)) '_' ...
                        num2str(dv(b)) '%'];
                    pv{(a-1)*ndv+b}=...
                        {'vs' [hng(a) 0 4 0 thk2n(2891-hng(a),kpkm);
                               2891 dv(b) 4 0 0],...
                         'vp' [hng(a) 0 4 0 thk2n(2891-hng(a),kpkm);
                               2891 dv(b) 4 0 0]};
                end
            end
        case 'GQU'
            hng=[2000 2400 2600 2650 2700 2750];
            dv=-80:20:80;
            nhng=numel(hng);
            ndv=numel(dv);
            name=cell(nhng*ndv,1);
            pv=cell(nhng*ndv,1);
            for a=1:nhng;
                for b=1:ndv;
                    name{(a-1)*ndv+b}=['GQu' num2str(hng(a)) '_' ...
                        num2str(dv(b)) '%'];
                    pv{(a-1)*ndv+b}=...
                        {'qu' [hng(a) 0 4 0 thk2n(2891-hng(a),kpkm);
                               2891 dv(b) 4 0 0]};
                end
            end
        case 'GPSQU'
            hng=[2000 2400 2600 2650 2700 2750];
            dv=-4:4;
            nhng=numel(hng);
            ndv=numel(dv);
            name=cell(nhng*ndv,1);
            pv=cell(nhng*ndv,1);
            for a=1:nhng;
                for b=1:ndv;
                    name{(a-1)*ndv+b}=['GPSQu' num2str(hng(a)) '_' ...
                        num2str(dv(b)) '%V_' num2str(100*dv(b)/5) '%Qu'];
                    pv{(a-1)*ndv+b}=...
                        {'vs' [hng(a) 0 4 0 thk2n(2891-hng(a),kpkm);
                               2891 dv(b) 4 0 0],...
                         'vp' [hng(a) 0 4 0 thk2n(2891-hng(a),kpkm);
                               2891 dv(b) 4 0 0],...
                         'qu' [hng(a) 0 4 0 thk2n(2891-hng(a),kpkm);
                               2891 100*dv(b)/5 4 0 0]};
                end
            end
        case 'RAUL'
            % Models from Raul Valenzuela Wong's PhD Thesis (1996) Chp. 3
            % Note: Density & Vp is NOT adjusted from PREM.
            %  A - decrease from PREM d" to 6.95
            %  B - increase from PREM d" to 7.35
            %  C - step at PREM d" to 7.3
            %  D - 2605, 2.95% to 7.19 linear
            %  E - 2605, 2.95% to 7.19 parabolic
            %  F02 - 2605, 2.95% to 7.19 error function
            %  G02 - 2605, 2.95% to 7.0 error function
            %  H - 2605, 2.95% to 6.9 error function
            %  I2505/2555/2605/2655/2705 (2605 is G02) - note values above and
            %      below the discontinuity are the same so they are not % of PREM
            %  J724/J734/J744 (J734 is G02) - values are at the bottom of the
            %      discontinuity (7.13 is the value at the top)
            %  K2655/2705 - like I models w/ jump to 7.44
            %  KE2655/2705 - like I models w/ jump to 7.4
            %  N690/695/700/706/710/713 - K2705 variations
            %  QE 39 7228 in PREM d" - N706 with various Qu
            %  QH 78 14456
            %  QK 156 28912
            %  QQ 624 99999.99
            %  QT 1248 99999.99
            %  L2705 - gradient version of K2705 (gradient from ~2630 to ~2730)
            %  M2705 - clipped version of K2705 (~7.34 at d", linear to 2790)
            name={'A' 'B' 'C' 'D' 'E' 'F02' 'G02' 'H' 'I2505' 'I2555' ...
                'I2655' 'I2705' 'J724' 'J744' 'K2655' 'K2705' 'KE2655' ...
                'KE2705' 'N690' 'N695' 'N700' 'N706' 'N710' 'QE' 'QH' ...
                'QK' 'QQ' 'QT' 'L2705' 'M2705'};
            pv={{'vs' [2741 0 4 0 thk2n(150,kpkm); 2891 6.95 0 0 0]} ...
                {'vs' [2741 0 4 0 thk2n(150,kpkm); 2891 7.35 0 0 0]} ...
                {'vs' [2741 0 4 0 0; 2741 7.3 0 0 thk2n(150,kpkm); 2891 7.3 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(305,kpkm); 2605 7.13 0 0 0; 2605 7.34 0 0 thk2n(286,kpkm); 2891 7.19 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(305,kpkm); 2605 7.13 0 0 0; 2605 7.34 0 2 thk2n(286,kpkm); 2891 7.19 0 0 0]} ... % this one needs to be fixed
                {'vs' [2300 0 4 0 thk2n(305,kpkm); 2605 7.13 0 0 0; 2605 7.34 0 2 thk2n(286,kpkm); 2891 7.19 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(305,kpkm); 2605 7.13 0 0 0; 2605 7.34 0 2 thk2n(286,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(305,kpkm); 2605 7.13 0 0 0; 2605 7.34 0 2 thk2n(286,kpkm); 2891 6.9 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(205,kpkm); 2505 7.13 0 0 0; 2505 7.34 0 2 thk2n(386,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(255,kpkm); 2555 7.13 0 0 0; 2555 7.34 0 2 thk2n(336,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(355,kpkm); 2655 7.13 0 0 0; 2655 7.34 0 2 thk2n(236,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.13 0 0 0; 2705 7.34 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(305,kpkm); 2605 7.13 0 0 0; 2605 7.24 0 2 thk2n(286,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(305,kpkm); 2605 7.13 0 0 0; 2605 7.44 0 2 thk2n(286,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(355,kpkm); 2655 7.13 0 0 0; 2655 7.44 0 2 thk2n(236,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.13 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(355,kpkm); 2655 7.13 0 0 0; 2655 7.4 0 2 thk2n(236,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.13 0 0 0; 2705 7.4 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 6.9 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 6.95 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.06 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.1 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'qu' [2705 0 4 0 0; 2705 39 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'qk' [2705 0 4 0 0; 2705 7228 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.06 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'qu' [2705 0 4 0 0; 2705 78 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'qk' [2705 0 4 0 0; 2705 14456 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.06 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'qu' [2705 0 4 0 0; 2705 156 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'qk' [2705 0 4 0 0; 2705 28912 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.06 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'qu' [2705 0 4 0 0; 2705 624 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'qk' [2705 0 4 0 0; 2705 99999.99 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.06 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'qu' [2705 0 4 0 0; 2705 1248 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'qk' [2705 0 4 0 0; 2705 99999.99 0 0 thk2n(186,kpkm); 2891 0 2 0 0] 'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.06 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.13 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0] 'vs' [2630 0 4 0 thk2n(100,kpkm); 2730 0 4 0 0]} ...
                {'vs' [2300 0 4 0 thk2n(405,kpkm); 2705 7.13 0 0 0; 2705 7.44 0 2 thk2n(186,kpkm); 2891 7 0 0 0] 'vs' [2705 7.34 0 0 thk2n(85,kpkm); 2790 0 4 0 0]}};
        otherwise
            error('seizmo:cmb_1dmodel_library:badMODNAME',...
                'Unknown MODNAME: %s',varargin{i});
    end
    
    % make sure name is a cell array
    if(ischar(name)); name=cellstr(name); end
    if(nargout<2)
        for j=1:numel(pv)
            cnt=cnt+1;
            models(cnt)=perturb_1dmodel(startmod,name{j},pv{j}{:});
        end
        varargout={models};
    else
        varargout={name pv};
    end
end

end


function [n]=thk2n(thk,kpkm)
n=ceil(1+thk*kpkm);
end


