function model=mkclr2model(clr,deltaz,mode);
% mkclr2model........convert CLR structure into discrete model
%
% call: model=mkclr2model(clr,deltaz);
%       model=mkclr2model(clr,deltaz,mode);
%
%       clr: continuously defined velocity model in CLR structure
%            as returned by MKREADCLR (see there for definition)
%       deltaz: depths or depth increment used in model sampling [km]
%               This might be a scalar or a vector:
%               scalar: a depth increment from which a depth list
%                       0:DELTAZ:clr.rp will be created
%               vector: a depth list which will be used directly.
%               NOTE that, irrespective of the depths, all
%               layer top and bottom depths as defined in CLR will
%               occur within the output. Discontinuities are preserved.
%       mode: sampling mode. available are:
%             'spherical': sampling will be equidistant in spherical coordinates
%             'flat': sampling will be equidistant in flat earth coorindates.
%                     Depths below a radius of DELTAZ will not be sampled (note
%                     that the planets center in flat coordinates is at infinite 
%                     depth!).
%                     This usually produces a huge amount of samples and thus a
%                     long computation time in all evaluations. On the other hand,
%                     results for core phases are better. Use with caution.
%             NOTE that the MODE is ignored if DELTAZ is a predefined list of depths,
%                  since it is meaningless to make a predefined list equidistant!
%             DEFAULT: 'spherical'
%
% result: model: discrete layer velocity model as returned by MKREADND
%                (see there for definition)
%
% Since the travel time computation routines work with discrete models
% only, the CLR storage format would be completely obsolete without a
% mechanism to convert CLR models into discrete ones.
%
% Martin Knapmeyer, 03.09.2003, 10.11.2003, 18.02.2004

%%% use MKEMPTYMODEL, MK25.09.2006


%%% determine sampling mode
if nargin==2
   mode='spherical';
end; % if nargin==2
mode=lower(mode);

%%% init result
model=mkemptymodel;
% model.name=[];
% model.rp=[];
% model.year=[];
% model.z=[];
% model.vp=[];
% model.vs=[];
% model.rho=[];
% model.qp=[];
% model.qs=[];
% model.conr=NaN;
% model.moho=NaN;
% model.d410=NaN;
% model.d520=NaN;
% model.d660=NaN;
% model.cmb=NaN;
% model.icb=NaN;
% model.dz=[];
% model.dname=[];

%%% sort CLR layers by depth
clr=mksortclr(clr);

%%% copy General Information
model.name=clr.name;
model.rp=clr.rp;
model.year=clr.year;

%%% sample Layers
switch mode
   case {'spherical'}
       if length(deltaz)==1
          %%% deltaz is a depth increment, create a list from it
          regsamples=0:deltaz:clr.rp;
       else
          %%% deltazis already a list - use it unchanged!
          regsamples=deltaz;
       end; % if length(deltaz)==1
   case {'flat'}
       if length(deltaz)==1
          %%% deltaz is a depth increment, create a list from it
          maxdep=model.rp-deltaz;
          [dmy,zmax]=mksfer2flat(0,maxdep,model.rp);
          regsamples=0:deltaz:zmax;
          [dmy,regsamples]=mkflat2sfer(regsamples,regsamples,model.rp);
       else
          %%% deltazis already a list - use it unchanged!
          regsamples=deltaz;
       end; % if length(deltaz)==1
   otherwise
       error(['MKCLR2MODEL: unknown sampling mode ' upper(mode)]);
end; % switch mode
depths=sort([regsamples,...
             clr.conr,...
             clr.moho,...
             clr.d410,...
             clr.d520,...
             clr.d660,...
             clr.cmb,...
             clr.icb,...
             clr.dz]);
depths=unique(depths); % doubling of depths at discontinuities is carried out in MKSAMPLECLR!
[depthout,model.vp]=mksampleclr(clr,depths,'vp');
model.z=depthout;
[depthout,model.vs]=mksampleclr(clr,depths,'vs');
[depthout,model.rho]=mksampleclr(clr,depths,'rho');
[depthout,model.qp]=mksampleclr(clr,depths,'qp');
[depthout,model.qs]=mksampleclr(clr,depths,'qs');

%%% sort by depth
[model.z,sorter]=sort(model.z);
model.vp=model.vp(sorter);
model.vs=model.vs(sorter);
model.rho=model.rho(sorter);
model.qp=model.qp(sorter);
model.qs=model.qs(sorter);

%%% copy Discontinuity Information
model.conr=clr.conr;
model.moho=clr.moho;
model.d410=clr.d410;
model.d520=clr.d520;
model.d660=clr.d660;
model.cmb=clr.cmb;
model.icb=clr.icb;
model.dz=clr.dz;
model.dname=clr.dname;