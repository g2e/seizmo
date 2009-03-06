function [newmodel,newdepths]=mkimprovemodel(oldmodel);
% mkimprovemodel........model analysis and sampling improvement
%
% call: [newmodel,newdepths]=mkimprovesmodel(oldmodel);
%
%       oldmodel: MODEL structure as returned by MKREADND or MKCLR2MODEL
%
% result: newmodel: extended MODEL structure: as OLDMODEL, but with an
%                   improved depth sampling and additional fields:
%
%                   .criticalrays: sub structure containing information on 
%                                  critical rays. It consists of the
%                                  following fields:
%
%                                  .z: critical depths [km]
%                                  .p: p wave ray parameter list [sec/deg]
%                                  .s: s wave ray parameter list [sec/deg]
%
%                                  .p(i) and .s(i) are the ray parameters
%                                  needed to reach depth .z(i).
%                                  The lists are NOT sorted by depth!
%                                 
%                                  These ray parameters correspond to
%                                  discontinuities and beginnings of
%                                  triplications. They define the kinks in
%                                  the epicentral distance as function of
%                                  rayparameter graph.
%
%                   If OLDMODEL already contains this field, it will be
%                   overwritten.
%
%
%
%         newdepths: explicit list of the depth samples that have been
%                    added to the model.
%                    This list is not sorted!
%
%         CRITICALDEPTHS and NEWDEPTHS together give a list of depth
%         samples that should be sufficient to describe all kinks and
%         corners of the distance(rayparm) function.
%
%
% The routine adds new depths samples as identified by MKDISCONTINUITYSMP
% to the depth sampling given in the model and identifies those parts of
% the model that are critical for the ray parameter sampling in the
% shooting process.
%
% Martin Knapmeyer, 23.06.2006



%%% a magical number
%%% vertical epsilon value for placing new samples around important parts
%%% of the velocity model
zepsilon=0.001; % [km]


%%% init result
newmodel=oldmodel;


%%%%%% identify discontinuities and construct depth values for new depth
%%%%%% samples based on these
[disconnewdepths,disccriticalz]=mkdiscontinuitysmp(oldmodel);



%%%%%% identify top and bottom depths of low velocity zones
% disp('--- MIMPROVEMODEL: unfertige LVZ-analyse auskommentiert! dran weiterarbeiten! ---');
[lvzextra,lvzcriticalz]=mklvzsmp(oldmodel);


%%%%%% identify kinks in v(z) that produce travel time triplications
%%%%%% Only the beginning of this type of triplication can be detected
%%%%%% directly from the velocity model.
[derivnewdepths,derivcriticalz]=mkderivsmp(oldmodel);


%%%%%% collect all new depths from discontinuity and derivateve analysis
%%%%%% note that MLLVZSMP returns velocities which have to be integrated
%%%%%% into the model separately, therefore the LVZ does not appear in
%%%%%% NEWDEPTHS
newdepths=[disconnewdepths; derivnewdepths;]; 
criticaldepths=[disccriticalz; derivcriticalz; lvzcriticalz]; % DO NOT SORT/UNIQUE THE CRITICAL DEPTHS LIST!!!
%disp(['MKIMPROVEMODEL: ' int2str(length(criticaldepths))...
%      ' critical depths and ' int2str(length(newdepths)) ' new depth samples.' ]);
  

%%%%%% interpolate model to extended depth samples list
%%%%%% This has to be done in spherical coordinates!
alldepths=sort([oldmodel.z; newdepths]); % unique would destroy discontinuities!!
newmodel=mkinterpmodel(newmodel,alldepths,'preserve');


%%%%%% integrate depths and velocities from LVZ analysis
%%%%%% That analysis gives not only depths but also velocities, since the
%%%%%% velocity interpolation for the LVZ has to be done in flat earth
%%%%%% coordinates, whereas the others are in spherical earth
% disp('--- MKIMPRVEMODEL: einfuegen neuer samples nach LVZ analyse auskommentiert! dran weiterarbeiten! ---')
lvzanz=length(lvzextra.z); % number of extra samples from LVZ analysis
for indy=1:lvzanz
    samplesbefore=find(newmodel.z<lvzextra.z(indy));
    samplesafter=(max(samplesbefore)+1):length(newmodel.z);
    if lvzextra.z(indy)~=newmodel.z(samplesafter(1))
        %%% insert new sample only if currently no sample at this depth exists
        newmodel.z=[newmodel.z(samplesbefore); lvzextra.z(indy); newmodel.z(samplesafter)];
        newmodel.vp=[newmodel.vp(samplesbefore); lvzextra.vp(indy); newmodel.vp(samplesafter)];
        newmodel.vs=[newmodel.vs(samplesbefore); lvzextra.vs(indy); newmodel.vs(samplesafter)];
        newmodel.rho=[newmodel.rho(samplesbefore); lvzextra.rho(indy); newmodel.rho(samplesafter)];
        newmodel.qp=[newmodel.qp(samplesbefore); lvzextra.qp(indy); newmodel.qp(samplesafter)];
        newmodel.qs=[newmodel.qs(samplesbefore); lvzextra.qs(indy); newmodel.qs(samplesafter)];
    else
        %disp(['MKIMPROVEMODEL: double sample at ' num2str(lvzextra.z(indy)) ' avoided - no sample insertion.']);
    end; % if lvzextra.z(indy)~=newmodel.z(samplesafter(1))
end; % for indy



%%%%%% add new field to model structure, containing critical ray
%%%%%% parameters!
%%% compute ray parameters from critical depths
[prayp,srayp,raypz]=mkdepth2rayp(newmodel,criticaldepths);
%%% build criticalrays sub srtructure
criticalrays.z=raypz;
criticalrays.p=prayp;
criticalrays.s=srayp;
%%% insert into model structure, overwriting existing one
newmodel.criticalrays=criticalrays;
