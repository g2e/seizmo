function model=mkemptymodel
% mkemptymodel.......create empty MODEL structure
%
% call: model=mkemptymodel;
%
% result: model: empty MODEL structure, to be filled with a discrete
%                velocity model.
%                This will be a structure with the following fields:
%
%                 model.z: depth [km below surface]
%                 model.vp: Vp
%                 model.vs: Vs
%                 model.rho: rho [g/ccm]
%                 model.qp: Qp
%                 model.qs: Qs
%                 model.conr: depth of conrad discontinuity
%                 model.moho: depth of moho
%                 model.d410: depth of Mantle-Transition Zone-discontinuity (the "410" on earth, where
%                             the olivine-alpha-beta-transition occurs)
%                 model.d520: depth of the olivine-beta-gamma transition (520km on earth)
%                 model.d660: depth of lower mantle discontinuity (the "660" on earth, where the 
%                             olivine-gamma-perovskite transition occurs)
%                 model.cmb: depth of core mantle boundary
%                 model.icb: depth of inner core boundary
%                 model.dz: depths of additional discontinuities
%                 model.dname: names of additional discvontinuities
%                 model.rp: planetary radius
%                           Defined usgin the !radius-Keyword. If this is nbot given,
%                           the largest value of z is assumed to be the planetary radius!
%                 model.name: name of model, determined by !name keyword
%                 model.year: year of publication of model, defined through the !year-Keyword.
%                             If this is not given, the year will be empty
%
%                 in order to identify the standard discontinuities within the .nd-file,
%                 use the following discontinuity names:
%                 "conrad" -> model.conr
%                 "moho" or "mantle" -> model.moho
%                 "transition zone" or "olivine alpha beta" -> model.d410
%                 "olivine beta gamma" -> model.d520
%                 "lower mantle" or "olivine gamma perovskite" -> model.d660
%                 "outer core" or "outer-core" -> core mantle boundary
%                 "inner core" or "inner-core" -> inner core boundary
%                 names are not case sensitive.
%
% Martin Knapmeyer, 25.09.2006
% Garrett Euler,    18.02.2006

%%% init model
model=struct(....
    'name',[],...
    'rp',[],...
    'year',[],...
    'z',[],...
    'vp',[],...
    'vs',[],...
    'rho',[],...
    'qp',[],...
    'qs',[],...
    'conr',NaN,...
    'moho',NaN,...
    'd410',NaN,...
    'd520',NaN,...
    'd660',NaN,...
    'cmb',NaN,...
    'icb',NaN,...
    'dz',[],...
    'dname',[]);

end
