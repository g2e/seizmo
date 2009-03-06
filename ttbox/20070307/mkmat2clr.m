function clrout=mkmat2clr(tag,name,year,planet,rp,z,vp,vs,rho,qp,qs,dz,dname);
% mkmat2clr.........convert layer polynomial matrices into clr structure
%
% call: clrout=mkmat2clr(tag,name,year,planet,rp,z,vp,vs,rho,qp,qs,dz,dname);
%     
%       tag: user tag string
%            this can be used to add some extra information, which will be read
%            from the file when loading it (comments are ignored)
%       name: model name (string)
%       planet:name of planet for which model is valid [string]
%       rp: radius of planet [km]
%       z: layer depth matrix (depth! D.E.P.T.H.!)
%          z is a n-by-2 matrix containing minimum and maximum depth of n layers
%          z(i,1) is the depth to the top of the i-th layer [km]
%          z(i,2) is the bottom of the i-th layer [km]
%       vp: P velocity layer polynomial matrix
%           vp is a n-by-m matrix containing the coefficients of layer polynomials
%           on f layers. The highest degree polynomial is of degree m-1.
%           vp(i,:) are the coefficients of the i-th layer in km/s.
%           vp(i,1) is the constant term, vp(i,2) the coefficient of the linear term,
%           vp(i,3) the coefficient of the qudratic term etc.
%       vs: as vp, but for S velocity
%       rho: as vp, but for density
%       qp: as vp, but for dimensionless Qp factor
%       qs: as vp, but for dimensionless Qs factor
%       dz: n elements vector containing depths of ALL discontinuities (standard
%           and non-standard discontinuities!) in km
%           dz(i) is the depth to the i-th discontinuity [km]
%       dname: string matrix containing the names of all discontinuities
%              dname(i,:) is the name of the discontinuity at depth dz(i).
%
%       If one of the parameters (probably density and Q factors) are not defined
%       in your model, insert NaNs into the matrices.
%       Since it must be clear which layer is undefined, NaN matrices must have
%       the same number of rows as z.
%
% result: clrout: a CLR structure as returned by MKREADCLR. (see there for definition)
%
% Depths of standard discontinuities are stored in the respective fields of the CLR
% structure. non-standard discontinuities remain in the dz and dname fields.
%
% Martin Knapmeyer, 03.09.2003


%%% init result
clrout.tag=[];
clrout.name=[];
clrout.year=[];
clrout.planet=[];
clrout.rp=[];
clrout.lyrcnt=0;
clrout.layers=[];
clrout.conr=NaN; %[];
clrout.moho=NaN; %[];
clrout.d410=NaN; %[];
clrout.d520=NaN; %[];
clrout.d660=NaN; %[];
clrout.cmb=NaN; %[];
clrout.icb=NaN; %[];
clrout.dz=[];
clrout.dname=[];


%%% store General Information
clrout.tag=tag;
clrout.name=name;
clrout.year=year;
clrout.planet=planet;
clrout.rp=rp;

%%% store Layer Information
clrout.lyrcnt=size(z,1);
for indy=1:clrout.lyrcnt
    newlayer.depth=z(indy,:)';
    newlayer.name='';
    newlayer.vp=vp(indy,:)';
    newlayer.vs=vs(indy,:)';
    newlayer.rho=rho(indy,:)';
    newlayer.qp=qp(indy,:)';
    newlayer.qs=qs(indy,:)';
    clrout.layers=[clrout.layers newlayer];
end; % for indy

%%% store Discontinuity Information
dzcnt=length(dz);
for indy=1:dzcnt
    currentname=deblank(dname(indy,:));
    switch lower(currentname)
        case {'conrad'}
             clrout.conr=dz(indy);
        case {'moho','mantle'}
             clrout.moho=dz(indy);
        case {'olivine alpha beta','transition zone'}
             clrout.d410=dz(indy);
        case {'olivine beta gamma'}
             clrout.d520=dz(indy);
        case {'olivine gamma perovskite','lower mantle'}
             clrout.d660=dz(indy);
        case {'outer core','outer-core'}
             clrout.cmb=dz(indy);
        case {'inner core','inner-core'}
             clrout.icb=dz(indy);      
        otherwise
           %%% non-standard-discontinuities
           clrout.dz=[clrout.dz dz(indy)];
           clrout.dname=strvcat(clrout.dname,deblank(dname(indy,:)));
    end; % switch currentname
end; % for indy