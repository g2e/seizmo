function dmy=mkwritend(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12);
% MKWRITEND.........write velocity model in named discontinuity format to disk
%
% call: mkwritend(model,pfad);
%       mkwritend(name,year,radius,z,vp,vs,rho,qp,qs,dz,dname,pfad);
%
%       model: model structure as returned by MKREADND (see there for definition)
%       pfad: path to write .nd file to
%
%       name: model name [string]
%       year: model publication year
%       radius: planetary radius [km]
%       z: layer depths (upper end) in km
%       vp: P wave velocity [km/s], has to be of same length as z
%       vs: S wave velocity [km/s], has to be of same length as z
%       rho: density [g/ccm], has to be of same length as z
%       qp: P wave Q factor, has to be of same length as z
%       qs: S wave Q factor, has to be of same length as z
%       dz: discontinuity depths
%       dname: discontinuity names
%
% results: dmy: -1 in case of error
%                0 otherwise
%
% At discontinuity depths, the model has to define to values for each parameter:
% one value for below and one value for above.
% values of dz have to correspond exactly to values in z.
%
% If any value is not defined, set NaN into the corresponding array element. This
% will be replaced by -1 in the file (for compatibility reasons with MKREADND)
% -1 is physically nonsense for all parameters and may therefore serve as
% identifier for invalid or undefined parameters.
%
% Martin Knapmeyer, 10.04.2002, 04.09.2003

%%% init result
dmy=0;


%%% interpret input
nin=nargin;
switch nin
    case {2}
       %%% transform input into 11 separate arrays
       name=p1.name;
       year=p1.year;
       radius=p1.rp;
       z=p1.z;
       vp=p1.vp;
       vs=p1.vs;
       rho=p1.rho;
       qp=p1.qp;
       qs=p1.qs;
       dz=[p1.dz,...
           p1.conr,...
           p1.moho,...
           p1.d410,...
           p1.d520,...
           p1.d660,...
           p1.cmb,...
           p1.icb];
       dname=strvcat(p1.dname,...
                     'conrad',...
                     'moho',...
                     'olivine alpha beta',...
                     'olivine beta gamme',...
                     'olvine gamma perovskite',...
                     'cmb',...
                     'icb');
       indies=find(~isnan(dz));
       dz=dz(indies);
       dname=dname(indies,:);
       [dz,sorter]=sort(dz);
       dname=dname(sorter,:);
       %%% savepath
       pfad=p2;
    case {12}
       name=p1;
       year=p2;
       radius=p3;
       z=p4;
       vp=p5;
       vs=p6;
       rho=p7;
       qp=p8;
       qs=p9;
       dz=p10;
       dname=p11;
       pfad=p12;
    otherwise
       error('MKWRITEND: illegal number of input parameters');
end; % switch nin


%%%%% set all NaNs to -1
% z
indy=find(isnan(z));
z(indy)=zeros(size(z(indy)))-1;
% vp
indy=find(isnan(vp));
vp(indy)=zeros(size(vp(indy)))-1;
% vs
indy=find(isnan(vs));
vs(indy)=zeros(size(vs(indy)))-1;
% rho
indy=find(isnan(rho));
rho(indy)=zeros(size(rho(indy)))-1;
% Qp
indy=find(isnan(qp));
qp(indy)=zeros(size(qp(indy)))-1;
% Qs
indy=find(isnan(qs));
qs(indy)=zeros(size(qs(indy)))-1;


%%%%% open file
[fid,msg]=fopen(pfad,'w');
if fid==-1
   warning(['MKWRITEND: ' msg]);
   dmy=-1;
else
   %%%%% write Miscellaneous Information
   fprintf(fid,'!name %s\n',name);
   fprintf(fid,'!year %s\n',year);
   fprintf(fid,'!radius %s\n',radius);

   %%%%% write parameter table
   writtend=1; % number of already written discontinuities
   anz=length(z); % number of depth layers
   for indy=1:anz
      %% write parameter set
      fprintf(fid,'%7.3f %7.5f %7.5f %7.5f %7.5f %7.5f\n',z(indy),vp(indy),vs(indy),rho(indy),qp(indy),qs(indy));
      if writtend<=length(dz)
         if z(indy)==dz(writtend) % writtend may have illegal values here, so do not combine the IFs!
            %% write discontinuity name
            %disp(['MKWRITEND: writing ' dname(writtend,:)]);
            fprintf(fid,'%s\n',dname(writtend,:));
            writtend=writtend+1;
         end; % if z(indy)
      end; % if writtend
   end; % for indy
   %%%%% close file
   fclose(fid);
end; % if fid==-1 else


%%%%% return result
dmy=dmy; % has been defined above.