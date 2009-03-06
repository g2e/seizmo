function res=mkmakeclr(restype);
% mkmakeclr........return empty CLR structure
%
% call: res=mkmakeclr(restype)
%
%           restype: string parameter that defines what to create:
%                    'clr': CLR structure with empty layer array
%                    'layer': empty layer structure, to be filled into
%                             the .layer array of a CLR structure
%
% res: result structure correpsonding to RESTYPE
%
% Martin Knapmeyer, 11.09.2003

%%% init result
res=[];

%%% construct result
switch lower(restype)
   case {'clr'}
      clr.name=[];
      clr.year=[];
      clr.planet=[];
      clr.rp=[];
      clr.lyrcnt=0;
      clr.layers=[];
      clr.conr=NaN;
      clr.moho=NaN;
      clr.d410=NaN;
      clr.d520=NaN;
      clr.d660=NaN;
      clr.cmb=NaN;
      clr.icb=NaN;
      clr.dz=[];
      clr.dname=[];
      clr.tag=[];
      res=clr;
   case {'layer'}
      layer.depth=[];
      layer.name=[];
      layer.vp=NaN;
      layer.vs=NaN;
      layer.rho=NaN;
      layer.qp=NaN;
      layer.qs=NaN;
      res=layer;
   otherwise
      error(['MKMAKECLR: unknown result type ' upper(restype)]);
end; % switch restype