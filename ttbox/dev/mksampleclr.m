function [depthout,sample]=mksampleclr(clr,depth,modifier);
% mksampleclr........evaluate CLR structure at given depth
%
% call: [depthout,sample]=mksampleclr(clr,depth,modifier);
%       [depthout,sample]=mksampleclr(clr,depth,modifier,mode);
%
%       clr: CLR structure as returned by MKREADCLR (see there for definition)
%       depth: depth at which layer polynomials are to be evaluated [km]
%              this might be a vector!
%       modifier: string that defines which quantity to evaluate
%                 possible values:
%                 'vp': P wave velocity [km/s]
%                 'vs': S wave velocity [km/s]
%                 'rho': density        [g/ccm]
%                 'qp': P wave Q-factor 
%                 'qs': S wave Q factor
%       mode: string defining the model evaluation mode. available are:
%             'strict': the model is evaluated as-is
%             'protected': negative values are replaced by zeroes
%                          Since negative values as meangingless for all
%                          quantities stored in CLR files, and since the use
%                          of negative velocities would crash TTBOX, this mode
%                          protects you and TTBOX from numerical artefacts of
%                          crude modelling.
%             DEFAULT: 'protected'
%             Mode strings are not case sensitive.
%
% result: depthout: depth for which CLR is evaluated.
%                   sample(i) is the model value at depth(i).
%                   DEPTHOUT may be longer than DEPTH since two values
%                   are returned at every discontinuity!
%         sample: value(s) of model at depth DEPTH, corresponding to MODIFIER.
%                 This might be a vector, since models are not unique at first
%                 order discontinuities. All Layer top and bottom depths will
%                 be sampled additionally. Discontinuities will be preserved
%                 and produce two samples each.
%                 NaN is returned if the chosen parameter is not defined.
%  
% Martin Knapmeyer, 03.09.2003, 19.09.2003


%%% init result
depthout=[];
sample=[];

%%% understand input
nin=nargin;
switch nin
   case {3}
      mode='protected';
   case {4}
      mode=lower(mode);
   otherwise
      error('MKSAMPLECLR: illegal number of input arguments');
end; % swtich nin

%%% loop through all layers and evaluate those
%%% which are defined at given DEPTH
for indy=1:clr.lyrcnt


    indies=find((min(clr.layers(indy).depth)<=depth)&...
               (max(clr.layers(indy).depth)>=depth));
               
    if ~isempty(indies)
    
    

       %%% choose polynomial
       switch lower(modifier)
           case {'vp'}
                polynom=clr.layers(indy).vp;
           case {'vs'}
                polynom=clr.layers(indy).vs;
           case {'rho'}
                polynom=clr.layers(indy).rho;
           case {'qp'}
                polynom=clr.layers(indy).qp;
           case {'qs'}
                polynom=clr.layers(indy).qs;
           otherwise
                error(['MKPLOTCLR: unknwon quantity ' upper(modifier)]);
       end; % switch modifier
       
       %%% evaluate polynomial
       x=(clr.rp-depth(indies))./clr.rp; % reduced radius
       quantity=polyval(polynom(end:-1:1),x);
       depthout=[depthout depth(indies)];
       sample=[sample quantity];
    
    end; % if (min...
    
end; % for indy


%%% apply MODE
switch mode
   case {'strict'}
      %%% nothing happens
   case {'protected'}
      %%% replace negative values by zero
      indies=find(sample<0);
      if ~isempty(indies)
         sample(indies)=sample(indies)*0;
      end; % if ~isempty
   otherwise
      error(['MKSAMPLECLR: unknown mode ' upper(mode)]);
end; % switch mode


%%% return result
if isempty(sample)
   sample=NaN;
else
   sample=sample(:);
   depthout=depthout(:);
end; %  if isempty