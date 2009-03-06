function handle=mkplotclr(clr,plotwhat,style);
% mkplotclr........plot continuous model from CLR structure
%
% call: handle=mkplotclr(clr,plotwhat);
%       handle=mkplotclr(clr,plotwhat,style);
%
%       clr: CLR structure as returned by MKREADCLR.
%       plotwhat: quantity to be plotted. This is a string which
%                 specifies which of the layer parameters given
%                 in CLR is to plot.
%                 possible values:
%                 'vp': P wave velocity
%                 'vs': S wave velocity
%                 'rho': density
%                 'qp': P wave Q-factor
%                 'qs': S wave Q factor
%       style: plot style. available options are:
%              'raw': plot only the polynomials.
%                     This is the default
%              'fancy': add discontinuities and dotted extensions
%                       of polynomilas downto clr.rp depth
%                       This is mainly for demonstration purpopses.
%              'nd': plot model in the style of MKPLOTMODEL.
%
% result: handle: handle to the line object of the model curve
%
% Plots layer polynomials for a chosen quantity, discontinuities and layer names
% into a figure. Layer polynomials are plotted for a) the entire depth range and
% b) the depth range for which they are valid.
%
% Martin Knapmeyer, 02.09.2003


%%% understand input
switch nargin
   case {2}
      style='raw';
   case {3}
      if strcmp(style,'nd')==1
         %%% produce a plot using MKPLOTMODEL
         model=mkclr2model(clr,50,'spherical');
         mkplotmodel(model);
         return;
      end; % if strcmp(style,'nd')==1
   otherwise
       error('MKPLOTCLR: illegal number of input args!');
end; % switch nargin

%%% prepare plot
if ~ishold
   cla;
end; % if ~ishold
set(gca,'ytickmode','auto');

%%% constuct layer-sequence ordered by depth
layerdepths=zeros(1,clr.lyrcnt);
for indy=1:clr.lyrcnt
    layerdepths(indy)=min(clr.layers(indy).depth);
end; % for indy
[layerdepths,layersequence]=sort(layerdepths);

%%% plot layers
alldepths=[];
allquantity=[];
for indy=1:clr.lyrcnt

    currentlyr=layersequence(indy);

    %%% choose polynomial
    switch lower(plotwhat)
        case {'vp'}
             polynom=clr.layers(currentlyr).vp(end:-1:1);
        case {'vs'}
             polynom=clr.layers(currentlyr).vs(end:-1:1);
        case {'rho'}
             polynom=clr.layers(currentlyr).rho(end:-1:1);
        case {'qp'}
             polynom=clr.layers(currentlyr).qp(end:-1:1);
        case {'qs'}
             polynom=clr.layers(currentlyr).qs(end:-1:1);
        otherwise
             error(['MKPLOTCLR: unknwon quantity ' upper(plotwhat)]);
    end; % switch plotwhat
    
    %%% evaluate & plot polynomial at all depths below layer top
    %%% note that coefficients in POLYVAL are sorted backwards...
    if strcmp(style,'fancy')
       depths=0:1:clr.rp;
       x=(clr.rp-depths)/clr.rp; % reduced radius
       quantity=polyval(polynom,x);
       hold on
       plot(quantity,depths,'k:');
       hold off
    end; % if strcmp
    
    %%% evaluate only at depths at which layer is
    depths=min(clr.layers(currentlyr).depth):1:max(clr.layers(currentlyr).depth);
    depths=[min(clr.layers(currentlyr).depth) depths max(clr.layers(currentlyr).depth)];
    x=(clr.rp-depths)/clr.rp; % reduced radius
    quantity=polyval(polynom,x);
    alldepths=[alldepths depths];
    allquantity=[allquantity quantity];
    
end; % for indy

%%% plot for depths at which definitions are valid
hold on
handle=plot(allquantity,alldepths,'b-');
hold off

%%% add standard discontinuities
if strcmp(style,'fancy')
   ax=axis;
   xrange=[min(0,ax(1)) ax(2)];
   hold on
   plot(xrange,[1 1]*clr.conr,'b--');
   plot(xrange,[1 1]*clr.moho,'b--');
   plot(xrange,[1 1]*clr.d410,'b--');
   plot(xrange,[1 1]*clr.d520,'b--');
   plot(xrange,[1 1]*clr.d660,'b--');
   plot(xrange,[1 1]*clr.cmb,'b--');
   plot(xrange,[1 1]*clr.icb,'b--');
   hold off


   %%% add non-standard discontinuities
   dzcnt=length(clr.dz);
   for indy=1:dzcnt
       hold on
       plot(xrange,[1 1]*clr.dz(indy),'b--');
       hold off
   end; % for indy



   %%% add layer boundaries and names
   for indy=1:clr.lyrcnt
       text(max(xrange)/2,...
            mean(clr.layers(indy).depth),...
            clr.layers(indy).name);
         
       %%% plot layer boundaries
       hold on
       plot(xrange,...
            [1 1]*min(clr.layers(indy).depth),'k:');
       plot(xrange,...
            [1 1]*max(clr.layers(indy).depth),'k:');
       hold off
   end; % for indy

end; % if strcmp

%%% beautify plot
box on
%grid on
axis ij
axis auto
ax=axis;
axis([min(0,ax(1)) ax(2) 0 clr.rp]);
ylabel('depth [km]');
ytick=get(gca,'YTick');
ytick=[ytick clr.rp];
set(gca,'YTick',ytick);
ylbl=strvcat(get(gca,'YTickLabel'),num2str(clr.rp));
set(gca,'YTickLabel',ylbl);
switch lower(plotwhat)
    case {'vp'}
         xlbl='P wave velocity [km/s]';
    case {'vs'}
         xlbl='S wave velocity [km/s]';
    case {'rho'}
         xlbl='Dinsity [g/ccm]';
    case {'qp'}
         xlbl='P wave Q-Factor';
    case {'qs'}
         xlbl='S wave Q-factor';
    otherwise
        error(['MKPLOTCLR: unknwon quantity ' upper(plotwhat)]);
end; % switch plotwhat
xlabel(xlbl);
title(['Model ' clr.name ' (' int2str(clr.year) '): Layer Polynomials and Model']);