function [tpcp,layertimes]=mkanalytictimes(phase,clrpath);
% mkanalytictimes........compute analytic traveltimes at 0deg and 180deg
%
% call: [tpcp,layertimes]=mkanalytictimes(phase,clrpath);
%
%            phase: seismic phase for which travel times are computed.
%                   This routine handles only a small catalogue of
%                   phases at predefined epicentral distances:
%                   'PcP' at 0deg 
%                   'PKiKP' at 0deg
%                   'PKIKP' or 'PKPdf' at 180deg
%                   'ScS' at 0deg
%                   'SKiKS' at 0deg
%                   'SKIKS' or 'SKSdf' at 180deg
%
%                   (in case of the S phases, the program uses Vp
%                    everywhere below the CMB)
%
%            clrpath: path to the CLR file describing the velocity model
%                     It is expected that CMB and ICB are explicitly
%                     defined in the model.
%                     
%
% result: tpcp: travel time in seconds
%               computed in model IASP91 for surface focus 
%
%         layertimes: the time the ray spends in each of the model layers
%                     (this is a two way time: ray goes down and up again)
%
%               For comparison: The IASP91 time for PcP at distance 0deg
%               according to Kennett (1991) is 8:31.28 or 511.28s.
%
% For vertical incidence, the travel time integral reduces to the integral
% over 1/v(r).
% This routine loads the model coefficients from a clr file and evaluates
% the integral layer by layer, but without any interpolation or approximation.
% Thus the result is an exact reference for usual travel time calculators.
% Analytical solutions for polynomial velocity laws of degree 0,1,2,3 are
% implemented and tested. Travel time for other degrees are solved by
% trapezoidal integration.
%
% In case of cubic velocity laws, the analytical solution is a little
% complicated. Due to numerical resolution, imaginary parts may arise.
% These are rounded off according to the quality of the roots of a third
% order polybnomial which have to be determined during computation. If
% the resulting travel time has imaginary parts, check the result by
% replacing the analytical solutions (calls to MKCONSTLAW, MKSQUARELAW
% and MKCUBELAW) against trapezoidal integration using MKPOLYLAW.
%
% Martin Knapmeyer, 06.04.2004, 22.04.2004, 27.04.2004



%%% load CLR file
clr=mkreadclr(clrpath,'silent');

%%% sort CLR layers by depth (makes loop over layers easier)
clr=mksortclr(clr);


%%% determine maximum depth of ray
switch phase
   case {'PcP','ScS'}
      epidist=0;
      maxdepth=clr.cmb; % go down to CMB
   case {'PKiKP','SKiKS'}
      epidist=0;
      maxdepth=clr.icb; % go down to ICB
   case {'PKIKP','PKPdf','SKIKS','SKSdf'}
      epidist=180;
      maxdepth=clr.rp;  % go down to center.
   otherwise
      error(['MKANALYTICTIMES: cannot handle ' phase]);
end; % switch phase

%%% some informational output
disp(['MKANALYTICTIMES: computing '...
      phase]);
disp(['                 for 0km source depth at '...
      num2str(epidist)...
      'deg epicentral distance']);

%%% loop over layers
%%% the travel time of PcP is twice the time a P wave needs to travel
%%% through all layers from focus (assumed at surface) to cor mantle
%%% boundary. So we compute the travel time for each of the layers
%%% separately (IASP91 has 9 layers above CMB) and then add all times
%%% and multiply the sum by two.
layertimes=zeros(1,clr.lyrcnt); % to collect times for individual layers
for layercnt=1:clr.lyrcnt % loop over ALL layers of model

    %%% evaluate only layers above max depth reached by phase
    if max(clr.layers(layercnt).depth)<=maxdepth
       
       %%% get coefficients of layer polynomial
       switch phase
           case {'PcP','PKIKP','PKiKP','PKPdf'}
                coeff=clr.layers(layercnt).vp;
           case {'ScS','SKIKS','SKiKS','SKSdf'}
                coeff=clr.layers(layercnt).vs;
                if min(clr.layers(layercnt).depth)>=clr.cmb
                   coeff=clr.layers(layercnt).vp;
                end; % if sum(coeff)
       end; % switch phase
       degree=length(coeff)-1; % polynomial degree
       
       %%% get depth range in terms of radius
       layertop=clr.rp-min(clr.layers(layercnt).depth);
       layerbottom=clr.rp-max(clr.layers(layercnt).depth);
       
       %%% evaluate integral according to polynomial degree
       switch degree
          case {0}
               layertimes(layercnt)=mkconstlaw(coeff,layerbottom,layertop,clr.rp);
          case {1}
               layertimes(layercnt)=mklinearlaw(coeff,layerbottom,layertop,clr.rp);
          case {2}
               layertimes(layercnt)=mksquarelaw(coeff,layerbottom,layertop,clr.rp);
          case {3}
               layertimes(layercnt)=mkcubiclaw(coeff,layerbottom,layertop,clr.rp);
          otherwise
               layertimes(layercnt)=mkpolylaw(coeff,layerbottom,layertop,clr.rp);
       end; % switch degree
       
       %%% print current layer's time to screen
       %disp(['Layer #' int2str(layercnt),...
       %      ' Degree: ' int2str(degree),...
       %      ' Depth: ' num2str(min(clr.layers(layercnt).depth)) ' - ',...
       %                num2str(max(clr.layers(layercnt).depth)),...
       %      'km Time: ' num2str(layertimes(layercnt)) 'sec']);
             
             
    end; % if max(clr.layer.depth)<=clr.cmb
    
end; % for layer

%%% sum times of individual layers to get tpcp
tpcp=sum(layertimes);

%%% write result to screen
disp(' ');
disp(['t(' phase ',' num2str(epidist) 'deg)='...
       num2str(tpcp) 's = ' mks2hms(tpcp) ' [h:m:s]']);

return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                              HELPERS                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function t=mkconstlaw(coeff,g,h,rp);
% mkconstlaw........analytic travel time for constant velocity
%
% call: t=mkconstlaw(coeff,g,h,rp);
%
%         coeff: coefficient list
%                the velocity law is v(r)=coeff(1)
%             g: layer bottom radius [km]
%             h: layer top radius [km]
%            rp: planetary radius [km]
%
% Martin Knapmeyer, 22.04.2004

%%% init result
t=0;

%%% decompose coefficient list into a form more similar to my paperwork
a=coeff(1);



%%% compute traveltime
t=2*(h-g)/a;


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t=mklinearlaw(coeff,g,h,rp);
% mklinearlaw.......analytic travel time for linear velocity law
%
% call: t=mkconstlaw(coeff,g,h,rp);
%
%         coeff: coefficient list
%                the velocity law is v(r)=coeff(1)+coeff(2)*x
%                with x being the normalized radius
%             g: layer bottom radius [km]
%             h: layer top radius [km]
%            rp: planetary radius [km]
%
% Martin Knapmeyer, 22.04.2004


%%% init result
t=0;

%%% decompose coefficient list into a form more similar to my paperwork
a=coeff(1);
b=coeff(2);


%%% compute traveltime
t=(rp/b)*log(((a*rp+b*h)^2)/...
             ((a*rp+b*g)^2));

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t=mksquarelaw(coeff,g,h,rp);
% mksquarelaw.......analytic travel time for quadratic velocity law
%
% call: t=mkconstlaw(coeff,g,h,rp);
%
%         coeff: coefficient list
%                the velocity law is v(r)=coeff(1)+coeff(2)*x+coeff(3)*x^2
%                with x being the normalized radius
%             g: layer bottom radius [km]
%             h: layer top radius [km]
%            rp: planetary radius [km]
%
% Martin Knapmeyer, 27.04.2004

%%% init result
t=0;

%%% decompose coefficient list into a form more similar to my paperwork
a=coeff(1);
b=coeff(2);
c=coeff(3);


%%% compute traveltime
d=(rp^2)*(4*a*c-b^2);
t=(-4*rp^2)/(sqrt(d))*...
  (-atan((2*h*c+rp*b)/sqrt(d))+...
    atan((2*g*c+rp*b)/sqrt(d)));

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t=mkcubiclaw(coeff,g,h,rp);
% mkcubiclaw........analytic travel time for cubic velocity law
%
% call: t=mkconstlaw(coeff,g,h,rp);
%
%         coeff: coefficient list
%                the velocity law is v(r)=coeff(1)+coeff(2)*x+coeff(3)*x^2+coeff(4)*x^3
%                with x being the normalized radius
%             g: layer bottom radius [km]
%             h: layer top radius [km]
%            rp: planetary radius [km]
%
%
% Martin Knapmeyer, 27.04.2004


%%%%%% init result
t=0;

%%%%%% decompose coefficient list into a form more similar to my paperwork
a=coeff(1);
b=coeff(2);
c=coeff(3);
e=coeff(4);


%%%%%% compute traveltime

%%%%%% analytic solution for roots of critical polynomial n*z^3+q*z+p

	%%% polynomial coefficients
	n=(...
	   27*(a^2)*(e^2)...
	   +4*e*(b^3)...
	   -18*e*a*c*b...
	   -(c^2)*(b^2)...
	   +4*a*(c^3)...
	  )*(rp^6);
	  
	q=(...
	   -3*e*b...
	   +(c^2)...
	  )*(rp^2);
	  
	p=-e;
	
	%%% numerical evaluation using matlab
	%%% this is not uses since it turned out that my analytical
	%%% solution is better.
	%rho=roots([n 0 q p]);
	
	%%% analytical soution for roots
	
	f=(...
	    (...
	      -108*p...
	      +12*sqrt(3)*sqrt(...
	                        4*(q^3)/n...
	                        +27*(p^2)...
	                      )...
	    )...
	    *(n^2)...
	  )^(1/3);
	
	z1=f/(6*n)-2*q/f;
	
	z2=-0.5*(z1+sqrt(-3)*(f/(6*n)+2*q/f));
	
	z3=-0.5*(z1-sqrt(-3)*(f/(6*n)+2*q/f));
	
	rho=[z1 z2 z3]';

	%%% check roots
	%%% if the roots are identified perfectly, threezeros will
	%%% contain three zeroes.
	threezeros=polyval([n 0 q p],rho);
	
	%%% since the roots are probably not perfect, the max deviation
	%%% from zero will be used as an epsilon to round off the final
	%%% result. This is a crude method to estimate the error, but
	%%% perhaps better than nothing.
	epsilon=max(abs(threezeros)); % max deviation from zero
	epsilon=10^(ceil(log10(epsilon))); % rounded to a power of ten
	if epsilon>0
      disp(['MKCUBELAW: eps=' num2str(epsilon) ' estimated roundoff error.']);
	end; % if epsilon>0
	
	
%%% evaluate sum over all roots
rootanz=length(rho);
timelist=rho*0;
for indy=1:rootanz
    q=rho(indy);
    f1=(...
        -2*(c^2)...
        +6*e*b...
       )...
       *q*(rp^2)...
       +3*e;
       
    f2=rp*(...
           (...
            -c*b...
            +9*e*a...
           )...
           *q*(rp^2)...
           +c...
          );
          
    numerator=f1*h+f2;
    denominator=f1*g+f2;
    
    %t=t+q*log(numerator/denominator);
    
    timelist(indy)=q*log(numerator/denominator);
    
end; % for indy

%%% total time
t=sum(timelist)*2*(rp^3);

%%% apply roundoff to real and imaginary part
tr=round(real(t)/epsilon)*epsilon;
ti=round(imag(t)/epsilon)*epsilon;

%%% rounded result
t=tr+sqrt(-1)*ti;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t=mkpolylaw(coeff,g,h,rp);
% mkpolylaw........analytic travel time for cubic velocity law
%
% call: t=mkpolylaw(coeff,g,h,rp);
%
%         coeff: coefficient list
%                the velocity law is 
%                 v(r)=coeff(1)+coeff(2)*x+coeff(3)*x^2+coeff(4)*x^3+...
%                with x being the normalized radius
%             g: layer bottom radius [km]
%             h: layer top radius [km]
%            rp: planetary radius [km]
%
% This routine applies trapezoidal numerical integration, which works
% for arbitrary polynomials.
%
% Martin Knapmeyer, 26.04.2004



%%%%%% init result
t=0;

%%% construct "grid" to integrate on
step=0.01;  % [km]  - this is essentially the depth sampling
r=g:step:h;

%%% disclaimer
disp(['MKANALYTICTIMES: trapezoid integration with '...
      num2str(step)...
      'km vertical resolution used for degree '...
      int2str(length(coeff)-1)...
      ' polynomials!']);

%%% transform radii given in r into normalized radii
x=r/rp;

%%% compute velocity at grid points
v=polyval(coeff(end:-1:1),x);

%%% integrate reciprocal velocities by trapezoidal rule
t=2*trapz(r,1./v);


return;







