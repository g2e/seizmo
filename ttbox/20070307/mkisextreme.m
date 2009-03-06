function yesno=mkisextreme(y);
% mkisextreme........check if point sequence constitutes an extremum
%
% call: yesno=mkisextreme(y);
%
%       y: three element vector of y=f(x) values
%
%          y is considered to be a function of some implied x.
%
% result: yesno: flag denoting if y(x(2)) constitutes a local maximum,
%                minimum or the edge of a plateau
%
%                possible values: 1: is an extremum or plateau edge
%                                 0: is not
%
%
% The result flag YESNO is zero if and only if y(2) is truly between y(1)
% and y(3), but not if it is a local minimum, a local maximum, or if it
% forms the edge of a plateau with one of the other two values.
% To decide if a maximum/minimum/plateau is present, the x values are not
% necessary, the y values are sufficient.
%
% This routine is not a generic extrema-detector, but written specifically
% for the purposes of MKIMPROVESAMPLING.
%
% Martin Knapmeyer, 12.06.2006, 21.06.2006


%%% here comes a magic number: not any numerical noise in the function
%%% schould be accepted as an extremum. PEPSILON defines a threshold:
%%% extreme must have function value ratios higher than this. The big
%%% question is: how to choose this number... 1.04 gave reasonable
%%% results in tests with IASP91, PREM and AK135.
pepsilon=1.04; %1.3;


%%% comparison with previous sample only: possibly just a jump
if ((y(2)/y(1)>pepsilon))... % y(2) is a local max or max-like plateau
   ||...
   ((y(1)/y(2)>pepsilon)) % y(2) is a local min or min-like plateau
  
   yesno=1;
   
   %disp(['MKISEXTREME: function ratios: ' num2str(y(2)/y(1)) ' ' num2str(y(1)/y(2))]);
  
else
    
   yesno=0;
    
end; % if ...




