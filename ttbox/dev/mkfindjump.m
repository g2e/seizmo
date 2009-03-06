function jumps=mkfindjump(x,y);
% mkfindjump.......find point at which y(x) has jumps
%
% call: jumps=mkfindjump(x,y);
%
%       x: aequidistant x coordinates
%       y: function value y(x) at given x values
%
% result: jumps: element index at which function y(x) jumps.
%                empty if there are no jumps.
%
% A jump is a type of discontinuity, at which limit of
% y=f(x) for x-epsilon is different from the limit for x+epsilon
% as epsilon approaches zero.
% It should be sufficient to compare an extrapolated function value
% with the actual value to identify jumps.
% This routine does not search for poles or holes in the definition.
% Poles with sign change will be detected as jumps. Poles without sign
% change will be found in some cases, in other cases the will be ignored
% - this depends on the sampling of the function.
% Holes (where y==NaN) will be ignored.
%
% It is assumed that x is sampled regularly, even at holes. This means: if your
% function is not defined in a certain interval, x has to be sampled regularly
% within this interval, and y should be set to NaN there.
%
% the idea behind the algorithm is:  use the n-th and (n+1)-th values
% of x and y to predicht the (n+2)-th value of y.
% Under the assumption that x is sampled regularly, the linear extrapolation
% from two given points becomes simple.
% At a jump, the predicted value will depart heavily from the true one.
%
% Due to the linear extrapolation, the routine will not find a jump within the
% first two samples of the function, and also not within the first two samples
% after a hole.
%
% Martin Knapmeyer, 27.10.2003

%tic; 

%%% init result
jumps=[];

%%% enough data?
xlen=length(x);
if length(y)~=xlen
   error('MKFINDJUMP: lengths of X and Y have to be equal!');
end; % if length(y)

%%% plot data
%plot(x,y,'k.');

%%% the n-th and (n+1)-th values of x and y
xn=x(1:(xlen-2));
xnplus1=x(2:(xlen-1));
yn=y(1:(xlen-2));
ynplus1=y(2:(xlen-1));

%%% predict (n+2)-th value
%%% at this point we assume regular sampling to
%%% simplify the linear extrapolation.
xnplus2=x(3:xlen);
deltay=(ynplus1-yn);
ynplus2=ynplus1+deltay;

% %%% control plot
% hold on
% plot(xnplus2,ynplus2,'go');
% hold off

%%% determine which deviation is associated with a jump
maxdiff=abs(deltay);
deviation=abs(ynplus2-y(3:xlen));
indies=find(deviation>maxdiff);
if ~isempty(indies)

   %%% cut out double jumps: a jump downward will  result in two indexes here
   jumps=indies(1);
   for indy=2:length(indies)
       if indies(indy)~=indies(indy-1)+1
          jumps=[jumps indies(indy)];
       end;
   end; % for indy
   jumps=jumps-1;
   
   %%% ynplus2 omits the first two elements of y. To get a valid index into y
   %%% we therefore have to add 2 to jumps.
   jumps=jumps+2;

   %%% control plot
   %hold on
   %ends=jumps;
   %starts=jumps+1;
   %plot(x(starts),y(starts),'r>');
   %plot(x(ends),y(ends),'r<');
   %hold off
   
   
   
else
   %%% else there is no jump, and the result remains empty.
end; % if ~isempty(indies)

%toctime=toc; disp(['MKFINDJUMP: elapsed time: ' num2str(toctime) 'sec.']);

return;




