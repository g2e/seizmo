function [h,m,s]=mks2hms(s,verbose);
% mks2hms........convert seconds with decimal fractions into hour, min, sec
%
% call: [h,m,s]=mks2hms(s);
%       [h,m,s]=mks2hms(s,'verbose');
%       hms=mks2hms(s);
%
%       s: number of seconds (with decimal fractions)
%       verbose: if present, the result is prettyprinted to the screen.
%
% result: h: hours
%         m: minutes
%         s: seconds
%         hms: a sting containing the time (only if the input is not a vector)
%
% Martin Knapmeyer, 26.04.2002


h=fix(s/3600);
s=s-h*3600;

m=fix(s/60);
s=s-m*60;

s=s;

if nargin==2
   for i=1:length(s)
      hstr=mktimestring(h(i));
      mstr=mktimestring(m(i));
      sstr=mktimestring(s(i));
      %disp([int2str(h(i)) ':' int2str(m(i)) ':' num2str(s(i))]);
      disp([hstr ':' mstr ':' sstr]);
   end; % for i
end; % if nargin

if nargout<=1
   %%% return string with prettyprinted time
   hstr=mktimestring(h);
   mstr=mktimestring(m);
   sstr=mktimestring(s);
   %disp([int2str(h(i)) ':' int2str(m(i)) ':' num2str(s(i))]);
   h=[hstr ':' mstr ':' sstr];
end; % if nargout==1

return

%%%%%%%%%%%% Helper functions

function s=mktimestring(x)

if x<10
   s=['0' num2str(x)];
else
   s=num2str(x);
end;

return;