function dmy=mkaxis(name,values);
% mkaxis................set axis range for a single axis
%
% call: mkaxis(name,values);
%
%              name: axis designator - decide which axis to set
%                    possible are:  'x'
%                                   'y'
%                                   'z'
%                    (not case sensitive)
%
%              values: vector containing axis range (min and max)
%                      the function sorts VALUES internally.
%
% result: dmy: always zero
%
% this is a shortcut to several set(gca,'XLim',value) calls.
%
% martin Knapmeyer 20.04.1998

name=lower(name);
values=sort(values);

ax=axis;

switch name
    case {'x'},
         set(gca,'xlim',values);
    case {'y'},
         set(gca,'ylim',values);
    case {'z'},
         set(gca,'zlim',values);
end; % switch
