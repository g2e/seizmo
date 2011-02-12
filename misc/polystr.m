function [str]=polystr(p,var)
%POLYSTR    Converts polynomial coefficients into string form
%
%    Usage:    str=polystr(p)
%              str=polystr(p,var)
%
%    Description:
%     STR=POLYSTR(P) creates a string from the polynomial coefficients in
%     P.  P must be in descending order (as from POLYFIT).  The string
%     includes 'x' as the default variable.
%
%     STR=POLYSTR(P,VAR) changes the variable in the polynomial equation to
%     the string VAR.
%
%    Notes:
%
%    Examples:
%     % Linear fit to some random points:
%     x=10*rand(1,100);
%     y=sin(x)+0.1*rand(1,100);
%     p=polyfit(x,y,5);
%     plot(x,y,'ko');
%     hold on
%     plot(sort(x),polyval(p,sort(x)));
%     hold off
%     title(['y=' polystr(p)])
%
%    See also: POLYFIT, POLYVAL

%     Version History:
%        Sep. 17, 2010 - initial version
%        Feb. 11, 2011 - minor mlint fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 20:00 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default variable
if(nargin<2 || isempty(var)); var='x'; end

% check inputs
if(~isreal(p) || ~isvector(p))
    error('seizmo:polystr:badInput',...
        'P must be a vector of real-valued polynomial coefficients!');
end
if(~isstring(var))
    error('seizmo:polystr:badInput',...
        'VAR must be a character string!');
end

% convert
np=numel(p);
for i=1:np
    if(i==1)
        str=sprintf('%g',p(i));
    else
        str=[str sprintf('%+g',p(i))];
    end
    if(np-i>1)
        str=[str var '^' num2str(np-i)];
    elseif(np-i>0)
        str=[str var];
    end
end

end
