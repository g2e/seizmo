function ynew=mk2ptinterp(x,y,xnew);
% mk2ptinterp........linear interpolation between 2 points
%
% call: ynew=mk2ptinterp(x,y,xnew);
%
%       x: two elements vector defining an x-axis interval
%
%       y: two elements vector defining the y-axis interval
%
%       xnew: scalar x value, for which ynew is to be computed
%             Since a straight line is assumed, XNEW does not need to be
%             between the two point defined in X. Extrapolation is easy here.
%
% result: ynew: interpolation ynew(xnew) of the straight line y(x) as
%               defined by input values X and Y.
%               If both X values are inf, NaN is returned.
%               if the gradient defined by X and Y is inf, NaN is returned.
%
% This routine computes a straight forward, non-vectorized linear
% interpolation between exactly two points.
%
% One exception has been implemented: if one of the X values is inf, the
% interpolation always returns the y associated with the other x. This was
% introduced to solve a problem with ScS vertex depth, but it makes sense
% in some asymptotic interpretaiton.
% 
%
% Martin Knapmeyer, 24.05.2006, 27.10.2006, 05.12.2006

% 27.10.2006: if x(2) is inf, return y1.
% 05.12.2006: if gradient is inf, return NaN

%%% make x(2)>x(1)
[x,sorter]=sort(x);
y=y(sorter);


infties=find(isinf(x));
switch length(infties)
    case {0}
        %%% both x values are finite, so do a simple interpolation
        if x(2)-x(1)~=0
            %%% evaluate line equation
            ynew=((x(2)-xnew)*y(1)+(xnew-x(1))*y(2))/(x(2)-x(1));
        else
            %%% gradient is inf, return NaN
            ynew=NaN;
        end; % if x(2)-x(1)~=0
    case {1}
        %%% one x is infinite, return the other y.
        if isinf(x(1))
           %%% x(1) is infinite, return y(2)
           ynew=y(2);
        else
           %%% x(2) is infinite, return y(1)
           ynew=y(1);
        end; % if isinf(x(1))
    case {2}
        %%% both x values are infinite, so wew can't do anything and return
        %%% a NaN.
        ynew=NaN;
end; % switch length(infties)


% %%% evaluate line equation
% ynew=((x(2)-xnew)*y(1)+(xnew-x(1))*y(2))/(x(2)-x(1));


% %%% control plot
% plot(x,y,'.-');
% hold on
% plot(xnew,ynew,'o');
% hold off


