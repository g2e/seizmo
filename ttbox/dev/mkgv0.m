function [g,v0]=mkgv0(v,z);
% mkgv0.........compute velocity gradient and baseline shift of linear velocity laws
%
% call: [g,v0]=mkgv0(v,z);
%
%       v: velocities at depths given in z [km/s]
%       z: depths at which velocities are given [km], below surface
%
% result: g: velocity gradients [(km/s)/km ]
%         v0: baseline coreection [km/s]
%
%         g and v0 are row vectors with the same number of elements as v and z have
%
% between the depths defined in Z, a linear velocity law is assumed: v(z)=v0+g*z
% This function computes parameters v0 and g of this velocity law.
%
% Martin Knapmeyer, 23.04.2002, 27.05.2005, de-vectorized 20.09.2005



g=(v(1)-v(2))/(z(1)-z(2));
v0=v(2)-z(2)*g;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             DONT'T DELETE THE FOLLOWING OLD CODE!                  %%%
%%%                                                                    %%%
%%%     old vectorized version, may be needed in future versions       %%%
%%%                           (MK20092005)                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%% input has to be in row vectors (why? MK24062005)
% v=v(:)';
% z=z(:)';
% 
% 
% %%% compute velocity gradients within all layers
% %%% to suppres "divide by zero"-warnings, gradient at 1st-order 
% %%% discontinuities is set to inf by hand
% %%% g(1) is the gradient from v(1) to v(2) at depths z(1) to z(2)
% anz=length(v); % so many parameters
% upperindy=1:(anz-1);
% lowerindy=2:anz;
% z1=z(upperindy); % upper edges of layers
% z2=z(lowerindy);   % lower edges of layers
% v1=v(upperindy); % velocities at upper edges
% v2=v(lowerindy);   % velocities at lower edges
% denominator=z1-z2;
% numerator=v1-v2;
% g=numerator+inf; % just an array full of inf
% nonzeros=find(denominator~=0);
% g(nonzeros)=numerator(nonzeros)./denominator(nonzeros); % don't divide by zero!
% 
% %%% the old version, which wastes an incredible amout fo time with
% %%% switching warning off and on
% % warning off    % to suppress ugly "divide by zero" warnings resulting from discontinuities
% % g=(v1-v2)./(z1-z2); % velocity gradients
% % warning on     % allow warnings later, because they'r unexpected
% 
% 
% %%% compute baseline correction
% %%% this is quite simple, when you have g
% v0=v2-z2.*g;
% 
% %%% return results
% %%% results are aleready in the output variables.