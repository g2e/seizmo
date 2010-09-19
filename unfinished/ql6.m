function [model]=ql6()

% Durek & Ekstrom 1996, A radial model of anelasticity consistent with
% long-period surface-wave attenuation
model.depth=[0 24.4 24.4 80 80 220 220 670 670 2891 2891 5150 5150 6371]';
model.qu=[300 300 191 191 70 70 165 165 355 355 0 0 104 104]';
model.qk=[1e9 1e9 943 943 943 943 943 943 1e9 1e9 1e9 1e9 1e9 1e9]';

end