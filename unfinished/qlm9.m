function [model]=qlm9()

% Lawrence & Wysession 2006, QLM9: A new radial quality factor (Qu) model
% for the lower mantle 
model.depth=...
    [0; 80; 80; 220; 220; 400; 400; 670; 670; 1000; 1000; 1350; 1350; ...
    1700; 1700; 2050; 2050; 2400; 2400; 2700; 2700; 2800; 2800; 2891];
model.qu=[600; 600; 80; 80; 143; 143; 276; 276; 362; 362; 325; 325; ...
    287; 287; 307; 307; 383; 383; 459; 459; 452; 452; 278; 278];

end
