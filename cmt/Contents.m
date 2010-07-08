% Seismology Toolbox - cmt
% Version 0.6.0-r150 Everest 8-July-2010
%
% Various Moment Tensor functions
%AR2HRV        - Convert moment tensor from Aki & Richards form to Harvard form
%AUXPLANE      - Returns strike-dip-slip of 2nd (auxiliary) focal plane
%ELEMENTARY_MT - Returns one of six elementary moment tensors
%HRV2AR        - Convert moment tensor from Harvard from to Aki & Richards form
%MT2SDR        x Convert moment tensor to strike-dip-rake
%MT2TPB        x 
%MT_62V        - Converts moment tensor from 6 Nx1 vectors to Nx6 array
%MT_DECOMP     x 
%MT_DIAG       x Returns diagonalized moment tensors
%MT_G2V        - Convert moment tensor from 3x3xN to Nx6
%MT_NORM       x Returns moment tensor normalized by its scalar moment
%MT_V26        - Converts moment tensor from Nx6 array to 6 Nx1 vectors
%MT_V2G        - Convert moment tensor from Nx6 to 3x3xN
%NS2SDR        x Returns strike-dip-rake for a given normal and slip vector
%NS2TPB        x 
%SDR2MT        x Convert strike-dip-rake to moment tensor in Aki & Richards form
%SDR2NS        x Returns normal & slip vectors (NEU) for a given strike-dip-rake
%STRIKEDIP     - Returns strike & dip given normal vector to plane
%TPB2MT        x 
%TPB2NS        x 
%
% functions to add after the above
%  checkmt
%  checkneu
%  checksdr
%  mt_info
%  rotatemt
%  plotmt
%  nodallines
%  plotmt3
%  nodallines3
%
%  getcmts (for a region)
%  getcmt
%  cmt2mt
