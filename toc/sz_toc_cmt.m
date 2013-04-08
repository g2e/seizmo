function []=sz_toc_cmt()
% Various Moment Tensor Functions
%<a href="matlab:help ar2hrv">ar2hrv</a>                         - Convert moment tensor from Aki & Richards form to Harvard form
%<a href="matlab:help auxplane">auxplane</a>                       - Returns strike-dip-slip of 2nd (auxiliary) focal plane
%<a href="matlab:help elementary_mt">elementary_mt</a>                  - Returns one of six elementary moment tensors
%<a href="matlab:help findcmt">findcmt</a>                        - Returns Global CMT that is closest to input time & position
%<a href="matlab:help findcmts">findcmts</a>                       - Returns Global CMTs in specified time & position ranges
%<a href="matlab:help hrv2ar">hrv2ar</a>                         - Convert moment tensor from Harvard from to Aki & Richards form
%<a href="matlab:help mapcmts">mapcmts</a>                        - Plots GlobalCMT moment tensors on a M_Map map
%<a href="matlab:help mo2hd">mo2hd</a>                          - Returns GlobalCMT half-duration estimate based on scalar moment
%<a href="matlab:help momentmag">momentmag</a>                      - Returns the moment magnitude for a moment tensor
%<a href="matlab:help mt2sdr">mt2sdr</a>                         - Converts moment tensors to strike-dip-rake using decomposition
%<a href="matlab:help mt2tpb">mt2tpb</a>                         - Returns the principal axes of moment tensors
%<a href="matlab:help mt_change">mt_change</a>                      - Alter moment tensor format (array or struct
%<a href="matlab:help mt_check">mt_check</a>                       - Checks moment tensor data (array or struct)
%<a href="matlab:help mt_decomp">mt_decomp</a>                      - Decompose moment tensor(s)
%<a href="matlab:help mt_diag">mt_diag</a>                        - Returns diagonalized moment tensors & principal axes
%<a href="matlab:help mt_norm">mt_norm</a>                        - Returns moment tensors normalized by their scalar moment
%<a href="matlab:help mt_undiag">mt_undiag</a>                      - De-diagonalizes moment tensor(s) into Harvard orientation
%<a href="matlab:help neu2vpa">neu2vpa</a>                        - Convert vector in North/East/Up to Value/Plunge/Azimuth
%<a href="matlab:help nodallines">nodallines</a>                     - Returns nodal lines for a focal mechanism
%<a href="matlab:help norm2strikedip">norm2strikedip</a>                 - Returns strike & dip given the normal to a fault plane
%<a href="matlab:help normslip2sdr">normslip2sdr</a>                   - Convert normal & slip vectors to strike, dip & rake angles
%<a href="matlab:help plotmt">plotmt</a>                         - Plot moment tensor(s) as beach ball(s)
%<a href="matlab:help plotmt3">plotmt3</a>                        - Plot moment tensor(s) as 3D beach ball(s)
%<a href="matlab:help radpat">radpat</a>                         - Calculates moment tensor radiation pattern
%<a href="matlab:help scalarmoment">scalarmoment</a>                   - Returns the scalar moment of a moment tensor
%<a href="matlab:help sdr2mt">sdr2mt</a>                         - Convert strike-dip-rake to moment tensor
%<a href="matlab:help sdr2null">sdr2null</a>                       - Calculates the null vector for a given strike, dip and rake
%<a href="matlab:help sdr2slip">sdr2slip</a>                       - Calculates the slip vector for a given strike, dip and rake
%<a href="matlab:help sdr2tpb">sdr2tpb</a>                        - Returns the principal axes of a focal mechanism
%<a href="matlab:help strikedip2norm">strikedip2norm</a>                 - Returns the normal vector to a fault plane
%<a href="matlab:help tpb2mt">tpb2mt</a>                         - Returns moment tensors given their principal axes
%<a href="matlab:help tpb2sdr">tpb2sdr</a>                        - Returns the strike, dip & rake given principal axes
%<a href="matlab:help vpa2neu">vpa2neu</a>                        - Convert vector in Value/Plunge/Azimuth to North/East/Up
%
% <a href="matlab:help seizmo">SEIZMO - Passive Seismology Toolbox</a>
end