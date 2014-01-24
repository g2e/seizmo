function []=sz_toc_models()
% 1D/3D Earth Model Functions
%<a href="matlab:help ak135">ak135</a>               - Returns the AK135 Earth model
%<a href="matlab:help available_1dmodels">available_1dmodels</a>  - Returns available 1D Earth model functions
%<a href="matlab:help available_3dmodels">available_3dmodels</a>  - Returns available 3D Earth mantle model functions
%<a href="matlab:help chk1dmodel">chk1dmodel</a>          - Validate if a struct is as defined by PREM/AK135/IASP91
%<a href="matlab:help cmb_1dmodel_library">cmb_1dmodel_library</a> - Generates a library of CMB 1D Earth models
%<a href="matlab:help flatten_1dmodel">flatten_1dmodel</a>     - Flattens a 1D Earth model
%<a href="matlab:help getcrust">getcrust</a>           - Returns crustal info at specified location(s)
%<a href="matlab:help get_scripps_value">get_scripps_value</a>   - Returns dv% value for a Scripps mantle model
%<a href="matlab:help iasp91">iasp91</a>              - Returns the IASP91 Earth Model
%<a href="matlab:help mantledv">mantledv</a>            - Returns the seismic velocity deviation for a mantle model
%<a href="matlab:help mantlemap">mantlemap</a>           - 3D mantle model map (aka depth slice)
%<a href="matlab:help mantleprofile">mantleprofile</a>       - 3D mantle model profile (aka radial slice)
%<a href="matlab:help perturb_1dmodel">perturb_1dmodel</a>     - Perturbs 1D Earth models
%<a href="matlab:help plot1dmodel">plot1dmodel</a>         - Plots 1D model properties
%<a href="matlab:help prem">prem</a>                - Returns the PREM Earth Model
%<a href="matlab:help prem_perfect">prem_perfect</a>        - Returns the PREM model
%<a href="matlab:help prem2_perfect">prem2_perfect</a>       - Returns the PREM2 model
%<a href="matlab:help qkqu2qp">qkqu2qp</a>             - Calculate corresponding Qalpha for given Qk, Qu, Vp, Vs
%<a href="matlab:help ql6">ql6</a>                 - Quality factor model of the Earth by Durek & Ekstrom 1996
%<a href="matlab:help qlm9">qlm9</a>                - Quality factor model of the Earth by Lawrence & Wysession 2006
%<a href="matlab:help qpqs2qk">qpqs2qk</a>             - Calculate corresponding Qkappa for given Qp, Qs, Vp, Vs
%<a href="matlab:help readcrust2">readcrust2</a>          - Reads Crust2.0 files, putting them into a struct
%<a href="matlab:help readcrust10">readcrust10</a>          - Reads Crust1.0 files, putting them into a struct
%<a href="matlab:help read_scripps_binary">read_scripps_binary</a> - Reads a Scripps mantle model unformatted binary
%<a href="matlab:help scripps2matlab">scripps2matlab</a>      - Creates a Matlab struct from a Scripps mantle model
%<a href="matlab:help write_1dmodel_nd">write_1dmodel_nd</a>    - Writes a 1D model struct in .nd format
%
% <a href="matlab:help seizmo">SEIZMO - Passive Seismology Toolbox</a>
end
