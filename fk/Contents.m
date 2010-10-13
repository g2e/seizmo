% Seismology Toolbox - fk
% Version 0.6.0-r170 Gunnbjørnfjeld 11-Oct-2010
%
% Frequency-Wavenumber analysis functions
%CHKFKARFSTRUCT    - Validate if a struct is as defined by FKARF
%CHKFKSTRUCT       - Validate if a struct is as defined by FKMAP/VOLUME/4D
%CHKGEOFKARFSTRUCT - Validate if a struct is as defined by GEOFKARF
%CHKGEOFKSTRUCT    - Validate if a struct is as defined by GEOFK functions
%FK4D              - Returns beamformer volumes in frequency-wavenumber-time space
%FKARF             - Returns the fk array response function for a seismic array
%FKCART2POL        - Converts a cartesian space based fk grid to polar space
%FKCORR            - Zero lag correlation coefficient between fk spectra
%FKDBINFO          - Returns the min/median/max dB for a FK struct
%FKFRAMESLIDE      - Slides through a sequence of fk maps plotting each one
%FKFREQSLIDE       - Slides through a fk volume plotting each frequency
%FKMAP             - Returns beamformer map in frequency-wavenumber space
%FKSUBVOL          - Extracts a frequency-based subvolume of a fk volume
%FKVOL2MAP         - Converts a fk volume to a fk map
%FKVOLUME          - Returns beamformer volume in frequency-wavenumber space
%FKXCHORZVOLUME    - Returns frequency-wavenumber space for horz. xc data
%FKXCVOLUME        - Returns energy map in frequency-wavenumber space for xc data
%GEOFKARF          - Returns the geofk array response function for a seismic array
%GEOFKARF2MAP      - Converts a geofk ARF volume to a geofk ARF map
%GEOFKARFCORR      - Computes correlation with ARFs in the parameter space
%GEOFKARFSLOWSLIDE - Slides through the slownesses of a geofkarf volume
%GEOFKCORR         - Zero lag correlation coefficient between geofk spectra
%GEOFKDBINFO       - Returns the min/median/max dB for a geofk struct
%GEOFKFRAMESLIDE   - Slides through a set of geofk maps plotting each one
%GEOFKFREQSLIDE    - Slides through a geofk volume plotting each frequency
%GEOFKFREQSLIDEX   - Marks peak while sliding through geofk frequencies
%GEOFKSLOWSLIDE    - Slides through a geofk volume plotting each slowness
%GEOFKSUBARF       - Extracts a subARF of a geofkarf volume
%GEOFKSUBVOL       - Extracts a subvolume of a geofk volume
%GEOFKVOL2MAP      - Converts a geofk volume to a geofk map
%GEOFKXCHORZVOLUME - Geographic FK beamforming of horizontals
%GEOFKXCVOLUME     - Geographic FK beamforming
%ISFKARFSTRUCT     - True if input is a fk struct as defined by FKARF
%ISFKSTRUCT        - True if input is a fk struct as defined by FKMAP/VOLUME/4D
%ISGEOFKARFSTRUCT  - True if input is a fk struct as defined by GEOFKARF
%ISGEOFKSTRUCT     - True if input is a fk struct as defined by GEOFK
%KXY2SLOWBAZ       - Converts wavenumbers in x & y to slowness and back-azimuth
%PLOTFKARF         - Plots an fk array response function
%PLOTFKAZIFREQ     - Plots beam intensity as a function of azimuth & frequency
%PLOTFKMAP         - Plots the frequency-wavenumber output from FKMAP
%PLOTGEOFKARF      - Plots a geofk array response
%PLOTGEOFKMAP      - Plots frequency-slowness-position response info
%SLOWBAZ2KXY       - Converts slowness and back-azimuth to wavenumbers in x & y
%SLOWNESS2DEG      - Converts horizontal slowness (s/deg) to distance (deg)
%SNYQUIST          - Returns the nyquist slowness for an array
%UPDATEFKMAP       - Quickly updates an existing fk plot with a new map
%UPDATEGEOFKARF    - Quickly updates an existing geofkarf plot with a new arf
%UPDATEGEOFKMAP    - Quickly updates an existing geofk plot with a new map
