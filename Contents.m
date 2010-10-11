% Seismology Toolbox - seizmo
% Version 0.6.0-r170 Gunnbjørnfjeld 11-Oct-2010
%
% Alphabetical list of SEIZMO functions
%ADD                  - Add a constant to SEIZMO records
%ADDARRIVALS          - Adds the indicated phases to the SEIZMO data header
%ADDRECORDS           - Add SEIZMO records
%ATTACH               - Attach data to SEIZMO records
%AMPH2RLIM            - Convert SEIZMO spectral records from AMPH to RLIM
%AUTOWINDOW           - Automatic windowing of SEIZMO records
%BINOPERR             - Controls behavior of SEIZMO binary functions
%BSEIZMO              - Arranges xy data into a SEIZMO data structure
%CHANGEBYTEORDER      - Change the byteorder of SEIZMO data records
%CHANGECLASS          - Change SEIZMO data storage (in memory)
%CHANGEHEADER         - Change SEIZMO data header values
%CHANGENAME           - Change the filename of SEIZMO records
%CHANGEPATH           - Change the filepath of SEIZMO records
%CHECKHEADER          - Check and fix header values of SEIZMO records
%CHECKOPERR           - Controls behavior of SEIZMO CHECKHEADER function
%COMPAREHEADER        - List SEIZMO headers side-by-side for easy comparison
%CONVOLVE             - Convolve SEIZMO records with a time function
%COPYHEADER           - Copy one record's header to all records
%CORRELATE            - Compute cross correlograms of SEIZMO data records
%CUT                  - Cut a window out of SEIZMO records
%DECONVOLVE           - Spectrally deconvolve a time function from SEIZMO records
%DELETERECORDS        - Deletes indicated records from SEIZMO data structure
%DETACH               - Detach data from SEIZMO records
%DFT                  - Performs a discrete fourier transform on SEIZMO data records
%DIFFERENTIATE        - Differentiate SEIZMO records
%DIVIDE               - Divide SEIZMO records by a constant
%DIVIDEOMEGA          - Integrate SEIZMO records in the frequency domain
%DIVIDERECORDS        - Divide SEIZMO records
%ENVELOPE             - Return envelopes of SEIZMO records
%FIXDELTA             - Fix sample spacing for SEIZMO records
%GENNAME              - Generic filenames for SEIZMO records
%GETARRIVAL           - Returns stored phase arrival time from SEIZMO data header
%GETCOMPONENTIDX      - Returns index array to separate dataset into components
%GETENUMDESC          - Get enum string description from SEIZMO header enum field
%GETENUMID            - Get enum string id from SEIZMO data header enum field
%GETHEADER            - Get SEIZMO data header values
%GETLGC               - Get logical string from SEIZMO logical header field
%GETNETWORKIDX        - Returns index array for separating dataset into networks
%GETNORM              - Return normalizers for SEIZMO records
%GETPOLYNOMIAL        - Get polynomial fit to SEIZMO records
%GETSITEIDX           - Returns index array for separating dataset by sites
%GETSPECTRALCMP       - Returns the indicated portion of spectral records
%GETSTATIONIDX        - Returns index array for separating dataset into stations
%GETSTREAMIDX         - Returns index array for separating dataset into streams
%GETVALUEFUN          - Applies a function to SEIZMO records and returns the value
%HILBRT               - Return Hilbert transform of SEIZMO records
%HORZPAIRS            - Returns indice arrays for pairing horizontal SEIZMO records
%IDFT                 - Performs an inverse discrete fourier transform on SEIZMO records
%IIRFILTER            - Apply an IIR filter to SEIZMO data records
%INSTANTFREQ          - Returns estimated instantaneous frequency of SEIZMO records
%INSTANTPHASE         - Return instantaneous phase of SEIZMO records
%INTEGRATE            - Integrate SEIZMO records
%INTERPOLATE          - Interpolate SEIZMO records to a new samplerate
%JOINRECORDS          - Join SEIZMO records into multiple-component record(s)
%KEEPAM               - Returns the amplitude component of spectral records
%KEEPIM               - Returns the imaginary component of spectral records
%KEEPPH               - Returns phase component of spectral records
%KEEPRECORDS          - Keeps indicated records in SEIZMO data structure
%KEEPRL               - Returns the real component of spectral records
%LISTHEADER           - List SEIZMO data headers
%MAKEUNEVEN           - Makes evenly-sampled SEIZMO records look like uneven ones
%MAPEVENTGRID         - Draws range/azimuth grid (relative to an event) on a map
%MAPSTATIONS          - Map station/earthquake locations of SEIZMO records
%MAT2RECORDS          - Distributes a record matrix back into a SEIZMO struct
%MERGE                - Merge SEIZMO records
%MIRRORFLIP           - Returns a mirror-flip of SEIZMO records
%MULTIPLY             - Multiply SEIZMO records by a constant
%MULTIPLYOMEGA        - Differentiate SEIZMO records in the frequency domain
%MULTIPLYRECORDS      - Multiply SEIZMO records
%NORMALIZE            - Normalizes SEIZMO records
%PLAYSEIZMO           - Plays SEIZMO records as sound
%PREWHITEN            - Prewhiten SEIZMO data records for spectral operations
%QUERYHEADER          - List SEIZMO headers side-by-side for easy comparison
%QUICKSNR             - Quick estimation of SNR for SEIZMO records
%RAISE                - Raises the dependent component of SEIZMO records by a given power
%READDATA             - Read SEIZMO data from datafiles
%READDATAWINDOW       - Read data window into SEIZMO data structure
%READHEADER           - Read header info from SEIZMO datafiles
%READSEIZMO           - Read datafiles into SEIZMO data structure
%RECORDFUN            - Perform basic function between SEIZMO records
%RECORDS2MAT          - Combines SEIZMO data records into a single numeric matrix
%REMOVEDEADRECORDS    - Removes constant SEIZMO records
%REMOVEDUPLICATES     - Remove duplicate SEIZMO records
%REMOVEMEAN           - Remove mean from SEIZMO records
%REMOVEPOLYNOMIAL     - Remove polynomial trend from SEIZMO records
%REMOVETREND          - Remove linear trend from SEIZMO records
%REVERSE              - Reverse SEIZMO records
%REVERSE_CORRELATIONS - Reverses correlations (switch master & slave)
%RLIM2AMPH            - Convert SEIZMO spectral records from RLIM to AMPH
%ROTATE               - Rotates SEIZMO records that are horizontal pairs
%SEIZMO2WAV           - Writes SEIZMO records as WAV files
%SEIZMOFUN            - Apply function to SEIZMO records
%SEIZMORESET          - Resets saved SEIZMO settings
%SEIZMOVERBOSE        - Turn verbose SEIZMO output on (TRUE) or off (FALSE)
%SETASIDAY            - Sets the reference times of records to the start of the day
%SLIDINGABSMEAN       - Returns sliding-window absolute-mean of SEIZMO records
%SLIDINGFUN           - Apply a sliding window function to SEIZMO records
%SLIDINGMEAN          - Returns sliding-window mean of SEIZMO records
%SLIDINGRMS           - Returns sliding-window root-mean-square of SEIZMO records
%SORTBYFIELD          - Sort SEIZMO records by a header or SEIZMO struct field
%SPLITPAD             - Splits and zero pads SEIZMO data records
%SPLITRECORDS         - Split up components into separate records
%SQUISH               - Downsample SEIZMO records by an integer factor
%STACK                - Stacks SEIZMO records
%STFT                 - Short Time Fourier Transform (aka Sliding Spectra) of SEIZMO data
%STRETCH              - Upsample SEIZMO records by an integer factor
%SUBTRACT             - Subtract a constant from SEIZMO records
%SUBTRACTRECORDS      - Subtract SEIZMO records
%SYNCHRONIZE          - Synchronizes the reference times of SEIZMO records
%SYNCRATES            - Resample SEIZMO records to a common sample rate
%TAPER                - Taper SEIZMO records
%TIMESHIFT            - Shift timing of SEIZMO records
%UNPREWHITEN          - Undo prewhitening of SEIZMO data records
%UNWRAPPHASE          - Unwraps the phase of SEIZMO records
%WHITEN               - Spectral whitening/normalization of SEIZMO data records
%WRITEHEADER          - Write SEIZMO data header info to datafiles
%WRITESEIZMO          - Write SEIZMO records to datafile
%
% Sub-Toolboxes
%<a href="matlab:help cmap">cmap</a>                - colormaps
%<a href="matlab:help cmt">cmt</a>                 - various moment tensor functions
%<a href="matlab:help event">event</a>               - functions related to earthquake info
%<a href="matlab:help fixes">fixes</a>               - functions to fix headers
%<a href="matlab:help fk">fk</a>                  - frequency-wavenumber beamforming functions
%<a href="matlab:help gui">gui</a>                 - simple gui functions
%<a href="matlab:help invert">invert</a>              - simple inversion functions
%<a href="matlab:help mattaup">mattaup</a>             - TauP functions
%<a href="matlab:help misc">misc</a>                - useful miscellaneous functions
%<a href="matlab:help models">models</a>              - 1D/3D Earth model functions
%<a href="matlab:help plotting">plotting</a>            - SEIZMO plotting routines
%<a href="matlab:help position">position</a>            - global positioning functions
%<a href="matlab:help response">response</a>            - instrument response functions
%<a href="matlab:help seizmo_internal">seizmo_internal</a>     - internal functions
%<a href="matlab:help shortnames">shortnames</a>          - shortcut names for functions
%<a href="matlab:help sphpoly">sphpoly</a>             - Spherical Polygon functions
%<a href="matlab:help synth">synth</a>               - functions for synthetics
%<a href="matlab:help time">time</a>                - absolute time functions
%<a href="matlab:help topo">topo</a>                - topography functions
%<a href="matlab:help tomo">tomo</a>                - simple tomography functions
%<a href="matlab:help tpw">tpw</a>                 - two plane wave functions
%<a href="matlab:help ttcorrect">ttcorrect</a>           - travel time correction functions
%<a href="matlab:help ww3">ww3</a>                 - WaveWatch III functions
%<a href="matlab:help xcalign">xcalign</a>             - cross correlation signal alignment functions
%
%
% Found a bug?  Something not quite right?
% Email me: ggeuler@seismo.wustl.edu
%
% Garrett Euler, Ph.D. Candidate
% Department of Earth & Planetary Sciences
% Washington University in Saint Louis
