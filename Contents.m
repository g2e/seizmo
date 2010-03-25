% Seismology Toolbox - seizmo
% Version 0.6.0-r120 Blackburn 20-Mar-2010
%
% Alphabetical list of SEIZMO functions
%ADD                 - Add a constant to SEIZMO records
%ADDARRIVALS         - Adds the indicated phases to the SEIZMO data header
%ADDRECORDS          - Add SEIZMO records
%ATTACH              - Attach data to SEIZMO records
%AMPH2RLIM           - Convert SEIZMO spectral records from AMPH to RLIM
%BINOPERR            - Controls behavior of SEIZMO binary functions
%BSEIZMO             - Arranges xy data into a SEIZMO data structure
%CHANGEBYTEORDER     - Change the byteorder of SEIZMO data records
%CHANGECLASS         - Change SEIZMO data storage (in memory)
%CHANGEHEADER        - Change SEIZMO data header values
%CHANGENAME          - Change the filename of SEIZMO records
%CHANGEPATH          - Change the filepath of SEIZMO records
%CHECKHEADER         - Check and fix header values of SEIZMO records
%CHECKOPERR          - Controls behavior of SEIZMO CHECKHEADER function
%COMPAREHEADER       - List SEIZMO headers side-by-side for easy comparison
%COMPAREHEADER2      - List SEIZMO headers side-by-side for easy comparison
%CONVOLVE            - Convolve SEIZMO records with a time function
%COPYHEADER          - Copy one record's header to all records
%CORRELATE           - Compute cross correlograms of SEIZMO data records
%CUT                 - Cut a window out of SEIZMO records
%DECONVOLVE          - Spectrally deconvolve a time function from SEIZMO records
%DELETERECORDS       - Deletes indicated records from SEIZMO data structure
%DETACH              - Detach data from SEIZMO records
%DFT                 - Performs a discrete fourier transform on SEIZMO data records
%DIFFERENTIATE       - Differentiate SEIZMO records
%DIVIDE              - Divide SEIZMO records by a constant
%DIVIDEOMEGA         - Integrate SEIZMO records in the frequency domain
%DIVIDERECORDS       - Divide SEIZMO records
%ENVELOPE            - Return envelopes of SEIZMO records
%FIXDELTA            - Fix sample spacing for SEIZMO records
%GETARRIVAL          - Returns stored phase arrival time from SEIZMO data header
%GETCOMPONENTIDX     - Returns index array to separate dataset into components
%GETENUMDESC         - Get enum string description from SEIZMO header enum field
%GETENUMID           - Get enum string id from SEIZMO data header enum field
%GETHEADER           - Get SEIZMO data header values
%GETLGC              - Get logical string from SEIZMO logical header field
%GETNETWORKIDX       - Returns index array for separating dataset into networks
%GETNORM             - Return normalizers for SEIZMO records
%GETPOLYNOMIAL       - Get polynomial fit to SEIZMO records
%GETSPECTRALCMP      - Returns the indicated portion of spectral records
%GETSTATIONIDX       - Returns index array for separating dataset into stations
%GETSTREAMIDX        - Returns index array for separating dataset into streams
%GETVALUEFUN         - Applies a function to SEIZMO records and returns the value
%HILBRT              - Return Hilbert transform of SEIZMO records
%HORZPAIRS           - Returns indice arrays for pairing horizontal SEIZMO records
%IDFT                - Performs an inverse discrete fourier transform on SEIZMO records
%IIRFILTER           - Apply an IIR filter to SEIZMO data records
%INSTANTPHASE        - Return instantaneous phase of SEIZMO records
%INTEGRATE           - Integrate SEIZMO records
%INTERPOLATE         - Interpolate SEIZMO records to a new samplerate
%JOINRECORDS         - Join SEIZMO records into multiple-component record(s)
%KEEPAM              - Returns the amplitude component of spectral records
%KEEPIM              - Returns the imaginary component of spectral records
%KEEPPH              - Returns phase component of spectral records
%KEEPRECORDS         - Keeps indicated records in SEIZMO data structure
%KEEPRL              - Returns the real component of spectral records
%LISTHEADER          - List SEIZMO data headers
%MAKEUNEVEN          - Makes evenly-sampled SEIZMO records look like uneven ones
%MAT2RECORDS         - Distributes a record matrix back into a SEIZMO struct
%MERGE               - Merge SEIZMO records
%MIRRORFLIP          - Returns a mirror-flip of SEIZMO records
%MULTIPLY            - Multiply SEIZMO records by a constant
%MULTIPLYOMEGA       - Differentiate SEIZMO records in the frequency domain
%MULTIPLYRECORDS     - Multiply SEIZMO records
%NORMALIZE           - Normalizes SEIZMO records
%PLOT0               - Plot SEIZMO data records in an evenly spaced record section
%PLOT1               - Plot SEIZMO data records in individual subplots
%PLOT2               - Overlay plot of SEIZMO data records
%PLOTCONFIG          - Returns default configuration structure for SEIZMO plotting
%PLOTDENDRO          - Plots correlation linkage and seismograms
%PLOTSTATIONS        - Plots station and earthquake locations of SEIZMO records
%PREWHITEN           - Prewhiten SEIZMO data records for spectral operations
%QUICKSNR            - Quick estimation of SNR for SEIZMO records
%RAISE               - Raises the dependent component of SEIZMO records by a given power
%READDATA            - Read SEIZMO data from datafiles
%READDATAWINDOW      - Read data window into SEIZMO data structure
%READHEADER          - Read header info from SEIZMO datafiles
%READSEIZMO          - Read datafiles into SEIZMO data structure
%RECORDFUN           - Perform basic function between SEIZMO records
%RECORDS2MAT         - Combines SEIZMO data records into a single numeric matrix
%RECORDSECTION       - Plots SEIZMO data records in a distance spaced record section
%REMOVEDEADRECORDS   - Removes constant SEIZMO records
%REMOVEDUPLICATES    - Remove duplicate SEIZMO records
%REMOVEMEAN          - Remove mean from SEIZMO records
%REMOVEPOLYNOMIAL    - Remove polynomial trend from SEIZMO records
%REMOVETREND         - Remove linear trend from SEIZMO records
%REVERSE             - Reverse SEIZMO records
%RLIM2AMPH           - Convert SEIZMO spectral records from RLIM to AMPH
%ROTATE              - Rotates SEIZMO records that are horizontal pairs
%SEIZMOFUN           - Apply function to SEIZMO records
%SEIZMORESET         - Resets saved SEIZMO settings
%SEIZMOVERBOSE       - Turn verbose SEIZMO output on (TRUE) or off (FALSE)
%SELECTRECORDS       - Select or delete SEIZMO data records graphically
%SETASIDAY           - Sets the reference times of records to the start of the day
%SLIDINGABSMEAN      - Returns sliding-window absolute-mean of SEIZMO records
%SLIDINGFUN          - Apply a sliding window function to SEIZMO records
%SLIDINGMEAN         - Returns sliding-window mean of SEIZMO records
%SLIDINGRMS          - Returns sliding-window root-mean-square of SEIZMO records
%SORTBYFIELD         - Sort SEIZMO records by a header or SEIZMO struct field
%SPLITPAD            - Splits and zero pads SEIZMO data records
%SPLITRECORDS        - Split up components into separate records
%SQUISH              - Downsample SEIZMO records by an integer factor
%STACK               - Stacks SEIZMO records
%STRETCH             - Upsample SEIZMO records by an integer factor
%SUBTRACT            - Subtract a constant from SEIZMO records
%SUBTRACTRECORDS     - Subtract SEIZMO records
%SYNCHRONIZE         - Synchronizes the reference times of SEIZMO records
%SYNCRATES           - Resample SEIZMO records to a common sample rate
%TAPER               - Taper SEIZMO records
%TIMESHIFT           - Shift timing of SEIZMO records
%UNPREWHITEN         - Undo prewhitening of SEIZMO data records
%UNWRAPPHASE         - Unwraps the phase of SEIZMO records
%WHITEN              - Spectral whitening/normalization of SEIZMO data records
%WRITEHEADER         - Write SEIZMO data header info to datafiles
%WRITESEIZMO         - Write SEIZMO records to datafile
%
% Sub-Toolboxes
%CMAP                - colormaps
%CMT                 - various moment tensor functions
%EVENT               - functions related to earthquake info
%FIXES               - functions to fix headers
%GUI                 - simple gui functions
%INVERT              - simple inversion functions
%MATTAUP             - TauP functions
%MISC                - useful miscellaneous functions
%POSITION            - global positioning functions
%RESPONSE            - instrument response functions
%SEIZMO_INTERNAL     - internal functions
%SHORTNAMES          - shortcut names for functions
%SYNTH               - functions for synthetics
%TIME                - absolute time functions
%TOPO                - topography functions
%TOMO                - simple tomography functions
%TPW                 - two plane wave functions

