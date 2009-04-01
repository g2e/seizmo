% Seismology Toolbox for Matlab and Octave
% Version 0.6.0-r68 (r2007b) 01-Apr-2009
%
% Alphabetical list of functions
%ADD                 - Add a constant to SEIZMO records
%ADDRECORDS          - Add SEIZMO records
%AMPH2RLIM           - Convert SEIZMO spectral records from AMPH to RLIM
%BINOPERR            - Controls behavior of SEIZMO binary functions
%BSEIZMO             - Arranges xy data into a SEIZMO data structure
%CHANGEBYTEORDER     - Change the byteorder of SEIZMO data records
%CHANGECLASS         - Change SEIZMO data storage (in memory)
%CHANGEHEADER        - Change SEIZMO data header values
%CHECKHEADER         - Check and fix header values of SEIZMO records
%COMBINERECORDS      - Combines SEIZMO data records into a single numeric record matrix
%CUT                 - Cut a window out of SEIZMO records
%DFT                 - Performs a discrete fourier transform on SEIZMO data records
%DIFFERENTIATE       - Differentiate SEIZMO records
%DISTRIBUTERECORDS   - Distributes a record matrix back into a SEIZMO struct
%DIVIDE              - Divide SEIZMO records by a constant
%DIVIDEOMEGA         - Integrate SEIZMO records in the frequency domain
%DIVIDERECORDS       - Divide SEIZMO records
%ENVELOPE            - Return envelopes of SEIZMO records
%FIXDELTA            - Fix sample spacing for SEIZMO records
%GETARRIVAL          - Returns stored phase arrival time from SEIZMO data header
%GETENUMDESC         - Get enum string description from SEIZMO header enum field
%GETENUMID           - Get enum string id from SEIZMO data header enum field
%GETHEADER           - Get SEIZMO data header values
%GETLGC              - Get logical string from SEIZMO logical header field
%GETMEDIAN           - Returns median of each SEIZMO record
%GETNCMP             - Return the number of dependent components for SEIZMO records
%GETNORM             - Return normalizers for SEIZMO records
%HILBRT              - Return Hilbert transform of SEIZMO records
%IDFT                - Performs an inverse discrete fourier transform on SEIZMO records
%IIRFILTER           - Apply an IIR filter to SEIZMO data records
%INTEGRATE           - Integrate SEIZMO records
%INTERPOLATE         - Interpolate SEIZMO records to a new samplerate
%LISTHEADER          - List SEIZMO data headers
%MERGE               - Merge SEIZMO records
%MULTIPLY            - Multiply SEIZMO records by a constant
%MULTIPLYOMEGA       - Differentiate SEIZMO records in the frequency domain
%MULTIPLYRECORDS     - Multiply SEIZMO records
%NORMALIZE           - Normalizes SEIZMO records
%PLOT0               - Plot SEIZMO data records in an evenly spaced record section
%PLOT1               - Plot SEIZMO data records in individual subplots
%PLOT2               - Overlay plot of SEIZMO data records
%PLOTCONFIG          - Returns default configuration structure for SEIZMO plotting
%PLOTDENDRO          - Plots correlation linkage and seismograms
%QUICKSNR            - Quick estimation of SNR for SEIZMO records
%READDATA            - Read SEIZMO data from datafiles
%READDATAWINDOW      - Read data window into SEIZMO data structure
%READHEADER          - Read header info from SEIZMO datafiles
%READSEIZMO          - Read datafiles into SEIZMO data structure
%RECORDFUN           - Perform basic function between SEIZMO records
%RECORDSECTION       - Plots SEIZMO data records in a distance spaced record section
%REMOVEDEADRECORDS   - Removes constant SEIZMO records
%REMOVEDUPLICATES    - Remove duplicate SEIZMO records
%REMOVEMEAN          - Remove mean from SEIZMO records
%REMOVETREND         - Remove linear trend from SEIZMO records
%REVERSE             - Reverse SEIZMO records
%RLIM2AMPH           - Convert SEIZMO spectral records from RLIM to AMPH
%SEIZMOFUN           - Apply function to SEIZMO records
%SELECTRECORDS       - Select or delete SEIZMO data records graphically
%SLIDINGABSMEAN      - Returns sliding-window absolute-mean of SEIZMO records
%SLIDINGFUN          - Apply a sliding window function to SEIZMO records
%SLIDINGRMS          - Returns sliding-window root-mean-square of SEIZMO records
%SORTBYFIELD         - Sort SEIZMO records by a header or SEIZMO struct field
%SQUISH              - Downsample SEIZMO records by an integer factor
%STRETCH             - Upsample SEIZMO records by an integer factor
%SUBTRACT            - Subtract a constant from SEIZMO records
%SUBTRACTRECORDS     - Subtract SEIZMO records
%SYNCRATES           - Resample SEIZMO records to a common sample rate
%TAPER               - Taper SEIZMO records
%TIMESHIFT           - Shift timing of SEIZMO records
%WRITEHEADER         - Write SEIZMO data header info to datafiles
%WRITESEIZMO         - Write SEIZMO records to datafile

