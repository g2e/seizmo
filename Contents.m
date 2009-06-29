% Seismology Toolbox for Matlab and Octave
% Version 0.6.0-r80 Ararat 29-Jun-2009
%
% Alphabetical list of SEIZMO functions
%ADD                 - Add a constant to SEIZMO records
%ADDRECORDS          - Add SEIZMO records
%AMPH2RLIM           - Convert SEIZMO spectral records from AMPH to RLIM
%BINOPERR            - Controls behavior of SEIZMO binary functions
%BSEIZMO             - Arranges xy data into a SEIZMO data structure
%CHANGEBYTEORDER     - Change the byteorder of SEIZMO data records
%CHANGECLASS         - Change SEIZMO data storage (in memory)
%CHANGEHEADER        - Change SEIZMO data header values
%CHANGENAME          - Change the filename of SEIZMO records
%CHANGEPATH          - Change the filepath of SEIZMO records
%CHECKHEADER         - Check and fix header values of SEIZMO records
%COPYHEADER          - Copy one record's header to all records
%CORRELATE           - Compute cross correlograms of SEIZMO data records
%CUT                 - Cut a window out of SEIZMO records
%DELETERECORDS       - Deletes indicated records from SEIZMO data structure
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
%GETMEDIAN           - Returns median of each SEIZMO record
%GETNCMP             - Return the number of dependent components for SEIZMO records
%GETNETWORKIDX       - Returns index array for separating dataset into networks
%GETNORM             - Return normalizers for SEIZMO records
%GETPOLYNOMIAL       - Get polynomial fit to SEIZMO records
%GETSPECTRALCMP      - Returns the indicated portion of spectral records
%GETSTATIONIDX       - Returns index array for separating dataset into stations
%GETSTREAMIDX        - Returns index array for separating dataset into streams
%HILBRT              - Return Hilbert transform of SEIZMO records
%IDFT                - Performs an inverse discrete fourier transform on SEIZMO records
%IIRFILTER           - Apply an IIR filter to SEIZMO data records
%INTEGRATE           - Integrate SEIZMO records
%INTERPOLATE         - Interpolate SEIZMO records to a new samplerate
%JOINRECORDS         - Join SEIZMO records into multiple-component record(s)
%KEEPAM              - Returns the amplitude component of spectral records
%KEEPIM              - Returns the imaginary component of spectral records
%KEEPPH              - Returns phase component of spectral records
%KEEPRECORDS         - Keeps indicated records in SEIZMO data structure
%KEEPRL              - Returns the real component of spectral records
%LISTHEADER          - List SEIZMO data headers
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
%PREWHITEN           - Prewhiten SEIZMO data records for spectral operations
%QUICKSNR            - Quick estimation of SNR for SEIZMO records
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
%SEIZMOFUN           - Apply function to SEIZMO records
%SELECTRECORDS       - Select or delete SEIZMO data records graphically
%SLIDINGABSMEAN      - Returns sliding-window absolute-mean of SEIZMO records
%SLIDINGFUN          - Apply a sliding window function to SEIZMO records
%SLIDINGRMS          - Returns sliding-window root-mean-square of SEIZMO records
%SORTBYFIELD         - Sort SEIZMO records by a header or SEIZMO struct field
%SPLITRECORDS        - Split up components into separate records
%SQUISH              - Downsample SEIZMO records by an integer factor
%STRETCH             - Upsample SEIZMO records by an integer factor
%SUBTRACT            - Subtract a constant from SEIZMO records
%SUBTRACTRECORDS     - Subtract SEIZMO records
%SYNCHRONIZE         - Synchronizes the reference times of SEIZMO records
%SYNCRATES           - Resample SEIZMO records to a common sample rate
%TAPER               - Taper SEIZMO records
%TIMESHIFT           - Shift timing of SEIZMO records
%UNPREWHITEN         - Undo prewhitening of SEIZMO data records
%WHITEN              - Spectral whitening/normalization of SEIZMO data records
%WRITEHEADER         - Write SEIZMO data header info to datafiles
%WRITESEIZMO         - Write SEIZMO records to datafile
%
%
% Time functions
%CAL2DOY             - Convert Year & Month & Day of Month to Year & Day of Year
%DOY2CAL             - Convert Year & Day of Year to Year & Month & Day of Month
%FIXDATES            - Assures calendar or day of year dates are in proper range
%FIXTIMES            - Cleans up times so they are in their proper ranges
%GETLEAPSECONDS      - Returns leapsecond info with caching for multiple calls
%GREGORIAN2MODSERIAL - Convert Gregorian dates to modified serial dates
%GREGORIAN2SERIAL    - Convert Gregorian dates to serial dates
%ISLEAPYEAR          - True if year is a leap year
%LEAPSECONDS         - Returns date and time strings of each leap second
%LEAPSINDAY          - Returns the number of leap seconds in the given dates
%MODSERIAL2GREGORIAN - Convert modified serial dates to Gregorian dates
%SERIAL2GREGORIAN    - Convert serial dates to Gregorian dates
%TAI2UTC             - Convert TAI time to UTC time
%TIMEDIFF            - Return number of seconds between times
%TOTALLEAPS          - Returns the accumulated leap seconds for the given dates
%UTC2TAI             - Convert UTC time to TAI time
%
%
% Position functions
%GEOCENTRIC2GEODETICLAT - Convert latitudes from geocentric to geodetic
%GEOCENTRIC2GEODETIC    - Converts coordinates from geocentric to geodetic
%GEOCENTRIC2XYZ         - Converts coordinates from geocentric to cartesian
%GEODETIC2GEOCENTRICLAT - Convert latitudes from geodetic to geocentric
%GEODETIC2GEOCENTRIC    - Converts coordinates from geodetic to geocentric
%GEODETIC2XYZ           - Converts coordinates from geodetic to cartesian
%GEODETICLAT2RADIUS     - Returns the radius at a geodetic latitude
%HAVERSINE              - Returns distance between 2 points using the Haversine formula
%SPHERICALFWD           - Finds a point on a sphere relative to another point
%SPHERICALINV           - Return distance and azimuth between 2 locations on sphere
%VINCENTYFWD            - Find destination point on an ellipsoid relative to a point
%VINCENTYINV            - Find distance and azimuth between 2 locations on ellipsoid
%XYZ2GEOCENTRIC         - Converts coordinates from cartesian to geocentric
%XYZ2GEODETIC           - Converts coordinates from cartesian to geodetic
%
%
% Miscellaneous Support functions
%CMOD                - Returns a centered modulus
%FILTER_BANK         - Makes a set of narrow-band bandpass filters
%FISHER              - Converts correlation coefficients to the Z statistic
%GETAPPLICATION      - Returns application running this script and its version
%GETWORDS            - Returns a cell array of words from a string
%IFISHER             - Converts Z statistics to correlation coefficients
%LTI2SUB             - Square matrix lower triangle linear index to subscripts
%MATCHSORT           - Replicates a sort operation using the returned permutation indices
%MCXC                - Multi-channel cross correlation with built-in peak picker
%NANMEAN             - Return mean excluding NaNs
%NANVAR              - Return variance excluding NaNs
%NATIVEBYTEORDER     - Returns native endianness of present platform
%NDSQUAREFORM        - Reshapes a multi-page distance matrix between square and triangle vector forms
%NEXTPOW2N           - Returns the next higher power of 2 for all array elements
%ONEFILELIST         - Compiles multiple filelists into one
%READSACPZ           - Reads in a SAC PoleZero file
%SAWMOD              - Returns a sawtooth modulus
%SLIDINGAVG          - Returns sliding-window average of data
%SORT2LI             - Transforms permutation indices from sort to linear indices
%STRNLEN             - Pad/truncate char/cellstr array to n character columns
%SUB2LTI             - Square matrix lower triangle linear index from subscripts
%SUB2UTI             - Square matrix upper triangle linear index from subscripts
%SUBMAT_EVAL         - Returns a submatrix reduced along indicated dimensions
%SUBMAT              - Returns a submatrix reduced along indicated dimensions
%SWAP                - Swap values
%UNIQSORT            - A fast unique sort for numeric vectors
%UTI2SUB             - Square matrix upper triangle linear index to subscripts
%WRITESACPZ          - Writes out a SAC PoleZero file
%XDIR                - Directory listing with recursion

