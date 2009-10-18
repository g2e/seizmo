% Seismology Toolbox for Matlab and Octave
% Version 0.6.0-r103 Ararat 17-Oct-2009
%
% Miscellaneous Support functions
%FILTER_BANK         - Makes a set of narrow-band bandpass filters
%FISHER              - Converts correlation coefficients to the Z statistic
%GAUSSIANTF          - Returns a gaussian time function
%GETAPPLICATION      - Returns application running this script and its version
%GETSUBFIELD         - Get substructure field contents
%GETWORDS            - Returns a cell array of words from a string
%IFISHER             - Converts Z statistics to correlation coefficients
%JOINWORDS           - Combines a cellstr into a space-separated string
%LATMOD              - Returns a latitude modulus
%LONMOD              - Returns a longitude modulus
%LTI2SUB             - Square matrix lower triangle linear indices to subscripts
%MATCHSORT           - Replicates a sort operation using the returned permutation indices
%MCXC                - Multi-channel cross correlation with built-in peak picker
%NANVARIANCE         - Return variance excluding NaNs
%NATIVEBYTEORDER     - Returns native endianness of present platform
%NDSQUAREFORM        - Reshapes between an n-d distance matrix and "vector"
%NEXTPOW2N           - Returns the next higher power of 2 for all array elements
%ONEFILELIST         - Compiles multiple filelists into one
%PRINT_TIME_LEFT     - Ascii progress bar
%READCSV             - Read in .csv formatted file as a structure
%SLIDINGAVG          - Returns sliding-window average of data
%SNR2MAXPHASEERROR   - Returns maximum narrow-band phase error based on SNR 
%SORT2LI             - Transforms permutation indices from sort to linear indices
%STRNLEN             - Pad/truncate char/cellstr array to n character columns
%SUB2LTI             - Square matrix lower triangle linear indices from subscripts
%SUB2UTI             - Square matrix upper triangle linear indices from subscripts
%SUBMAT              - Returns a submatrix reduced along indicated dimensions
%SUBMAT_EVAL         - Returns a submatrix using eval
%SWAP                - Swap values
%TAPERFUN            - Returns a taper as specified
%TRIANGLETF          - Returns a triangle time function
%UNSORT              - Undoes a sort operation using the returned sort indices
%UTI2SUB             - Square matrix upper triangle linear indices to subscripts
%WRITECSV            - Write out .csv formatted file from a structure
%XDIR                - Cross-app compatible directory listing with recursion

