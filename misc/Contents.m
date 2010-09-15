% Seismology Toolbox - misc
% Version 0.6.0-r165 Fuji 15-Sept-2010
%
% Miscellaneous Support functions
%AXEXPAND             - Expands axes by factor
%AXMOVE               - Moves a set of axes by the specified amount
%AXSTRETCH            - Stretches a set of axes as a group (resizing+moving)
%BUTTORD2             - Butterworth filter order selection. (Honors passband corners)
%CIRCLE               - Returns points on a circle in cartesian space
%CLONEFIGURE          - Makes a clone of a figure
%DDENDROGRAM          - Generate dendrogram plot (with some extra color options)
%EXPANDSCALARS        - Expands scalars to match size of array inputs
%FIG2PRINT            - Adjusts figures to be as printed (aka print preview-ish)
%FILTER_BANK          - Makes a set of narrow-band bandpass filters
%FISHER               - Converts correlation coefficients to the Z statistic
%GAUSSIANTF           - Returns a gaussian time function
%GETAPPLICATION       - Returns application running this script and its version
%GETSUBFIELD          - Get substructure field contents
%GETWORDS             - Returns a cell array of words from a string
%IFISHER              - Converts Z statistics to correlation coefficients
%IIRDESIGN            - Designs an iir filter with the given constraints
%INVERTCOLOR          - Inverts colors given as rgb triplet or as short/long names
%INTERPDC1            - 1D interpolation (table lookup) with discontinuity support
%ISEQUALNUMELORSCALAR - True if all inputs have equal numel or are scalar
%ISEQUALSIZEORSCALAR  - True if all input arrays are equal size or scalar
%ISORTHOGONAL         - TRUE if orientations are orthogonal
%ISPARALLEL           - TRUE if orientations are parallel
%ISSTRING             - True for a string (row vector) of characters
%JOINWORDS            - Combines a cellstr into a space-separated string
%LIND                 - Returns a linear indices matrix
%LBWH2LRBT            - Convert left-bottom-width-height to left-right-bottom-top
%LRBT2LBWH            - Convert left-right-bottom-top to left-bottom-width-height
%LTI2SUB              - Square matrix lower triangle linear indices to subscripts
%MAKESUBPLOTS         - Makes subplots in current figure returning axes handles
%MAPLOCATIONS         - Map station/event locations
%MAT2VEC              - Converts matrices to column vectors
%MATCHSORT            - Replicates a sort operation using the returned permutation indices
%MCXC                 - Multi-channel cross correlation with built-in peak picker
%MOVEKIDS             - Moves the specified children plot objects
%NAME2RGB             - Converts short/long color names to RGB triplets
%NANMEAN              - Return mean excluding NaNs
%NANVARIANCE          - Return variance excluding NaNs
%NATIVEBYTEORDER      - Returns native endianness of present platform
%NDSQUAREFORM         - Reshapes between an n-d distance matrix and "vector"
%NEXTPOW2N            - Returns the next higher power of 2 for all array elements
%NOCOLORBARS          - Removes colorbars associated with the specified axes
%NOINVERT             - Turns off hardcopy black/white inversion
%NOLABELS             - Removes tick and axis labels
%NOTICKS              - Removes ticks and tick labels from axes
%NOTITLES             - Removes titles from specified axes
%ONEFILELIST          - Compiles multiple filelists into one
%PPDCVAL              - Evaluate piecewise polynomial (w/ discontinuity support)
%PRINT_TIME_LEFT      - Ascii progress bar
%READCSV              - Read in .csv formatted file as a structure
%READTXT              - Reads in a text file as a single string
%RRAT                 - Relative rational approximation
%SETFONTS             - Sets font props for all text objects in the specified axes
%SLIDINGAVG           - Returns sliding-window average of data
%SNR2MAXPHASEERROR    - Returns maximum narrow-band phase error based on SNR 
%SORT2LI              - Transforms permutation indices from sort to linear indices
%SSIDX                - Scalar struct database indexing
%STAR69               - Returns who called the current function
%STRNLEN              - Pad/truncate char/cellstr array to n character columns
%SUB2LTI              - Square matrix lower triangle linear indices from subscripts
%SUB2UTI              - Square matrix upper triangle linear indices from subscripts
%SUBMAT               - Returns a submatrix reduced along indicated dimensions
%SUBMAT_EVAL          - Returns a submatrix using eval
%SUPERCOLORBAR        - Makes a colorbar spanning multiple axes
%SUPERTITLE           - Makes a title spanning multiple axes
%SUPERXLABEL          - Makes an x-axis label spanning multiple axes
%SUPERYLABEL          - Makes a y-axis label spanning multiple axes
%SWAP                 - Swap values
%TAPERFUN             - Returns a taper as specified
%TRIANGLETF           - Returns a triangle time function
%UNSORT               - Undoes a sort operation using the returned sort indices
%UNIXCOMPRESSAVI      - Compress an AVI file in Unix with "MEncoder"
%UTI2SUB              - Square matrix upper triangle linear indices to subscripts
%VECNORM              - Returns vector norms
%WRITECSV             - Write out .csv formatted file from a structure
%XDIR                 - Cross-app compatible directory listing with recursion

