function [conf]=defaultalignconf()
%DEFAULTALIGNCONF    Returns default configuration structure for align

% MAKE ID THE CURRENT DATE/TIME
id=clock; id=fix(id(6)+id(5)*1e2+id(4)*1e4+id(3)*1e6+id(2)*1e8+id(1)*1e10);

% SETTING CONFIGURE STRUCTURE DEFAULTS
conf=struct(...
...% MAIN
'INPUTCONFFILE','default',...   % CONFIGURATION FILENAME
'RUNID',id,...                  % SEMI-UNIQUE RUN ID FOR OUTPUT FILE NAMING
'DATEDIRPATH','.',...           % PATH TO DIRECTORY CONTAINING DATE DIRECTORIES
'DATEDIR','.',...               % DATE DIRECTORY CONTAINING INPUT DATA FILES
'OUTPUTDIR','.',...             % DIRECTORY TO PUT OUTPUT FILES
...
...% TRAVEL TIMES
'PHASE','P',...                 % SEISMIC PHASE NAME
'PULLARR',1,...                 % HOW TO GET TT (0=TTBOX,1=PULLARR,2=USE TTFIELD)
...
...% BAND LIMITS
'STOPSHORT',[],...              % STOP BAND SHORT PERIOD LIMIT
'PASSSHORT',10,...              % SHORTEST PERIOD OF CURRENT INTEREST
'PASSLONG',100,...              % LONGEST PERIOD OF CURRENT INTEREST
'STOPLONG',[],...               % STOP BAND LONG PERIOD LIMIT
...
...% DISTANCE/AZIMUTHAL LIMITS
'DISTCUT1',-1,...               % CLOSEST DISTANCE ALLOWED
'DISTCUT2',181,...              % FARTHEST DISTANCE ALLOWED
'AZICUT1',-181,...              % LOWEST AZIMUTH ALLOWED
'AZICUT2',361,...               % HIGHEST AZIMUTH ALLOWED
...
...% TEMP HEADER STORAGE
'TTFIELD','unused6',...         % WORKING LOCATION OF TT IN HEADER
'SNRFIELD','unused7',...        % WORKING LOCATION OF SNR IN HEADER
'AMPFIELD','unused8',...        % WORKING LOCATION OF AMPLITUDES IN HEADER
'POLFIELD','unused9',...        % WORKING LOCATION OF POLARITIES IN HEADER
...
...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
...%%%%                     START PREP SECTION                         %%%%
...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
...
...% PRE-FILTER WINDOW
'PREWINDOW',0,...               % PREFILTER WINDOWING
'PREWINREF1','b',...            % WINDOW START REFERENCE
'PREWINOFF1',0,...              % WINDOW START REFERENCE OFFSET
'PREWINREF2','e',...            % WINDOW END REFERENCE
'PREWINOFF2',0,...              % WINDOW END REFERENCE OFFSET
'PREWINFILL',1,...              % WINDOW FILL
'PREWINFILLER',0,...            % WINDOW FILLER
...
...% DOUBLE PRECISION
'DP',1,...                      % MAKE DATA DOUBLE PRECISION
...
...% PRE-FILTER CLEANING
'RDRIFT0',3,...                 % TREND REMOVAL (0=NONE,1=RMEAN,2=RSLOPE,3=RDRIFT)
'RDEAD0',1,...                  % REMOVE DEAD LOGICAL (0 = KEEP DEAD)
'ROFFSET0',0,...                % OFFSET REMOVER LOGICAL (0 = NO)  *** NOT IMPLEMENTED YET
'RGLITCHES0',0,...              % GLITCH REMOVER LOGICAL (0 = NO)  *** NOT IMPLEMENTED YET
...
...% PRE-FILTER TAPERING
'TAPER0',1,...                  % TAPERING LOGICAL (0 = NO TAPER)
'TAPERHW0',0.05,...             % TAPER HALFWIDTH (FRACTION OF WINDOW)
'TAPERTYPE0','hann',...         % TAPER TYPE (LOOK AT MATLAB FUNCTION WINDOW)
'TAPEROPT0',[]',...             % TAPER OPTION (LOOK AT MATLAB FUNCTION WINDOW)
...
...% PRE-FILTER RESAMPLE
'RESAMPLE0',1,...               % RESAMPLE LOGICAL (0 = LEAVE SAMPLE RATES ALONE)
'INTERP0',1,...                 % INTERPOLATE LOGICAL (0 = USE RESAMPLE)
'INTERPMETH0','spline',...      % INTERPOLATION METHOD
'RATE0',0.5,...                 % NEW SAMPLING FREQUENCY
...
...% FILTER
'FILTER',1,...                  % FILTER LOGICAL (0 = NO FILTER)
'TYPE','bandpass',...           % FILTER TYPE (low,high,notch,bandpass)
'STYLE','butter',...            % FILTER STYLE (butter,cheby1,cheby2,ellip)
'ORDER',0,...                   % FILTER ORDER (1 TO ~13) ([] or 0 = AUTO)
'NPASS',1,...                   % NUMBER OF PASSES (1 OR 2)
'LIMITS',[],...                 % FILTER LIMITS (IN Hz) (LEAVE BLANK FOR AUTO)
'RIPPLE',[],...                 % RIPPLE VALUES (IN DB)
...
...% POST-FILTER TAPERING 1
'TAPER1',1,...                  % TAPERING LOGICAL (0 = NO TAPER)
'TAPERHW1',0.05,...             % TAPER HALFWIDTH (FRACTION OF WINDOW)
'TAPERTYPE1','hann',...         % TAPER TYPE (LOOK AT MATLAB FUNCTION WINDOW)
'TAPEROPT1',[]',...             % TAPER OPTION (LOOK AT MATLAB FUNCTION WINDOW)
...
...% POST-FILTER CLEANING
'RDRIFT1',3,...                  % TREND REMOVAL (0=NONE,1=RMEAN,2=RSLOPE,3=RDRIFT)
'RDEAD1',1,...                   % REMOVE DEAD LOGICAL (0 = KEEP DEAD)
'ROFFSET1',0,...                 % OFFSET REMOVER LOGICAL (0 = NO)  *** NOT IMPLEMENTED YET
'RGLITCHES1',0,...               % GLITCH REMOVER LOGICAL (0 = NO)  *** NOT IMPLEMENTED YET
...
...% POST-FILTER TAPERING 2
'TAPER2',1,...                  % TAPERING LOGICAL (0 = NO TAPER)
'TAPERHW2',0.05,...             % TAPER HALFWIDTH (FRACTION OF WINDOW)
'TAPERTYPE2','hann',...         % TAPER TYPE (LOOK AT MATLAB FUNCTION WINDOW)
'TAPEROPT2',[]',...             % TAPER OPTION (LOOK AT MATLAB FUNCTION WINDOW)
...
...% ANALYSIS RESAMPLE
'RESAMPLE1',1,...               % RESAMPLE LOGICAL (0 = LEAVE SAMPLE RATES ALONE)
'INTERP1',1,...                 % INTERPOLATE LOGICAL (0 = USE RESAMPLE)
'INTERPMETH1','spline',...      % INTERPOLATION METHOD
'RATE1',5,...                   % NEW SAMPLING FREQUENCY
...
...% ANALYSIS GROUND UNITS
'GUNITS','velo',...             % GROUND UNITS TO USE IN ANALYSIS (disp,velo,accel)
...
...% ANALYZE ENVELOPES
'ENV',0,...                     % ENVELOPE LOGICAL (0 = NO ENVELOPE)
...
...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
...%%%%                      END PREP SECTION                          %%%%
...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
...
...% QUICK SNR
'QCKSNR',1,...                  % QUICK SNR LOGICAL (0 = SKIP)
'QCKSNRCUT',0,...               % SNR CUTOFF (0 = NO CUT)
'NOISWIN',[-100 -20],...        % NOISE WINDOW
'SIGWIN',[-20 50],...           % SIGNAL WINDOW (ALSO DEFAULT NON-INTERACTIVE WINDOW)
...
...% WINDOWING
'WINDOW',1,...                  % WINDOWING LOGICAL (0 = NO WINDOW)
'USERWIN',1,...                 % INTERACTIVE WINDOWING (0 = NO)
'INTWIN',[-200 200],...         % INITIAL WINDOW IN INTERACTIVE MODE
'FILL',1,...                    % ADD FILLER IF WINDOW GOES OUT OF BOUNDS
'FILLER',0,...                  % FILLER TO ADD
...
...% TAPERING
'TAPER',1,...                   % TAPERING LOGICAL (0 = NO TAPER)
'USERTAPER',1,...               % INTERACTIVE TAPERING (0 = NO)
'TAPERHW',0.05,...              % TAPER HALFWIDTH (FRACTION OF WINDOW)
'TAPERTYPE','hann',...          % TAPER TYPE (LOOK AT MATLAB'S WINDOW FUNCTION)
'TAPEROPT',[]',...              % TAPER OPTION (LOOK AT MATLAB'S WINDOW FUNCTION)
...
...% CORRELATE
'UCORR',1,...                   % INTERACTIVE CORRELATION (0 = NO)
'NPEAKS',3,...                  % NUMBER OF CORRELEGRAM PEAKS TO PICK
'SPACING',[],...                % NUMBER OF SAMPLES BETWEEN PEAKS (LEAVE BLANK FOR AUTO)
'ABSXC',true,...                % PICK ABSOLUTE MAXIMAS (1 = YES)
'NORMXC',true,...               % NORMALIZE CORRELOGRAMS BY AUTOCORRELATIONS (1 = YES)
'POW2PAD',1,...                 % ZERO PADDING CONTROL (POWER OF 2 ADJUSTMENT)
...
...% CLUSTERING (LOCALIZED NOT IMPLEMENTED)
'CLUSTER',1,...                 % CLUSTERING (0=NONE,1=LOCALIZED,2=GLOBAL,3=BOTH)
'USERCLUSTER',1,...             % INTERACTIVE GLOBAL CLUSTER ANALYSIS (0 = NO)
'CMETHOD','average',...         % GLOBAL CLUSTERING METHOD (BETTER LEAVE THIS ALONE)
'LOCALDIST',0.1,...             % DISSIMILARITY LIMIT FOR LOCALIZED CLUSTERING (0 TO 1)
'GLOBALDIST',0.35,...           % DISSIMILARITY LIMIT FOR GLOBAL CLUSTERING (0 TO 1)
'LOCALMINSIG',2,...             % MINIMUM POPULATION A SIGNAL MUST BE LOCALLY SIMILAR TO (2+)
'GLOBALMINSIG',0.1,...          % GLOBAL CLUSTER POPULATION LIMIT (INTEGER OR FRACTION
...                             % OF TOTAL POPULATION (0 TO 1 OR 3+)
...
...% SALVAGE (NOT IMPLEMENTED)
'RAISE',1,...                   % SALVAGE LOGICAL (0 = LEAVE REJECTED ALONE)
'NOSMALLGROUPS',0,...           % SALVAGE SINGLES & DOUBLES ONLY OR CHECK SMALL GROUPS TOO
'USERCHECK',1,...               % LET USER LOOK OVER FINAL CLUSTERS
'POLISH',1,...                  % APPLY CORRECTIONS TO WRITTEN DATA
...
...% SOLVER
'WEIGHTPOW',2,...               % WEIGHTING POWER OF CORRELATIONS (0 = NO WEIGHTING)
'MRI',20,...                    % MAX REFINEMENT ITERATIONS
'NSTDDEV',2,...                 % NUMBER OF STD DEVIATIONS FOR ERROR CALCULATIONS
...
...% STACK PICKS (NOT IMPLEMENTED)
'USERSTACK',1,...               % INTERACTIVE PICKING OF STACKED TRACES
'TIESTACKS',1,...               % TIE PICKED CLUSTERS TOGETHER
...
...% OUTLIER ANALYSIS
'USERCUT',0,...                 % INTERACTIVE OUTLIER ANALYSIS LOGICAL
'NSTDCUT',2,...                 % NUM STD DEVS FOR DEFAULT CUT
'TIMECUT',0,...                 % RUN CUT BASED ON REL ARR VS DISTANCE
'AMPCUT',0,...                  % RUN CUT BASED ON AMPLITUDE VS DISTANCE
'STDCUT',0 ...                  % RUN CUT BASED ON STD DEV OF REL ARR
...
);

end
