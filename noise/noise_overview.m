function []=noise_overview()
%NOISE_OVERVIEW    Outlines Noise Correlation Analysis Workflow
%
% Making NCFs:
%    0. Type "help sz_toc_noise" at the Matlab prompt to bring up the list
%       of noise routines and their short descriptions. You can click on
%       the functions in that list to bring up more detailed help or just
%       type "help XXXXXXXX" where XXXXXXXX is the function name.
%    1. Depends on your dataset.
%       a. If you have a seed volume:
%          1. Extract the seismograms into a directory as SAC files using
%             RDSEED (http://www.iris.edu/forms/rdseed_request.htm).
%             Help using RDSEED: http://www.iris.edu/manuals/rdseed.htm
%          2. Fix the headers (in Matlab):
%              w(fix_rdseed_v48(r('mydir')))
%             This isn't necessary but will silence some warnings.
%       b. Otherwise you just need to get all the records as SAC files
%          under the same directory and make sure following header fields
%          are set correctly:
%           KNETWK, KSTNM, KHOLE, KCMPNM, NZYEAR, NZJDAY, NZHOUR, NZMIN,
%           NZSEC, NZMSEC, B, E, NPTS, DELTA, CMPINC, CMPAZ, STLA, STLO,
%           STDP, STEL, LEVEN (MUST BE EVENLY SPACED!!!), IFTYPE, IZTYPE
%    2. Create 3hr timesections to be processed separately (in Matlab):
%        noise_setup('mydir','setup');
%       Note that if your data files are under subdirectories you will need
%       to use something like the following instead to force a recursive
%       search for files:
%        noise_setup('mydir','setup','f','**/');
%    3. Process the timesections (in Matlab):
%        noise_process('setup','ncf');
%    4. Stack the noise correlations (in Matlab):
%        noise_stack('ncf','stacks','zz');
%       And for horizontals:
%        noise_stack('ncf','stacks',{'rr' 'rt' 'tr' 'tt'});
%
%    Steps 2-4 can be modified significantly by passing options to the
%    corresponding function.  See their respective help documentation for
%    details.
%
% Making a pretty ncfs vs distance plot or movie:
%    0. Do the above.
%    1. Read in the particular NCFs you are interested in:
%        data=r('stacks/zz_all_full/*/');
%    2. Make a plot or movie:
%       a. aziwindow(data) or aziwindow(data,[0 30])
%       b. mov=azisweep(data);
%          movie2avi(mov,'mymovie.avi');
%          % If you are using linux and have mencoder installed:
%          unixcompressavi('mymovie.avi')
%
% There are also functions to make fk datasets from the NCFs and to make
% plots for those. That should be pretty easy to figure out so I leave
% that for you to explore for now -- see "help sz_toc_noise" for a list.
%
% <a href="matlab:help sz_toc_noise">Seismic Noise Analysis Functions</a>
% <a href="matlab:help seizmo">SEIZMO - Passive Seismology Toolbox</a>
end
