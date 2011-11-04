function []=daydirs_workflow()
%DAYDIRS_WORKFLOW    Outlines SEIZMO Noise Correlation Processing Workflow
%
% Making NCFs & EGFs:
%    0. Type "help sz_toc_noise" at the Matlab prompt to bring up the list
%       of noise routines and their short descriptions. You can click on
%       the functions in this list to bring up more detailed help or just
%       type "help daydirs_XXXXXXXX" where XXXXXXXX is the function name.
%    1. Depends on your dataset.
%       a. If you have a seed volume:
%          1. Extract the seismograms into a directory as SAC files -- see
%             rdseed for instructions on how to do that.
%          2. Fix the headers (in Matlab):
%              w(fix_rdseed_v48(r('mydir')))
%             This isn't necessary but will silence some warnings.
%          3. If you are working with horizontals, make sure they are all
%             in the North or East directions. You can use (in Matlab):
%              data=r('mydir');
%              w(rotate(data,'to',0),'dir','rotated');
%              w(data(vertcmp(data)),'dir','rotated');
%          4. Use DAYDIRS_MAKE on that directory (in Matlab):
%              daydirs_make('rotated','STEP_1-DAYS')
%       b. Otherwise things are more difficult. The routines expect you to
%          have day-long SAC files organized in a year/day/file filesystem.
%    2. Create 25 hour seismograms by running DAYDIRS_MERGECUT_25HRS on the
%       dataset (in Matlab):
%        daydirs_mergecut_25hrs('STEP_1-DAYS','STEP_2-25HRS')
%       This is how I did my processing but many people do this differently
%       now (like 15 minute, 1 hour or 1 day files with varying overlap). I
%       am working on making these choices possible -- stay tuned.
%    3. Resample to 1sps or some other samplerate if you choose. This is
%       optional if all your data have the same samplerate. Let me just say
%       that 1sps is great for diskspace and typical noise studies.
%        daydirs_resample('STEP_2-25HRS','STEP_3-1SPS')
%    4. Remove the instrument response. If your response is already removed
%       then you can skip this step. You might also skip this step if all
%       your instruments have the same response.
%        daydirs_rinst('STEP_3-1SPS','STEP_4-RINST')
%    5. Whiten the seismograms. Basically we are reducing the influence of
%       earthquakes and glitches. This does both time-domain and
%       frequency-domain whitening. See the code for details.
%        daydirs_normalize('STEP_4-RINST','STEP_5-NORM')
%    6. Correlate the noise:
%        daydirs_correlate('STEP_5-NORM','STEP_6-XCORR')
%    7. Rotate horizontal correlations. Remember they need to already be in
%       the North/East orientations (I could fix this but I'm lazy).
%        daydirs_rotcorr('STEP_6-XCORR','STEP_7-ROTXC')
%    8. Make stacks. Again implementation varies by seismologist so this
%       may not do what you want (or it may do way more than you want). See
%       "help daydirs_stackcorr" for more details. This function is also
%       more complicated to use (right now) because it needs all the
%       station names as an input and it only does one pairing type at a
%       time. Note there are no verticals in the STEP_7-ROTXC directory.
%        stns={'STN1' 'BLAH' 'CMB' ... 'STNXX'}; % NEED TO PUT YOUR STNS!!!
%        daydirs_stackcorr('STEP_6-CORR','STEP_8-STACK',stns,'Z','Z')
%       If you are also doing horizontals:
%        daydirs_stackcorr('STEP_7-ROTXC','STEP_8-STACK',stns,'R','R')
%        daydirs_stackcorr('STEP_7-ROTXC','STEP_8-STACK',stns,'R','T')
%        daydirs_stackcorr('STEP_7-ROTXC','STEP_8-STACK',stns,'T','R')
%        daydirs_stackcorr('STEP_7-ROTXC','STEP_8-STACK',stns,'T','T')
%
% Making those pretty noise correlation vs distance plots or movies:
%    0. Do the above.
%    1. Read in the particular NCFs or EGFs you are interested in. Example:
%       data=r('STEP_8-STACK/ZZ/GREEN_FULLSTACK');
%    2. Make a plot or movie:
%       a. aziwindow(data) or aziwindow(data,[0 30])
%       b. mov=azisweep(data);
%          movie2avi(mov,'mymovie.avi');
%          % If you are using linux and have mencoder installed:
%          % unixcompressavi('mymovie.avi')
%
% There are also functions to make fk datasets from the NCFs and to make
% plots for those. That should be pretty easy to figure out so I leave
% that for you to explore for now -- see "help sz_toc_noise" for a list.
%
% <a href="matlab:help sz_toc_noise">Seismic Noise Analysis Functions</a>
% <a href="matlab:help seizmo">SEIZMO - Passive Seismology Toolbox</a>
end
