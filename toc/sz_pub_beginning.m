%% SEIZMO - A Matlab(R) & Octave Toolbox for Passive Seismology
% 
%% What is SEIZMO?
% <http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/ SEIZMO> is a
% <http://www.mathworks.com Matlab(R)> and
% <http://www.gnu.org/software/octave/ Octave> based toolbox encompassing a
% collection of <matlab:doc('alphabetical_list') over 700 functions> that
% provide a framework for seismic data preparation, quality control, and
% analysis akin to that of
% <http://www.iris.edu/software/sac/manual.html Seismic Analysis Code> with
% a focus on passive source seismology at regional and global scales (where
% earthquakes and microseisms are the sources).

%% Why use SEIZMO?
%
% <http://www.iris.edu/software/sac/manual.html SAC> is syntactically
% simpler than SEIZMO and starts significantly faster than Matlab or Octave
% -- so why use SEIZMO+Matlab/Octave?
%
% # the ease of manipulating variables
% # the ease of minipulating plots
% # the simplicity of mathematical syntax
% # the huge library of mathematical functions
% # the extensive & browsible documentation of Matlab
%
% These are great reasons but SAC already has a method for interfacing with
% Matlab.  So where does SEIZMO come in?
%
% # SAC doesn't work with Octave (and Matlab is expensive).
% # SAC's Matlab functions are few.
% # SEIZMO can do most SAC operations natively in Matlab/Octave.
% # SEIZMO simplifies the more complicated processing steps (
%     <matlab:helpwin('rotate') rotation>, <matlab:helpwin('meld') merging>,
%     <matlab:helpwin('removesacpz') response removal>, etc) so you can focus
%     on new science.
% # SEIZMO has <matlab:doc('sz_toc_models') 1D & 3D Models>,
%     <matlab:doc('sz_toc_cmt') CMTs>,
%     <matlab:doc('sz_toc_mattaup') raypaths>,
%     <matlab:doc('sz_toc_cmb') waveform cluster analysis>,
%     <matlab:doc('sz_toc_noise') noise analysis>,
%     <matlab:doc('sz_toc_fk') beamforming analysis> and
%     <matlab:doc('categorical_list') so much more> built in.
%
% Convinced that SEIZMO+Matlab/Octave is a viable alternative to SAC,
% SAC+Matlab, etc?  Read on!

%% How to read a file into SEIZMO
% Currently SEIZMO only supports one type of seismic data format: SAC
% binary.  It is a letdown but I promise the rest gets better!  Reading in
% the SAC files is done with <matlab:helpwin('readseizmo') |readseizmo|> or
% the <matlab:doc('sz_toc_shortnames') shortform>: |r|.  For
% example, to read in all the SAC files in a directory |data| into the
% Matlab variable called |dataset|:
dataset=readseizmo('data');

%%
% Or using the shortcut form (replacing |readseizmo| with just |r|):
dataset=r('data');

%%
% The data can then be plotted using one of SEIZMO's
% <matlab:doc('sz_toc_plotting') plotting commands>:
recordsection(dataset)

%% How to save a file from SEIZMO
% There are 2 options when saving SEIZMO datasets: as a |.mat| file or as
% SAC files.  The MAT file option allows you to save the entire dataset as
% well as any other variables you request into a single file.  The downside
% is that this file is not readable by other seismology programs such as
% SAC or
% <http://www.passcal.nmt.edu/content/pql-ii-program-viewing-data PQLII>.
save myfile.mat dataset

%%
% To save the records in the dataset individually as SAC files use the
% command <matlab:doc('writeseizmo') |writeseizmo|> (here we will use the
% shortcut |w|).  We also give a parameter & value pair to change the path
% of the output files to |data-new|.
w(dataset,'path','data-new');

%%
% Listing the 2 directories shows all is well:
ls -n data/ data-new/

%% How to create a SEIZMO dataset from a matrix
% Creating a SEIZMO dataset from a matrix is done with the |bseizmo|
% command.  For example, we can create a vector of points with random
% values and pass those to <matlab:doc('bseizmo') |bseizmo|>.  The output
% is a SEIZMO dataset (we will explore the SEIZMO data format in detail in
% the following sections).  We then can plot the matrix and SEIZMO dataset
% to verify.
x=rand(1000,1);
data=bseizmo(x);

% plotting...
figure;
subplot(2,1,1);
plot(x);
title('normal Matlab matrix');
subplot(2,1,2);
plot1(data,'axis',gca);
title('now in SEIZMO')

%% The SEIZMO struct
% In this section we explore the main parts of a SEIZMO dataset.  First,
% the dataset is actually stored as a struct in Matlab.  A "struct" is
% simply a structured set of variables.  The SEIZMO struct contains 10
% variables or "fields" that organize all the info in a record.  Those 10
% fields are:
%
% * |path|      - directory of file
% * |name|      - file name
% * |filetype|  - type of file
% * |version|   - version of filetype
% * |byteorder| - byte-order of file (ieee-le or ieee-be)
% * |head|      - header data
% * |hasdata|   - logical indicating if data is read in
% * |ind|       - independent component data (for uneven)
% * |dep|       - dependent component data
% * |misc|      - place for miscellaneous record info

%%
% To list the fields and their values for a record, enter the dataset name
% and record index at the commandwindow.  For instance, record 3 of the
% dataset from our reading SAC files example (we saved the dataset to a
% variable called |dataset|) can be explored by entering:
dataset(3)

%%
% We quickly assess that this data is from the |data| folder and more
% specifically is the file |SAC.XB.CM03.02.BHZ.00| in that folder.  The
% filetype is *SAC binary version 6* as expected (all records will probably
% have those entries until I add support for other filetypes).  The |head|
% field contains all the metadata in a |302x1| double-precision array.  The
% |hasdata| field is |1| or logically |TRUE| and denotes that we have read
% in the data, which is stored in the |dep| field (ie dependent data
% component) as a |4500x1| single-precision array.  The |ind| field (ie
% independent data component) is empty, indicating that the data is evenly
% sampled in time because we do not need to store the timing of every
% sample for such a record.  The |misc| field contains nothing at this
% point but may be populated later by other SEIZMO functions for keeping
% track of related information (like the instrument response).

%% Copying a SEIZMO struct
% To copy a dataset, assign it to a new variable:
dataset2=dataset

%%
% If you are only interested in a single record, you can save that record
% to a new dataset by assigning that record to a new variable:
new=dataset(3)

%%
% You may then assign it back by switching the two:
dataset(3)=new;

%% Altering SEIZMO struct fields
% Now that you have a taste of the SEIZMO struct, the next step is to learn
% how to apply that knowledge by changing the struct fields.  While I do
% not recommend altering the |filetype| & |version| fields as there is no
% reason, altering other fields is convenient.  For example, changing field
% values allows you to alter the filename & path when the record is
% written out as a SAC file.  Altering record 3 of the dataset:
dataset(3).name='mynewname.sac';
dataset(3).path='data-new';

%%
% Display the record fields to check:
dataset(3)

%%
% And writing out to a SAC file:
w(dataset(3));
ls -n data-new/*.sac

%%
% In summary, to adjust a field for a record the format is as follows:
%
% |datasetname(recordindex).field=value|
%
% If there is only 1 record in a dataset then the |(recordindex)| may be
% omitted.

%%
% Another easy struct adjustment example is to change the data of a
% record.  Say you wanted to add some white noise to the record.  First,
% make a copy so you can compare the noisy signal to the original.  Then
% add some noise and plot the two records in an overlay using
% <matlab:doc('plot2') |plot2|>:
noisy=new;
noisy.dep=noisy.dep+(rand(4500,1)-0.5);
plot2(noisy,new)

%%
% In this case, the noise is too weak to affect the character of the
% signal.  Enhance the noise |2000x| for a pronounced effect:
noisy=new;
noisy.dep=noisy.dep+2e3*(rand(4500,1)-0.5);
plot2(noisy,new)

%% Extracting the data
% The quickest way to access the data in a SEIZMO record is to use a
% command with the following form:
%
% |mymatrix=datasetname(recordindex).dep;|
%
% For example, to extract the first 10 values from the noisy record:
a=noisy.dep(1:10)

%% Viewing header info
% There are 3 different header (metadata) viewers included in SEIZMO:
%
% * <matlab:doc('listheader') |listheader|>    - List SEIZMO data headers
% * <matlab:doc('compareheader') |compareheader|> - List SEIZMO headers in field x recond table form
% * <matlab:doc('queryheader') |queryheader|>   - List SEIZMO headers in record x field table form
%
% The difference between |compareheader| and |queryheader| is a
% transposition of the table.  I personally prefer |queryheader| aka |qh|.
%
% To list some header fields of the first 3 records in |dataset|:
lh(dataset(1:3),'delta','b','e','stla','stlo')

%%
% Compare that to |queryheader| output:
qh(dataset(1:3),'delta','b','e','stla','stlo')

%%
% Wildcards are also allowed:
qh(dataset(1:3),'l*')

%% Extracting header info
% <matlab:doc('getheader') |getheader|> (aka |gh|) allows for exporting
% header values to matlab variables.  For instance if you wanted to extract
% the beginning time of each record, you would ask for the |b| header
% field:
values=gh(dataset,'b')

%%
% Say you wanted to output a string header field such as |kstnm|.
% |getheader| returns a cell-string array in this case, which allows for
% simpler access to each record's string.  To access the 3rd record's
% |kstnm| value, index just like you would into the dataset:
values=gh(dataset,'kstnm')
values(3)

%%
% Converting a cell-string array to a character array is easily done using
% the |char| command:
values=char(values)

%%
% You can also return multiple header fields in the same call:
[delta,e]=gh(dataset,'delta','e')

%%
% Enumerated header fields are a little more complex. An integer is stored
% in the header location for an enum field.  This is what |gh| will return
% by default.  This integer corresponds to a specific string in a lookup
% table that SEIZMO keeps internally.  SEIZMO's strings match those in SAC
% and include a few extensions.  To return the id or description strings
% for an enum field, you have to add a modifier to the field string.
% Please note that these also return cell-string arrays.
gh(dataset(3),'idep')
gh(dataset(3),'idep id')
gh(dataset(3),'idep desc')

%%
% The main thing to remember about logical fields is that the SAC format
% allows for these fields to be undefined (ie set as -12345).  Thus a
% non-zero value returned by |gh| does not necessarily indicate the logical
% is |TRUE|.  It is typically safe to make that assumption though.

%% Altering header info
% Changing a header field to a new value is facilitated by
% <matlab:doc('changeheader') |changeheader|>
% aka |ch|.  An example of changing the |kt0| & |t0| header fields of all
% records in |dataset| to the same value:
dataset=ch(dataset,'kt0','nothing','t0',8000);

%%
% Note that the |ch| output was assigned back to the input dataset.  This
% isn't necessary and can be assigned to a new variable instead.
% To view the 'nothing' markers in an ammended <matlab:doc('plot1') |plot1|> call:
plot1(dataset(1:4),'showmarkers',true);

%%
% You may use an array to give each record a diffent value.  To add a
% little randomness to the marker positions:
dataset=ch(dataset,'t0',1e4+1e4*(rand(29,1)-0.5));
plot1(dataset(1:4),'showmarkers',true);

%% Beyond the basics
% Next: <matlab:showdemo('sz_pub_intermediate') Processing in SEIZMO>
