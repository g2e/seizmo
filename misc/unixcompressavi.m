function [varargout]=unixcompressavi(filein,fileout,codec,options)
%UNIXCOMPRESSAVI    Compress an AVI file in Unix with "MEncoder"
%
%    Usage:    unixcompressavi
%              unixcompressavi(filename)
%              unixcompressavi(filein,fileout)
%              unixcompressavi(filein,[],codec)
%              unixcompressavi(filein,fileout,codec)
%              unixcompressavi(filein,fileout,codec,options)
%              output=unixcompressavi(...)
%
%    Description:
%     UNIXCOMPRESSAVI shows the available MEncoder video "codecs" for
%     compressing an AVI file.  Codecs that I have used successfully are
%     'lavc' (the default) & 'x264'.
%
%     UNIXCOMPRESSAVI(FILENAME) compresses an AVI file FILENAME using
%     MEncoder & the AVC codec.  This results in massive reductions in file
%     size (99% or better is typical when going from an uncompressed AVI).
%     FILENAME may be a single filename string or a cell array of such
%     (to encode multiple files).  DOES NOT SUPPORT WILDCARDS!!!
%
%     UNIXCOMPRESSAVI(FILEIN,FILEOUT) takes the raw AVI file FILEIN,
%     compresses it using the default AVC codec and writes it to a new file
%     FILEOUT.  FILEIN & FILEOUT may be single filenames or equal sized
%     cell arrays of filenames.
%
%     UNIXCOMPRESSAVI(FILEIN,[],CODEC) or
%     UNIXCOMPRESSAVI(FILEIN,FILEOUT,CODEC) allows specifying the codec
%     used (see first Usage form for available codecs).  The default is
%     'lavc' and is known to make avi files that are compatible with Linux
%     & Windows.  You will want to do a 'man mencoder' at a bash terminal
%     if you think you want to adjust this.  Altering CODEC will require
%     adjusting OPTIONS (next Usage form) as they are set up for the
%     default codec.
%
%     UNIXCOMPRESSAVI(FILEIN,FILEOUT,CODEC,OPTIONS) allows altering the
%     encoding options being passed to the encoder.  The default OPTIONS
%     string is the following:
%      vbitrate=2560000:vcodec=msmpeg4v2:vqblur=1.0: (continues)
%      mbd=2:cmp=2:subcmp=2:dia=2:last_pred=3:mv0
%     This does a good job at compressing the file while not pixelating the
%     video nearly as much as Mencoder's defaults would.  This will be need
%     to be altered if CODEC is not 'lavc'!
%
%     OUTPUT=UNIXCOMPRESSAVI(...) exports the output from the last encoding
%     run to the string OUTPUT.  Note this will only contain info on the
%     last file to be encoded and no other.
%
%    Notes:
%     - UNIXCOMPRESSAVI will throw an error if MEncoder indicates there is
%       a problem.
%     - Temporary filenames are created using TEMPNAME.  This may be
%       useful in tracking down a temporary file if something goes wrong.
%
%    Examples:
%     % I use this for compressing movies made with my fk routines:
%     s3d=fkxcvolume(xcdata,50,201,[1/100 1/3]);
%     mov=fkfreqslide(s3d,0);
%     movie2avi(mov,'xcfk_vs_freq.avi');
%     unixcompressavi('xcfk_vs_freq.avi');
%     unix('mplayer xcfk_vs_freq.avi');
%
%     % To encode a list of avi files:
%     files=xdir('*.avi');
%     unixcompressavi({files.name});
%
%    See also: MOVIE2AVI, AVIINFO, MOVIE, GETFRAME, UNIX

%     Version History:
%        May  15, 2010 - initial version
%        May  18, 2010 - slight doc touch
%        May  24, 2010 - fixed single arg bug & one other
%        May  25, 2010 - default encode string should improve windows
%                        compatibility, added .avi to tmp name, doc update,
%                        2pass encoding, options argument
%        Apr.  3, 2012 - minor doc update
%        Apr. 26, 2012 - add LD_LIBRARY_PATH="" to unix calls to avoid
%                        linker issues (glibcxx_3.4.11 no found...)
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 17:45 GMT

% todo
% - options for other codecs?

% require unix
if(~isunix)
    error('seizmo:unixcompressavi:notUnix',...
        'UNIXCOMPAVI can only be called on a Unix OS!');
end

% require mencoder
[s,w]=unix('which mencoder');
if(s)
    error('seizmo:unixcompressavi:noMencoder',...
        'UNIXCOMPAVI requires mencoder (aka mplayer) to be installed!');
end

% display available codecs if no input
if(~nargin); unix('mencoder -ovc help'); return; end

% check nargin
error(nargchk(1,4,nargin));

% directory separator
fs=filesep;

% default codec
if(nargin==1 || isempty(fileout)); fileout=filein; end
if(nargin<3 || isempty(codec)); codec='lavc'; end
if(nargin<4 || isempty(options))
    options=['vbitrate=2560000:vcodec=msmpeg4v2:vqblur=1.0:' ...
        'mbd=2:cmp=2:subcmp=2:dia=2:last_pred=3:mv0'];
end

% check inputs
if(ischar(filein)); filein=cellstr(filein); end
if(ischar(fileout)); fileout=cellstr(fileout); end
if(ischar(codec)); codec=cellstr(codec); end
if(ischar(options)); options=cellstr(options); end
szin=numel(filein); szout=numel(fileout);
szenc=numel(codec); szopt=numel(options);
if(~iscellstr(filein) || ~iscellstr(fileout) ...
        || ~iscellstr(codec) || ~iscellstr(options))
    error('seizmo:unixcompressavi:badInput',...
        'Inputs must all be strings!');
elseif((szin~=1 || szenc~=1 || szopt~=1) && szout==1)
    error('seizmo:unixcompressavi:writeToSingleFile',...
        'Single output file for multiple input files/codecs not allowed!');
elseif(~isequalsizeorscalar(filein,fileout,codec,options))
    error('seizmo:unixcompressavi:badInput',...
        ['FILEIN/FILEOUT/CODEC/OPTIONS must be simple strings ' ...
        'or equal-sized cellstr arrays!']);
end

% expand scalars
n=max([szin szout szenc szopt]);
if(isscalar(filein)); filein(1:n,1)=filein; end
if(isscalar(fileout)); fileout(1:n,1)=fileout; end
if(isscalar(codec)); codec(1:n,1)=codec; end
if(isscalar(options)); options(1:n,1)=options; end

% loop over filein checking each exists
for i=1:n
    if(~isabspath(filein{i})); filein{i}=[pwd fs filein{i}]; end
    if(~exist(filein{i},'file'))
        error('seizmo:unixcompressavi:fileNotFound',...
            'Could not locate file: %s',filein{i});
    end
end

% codec to opt
codecs={'lavc' 'faac' 'xvid' 'x264' 'nuv'};
opt={'-lavcopts' '-faacopts' '-xvidencopts' '-x264encopts' '-nuvopts'};

% encode files
% - note we are using a 2-pass encoding scheme for better quality
% - writes to a temporary filename to avoid clobber (and total file loss)
tmp=[tempname '.avi'];
for i=1:n
    % what is our option string?
    copt=opt{strcmpi(codecs,codec{i})};
    
    % 2pass encoding
    disp(['LD_LIBRARY_PATH="" && mencoder ' filein{i} ' -nosound -o /dev/null' ...
        ' -ovc ' codec{i} ' ' copt ' ' options{i} ':vpass=1']);
    [s,w]=unix(['LD_LIBRARY_PATH="" && mencoder ' filein{i} ' -nosound -o /dev/null' ...
        ' -ovc ' codec{i} ' ' copt ' ' options{i} ':vpass=1']);
    if(s); error('seizmo:unixcompressavi:mencoderError',w); return; end
    disp(['LD_LIBRARY_PATH="" && mencoder ' filein{i} ' -nosound -o ' tmp ...
        ' -ovc ' codec{i} ' ' copt ' ' options{i} ':vpass=2']);
    [s,w]=unix(['LD_LIBRARY_PATH="" && mencoder ' filein{i} ' -nosound -o ' tmp ...
        ' -ovc ' codec{i} ' ' copt ' ' options{i} ':vpass=2']);
    if(s); error('seizmo:unixcompressavi:mencoderError',w); return; end
    [s,w1]=unix(['mv ' tmp ' ' fileout{i}]);
    if(s); error('seizmo:unixcompressavi:moveFileError',w1); return; end
    [s,w1]=unix('rm divx2pass.log');
    if(s); error('seizmo:unixcompressavi:removeFileError',w1); return; end
end

% output
if(nargout); varargout{1}=w; end

end
