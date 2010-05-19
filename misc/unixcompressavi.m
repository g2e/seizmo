function [varargout]=unixcompressavi(filein,fileout,codec)
%UNIXCOMPRESSAVI    Compress an AVI file in Unix with "MEncoder"
%
%    Usage:    unixcompressavi
%              unixcompressavi(filename)
%              unixcompressavi(filein,fileout)
%              unixcompressavi(filein,[],codec)
%              unixcompressavi(filein,fileout,codec)
%              output=unixcompressavi(...)
%
%    Description: UNIXCOMPRESSAVI shows the available MEncoder video
%     "codecs" for compressing an AVI file.  Codecs that I have used
%     successfully are 'lavc' (the default) & 'x264'.
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
%     'lavc' and there is little reason to adjust this unless you are
%     making some particularly lengthy movies (in that case 'x264' may be
%     a better choice).
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
%     I use this for compressing movies made with my fk routines:
%      s3d=fkxcvolume(xcdata,50,201,[1/100 1/3]);
%      mov=fkfreqslide(s3d,0);
%      movie2avi(mov,'xcfk_vs_freq.avi');
%      unixcompressavi('xcfk_vs_freq.avi');
%      unix('vlc xcfk_vs_freq.avi');
%
%    See also: MOVIE2AVI, AVIINFO, MOVIE, GETFRAME, UNIX

%     Version History:
%        May  15, 2010 - initial version
%        May  18, 2010 - slight doc touch
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  18, 2010 at 23:15 GMT

% todo

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
error(nargchk(1,3,nargin));

% default codec
if(nargin==1 || isempty(fileout)); fileout=filein; end
if(nargin==2 || isempty(codec)); codec='lavc'; end

% check inputs
if(ischar(filein)); filein=cellstr(filein); end
if(ischar(fileout)); fileout=cellstr(fileout); end
if(ischar(codec)); codec=cellstr(codec); end
szi=numel(filein); szo=numel(fileout); szc=numel(codec);
if(~iscellstr(filein) || ~iscellstr(fileout) || ~iscellstr(codec))
    error('seizmo:unixcompressavi:badInput',...
        'Inputs must all be strings!');
elseif((szi~=1 || szc~=1) && szo==1)
    error('seizmo:unixcompressavi:writeToSingleFile',...
        'Single output file for multiple input files/codecs not allowed!');
elseif(~isequalsizeorscalar(filein,fileout,codec))
    error('seizmo:unixcompressavi:badInput',...
        'FILEIN/FILEOUT/CODEC must be strings or equal-sized cellstr!');
end

% expand scalars
n=max([szi szo szc]);
if(isscalar(filein)); filein(1:n,1)=filein; end
if(isscalar(fileout)); fileout(1:n,1)=fileout; end
if(isscalar(filecodec)); codec(1:n,1)=codec; end

% loop over filein checking each exists
for i=1:n
    if(~exist(filein{i},'file'))
        error('seizmo:unixcompressavi:fileNotFound',...
            'Could not locate file: %s',filein{i});
    end
end

% encode files
tmp=tempname;
for i=1:n
    % write to a temporary name so we don't clobber if filein==fileout
    [s,w]=unix(['mencoder ' filein{i} ' -o ' tmp ' -ovc ' codec{i}]);
    if(s); error('seizmo:unixcompressavi:mencoderError',w); return; end
    [s,w1]=unix(['mv ' tmp ' ' fileout{i}]);
    if(s); error('seizmo:unixcompressavi:moveFileError',w1); return; end
end

% output
if(nargout); varargout{1}=w; end

end
