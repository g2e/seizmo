function []=writesacpz(file,z,p,k,o)
%WRITESACPZ    Writes out a SAC PoleZero file
%
%    Usage:    writesacpz(file,z,p,k)
%              writesacpz(file,z,p,k,overwrite)
%
%    Description:
%     WRITESACPZ(FILE,Z,P,K) writes the SAC PoleZero file FILE using the
%     zeros in Z, the poles in P, and the constant in K.  See the Notes
%     section for file format details.  Z and P must be numeric vectors
%     composed of reals and/or complex conjugate pairs.  K must be a real
%     scalar.
%
%     WRITESACPZ(FILE,Z,P,K,OVERWRITE) quietly overwrites pre-existing SAC
%     PoleZero files without confirmation when OVERWRITE is set to TRUE.
%     By default OVERWRITE is FALSE.
%
%    Notes:
%     - The format of a SAC PoleZero file is free format and is keyword
%       driven.  The keywords are 'ZEROS' 'POLES' and 'CONSTANT'.  Specify
%       the number of zeros by using the keyword 'ZEROS' followed by an
%       integer.  Subsequent lines are taken as the locations of zeros
%       until the next keyword is given.  The locations should be two
%       numbers specifying the value of the real and imaginary components.
%       Zeros located at the origin do not need to be given (by default
%       they are all assumed to be at the origin).  Poles may be specified
%       in the same manner but with the keyword 'POLES'.  Specify a scaling
%       constant by giving the keyword 'CONSTANT' followed by the number.
%       By default the scaling constant is 1.  An example:
%           ZEROS 3
%           POLES 5
%           -0.0370  0.0370
%           -0.0370  -0.0370
%           -251.3000  0.0000
%           -131.0000  467.3000
%           -131.0000  -467.3000
%           CONSTANT 5.588419e+16
%       Note that all the 3 zeros are at the origin, while all 5 poles are
%       not.  The first number in the lines listing the pole locations give
%       the real component and the second number gives the imaginary
%       component (thus there are 2 complex conjugate pairs).  In this case
%       the last line gives the multiplicative factor.  The order of these
%       three sections does not matter.
%
%    Examples:
%     % Read in a SAC PoleZero file, alter the constant, and write out:
%     [z,p,k]=readsacpz('SAC_PZs_XB_CM32_BHZ_02');
%     k=k*correction_factor;
%     writesacpz('SAC_PZs_XB_CM32_BHZ_02',z,p,k);
%
%    See also: READSACPZ, GETSACPZ, APPLYSACPZ, REMOVESACPZ, MAKESACPZDB,
%              PARSE_SACPZ_FILENAME, GENSACPZNAME, FIX_OLD_SACPZ,
%              WRITESACPZ_RDSEED, READSACPZ_RDSEED, ISSACPZ_RDSEED

%     Version History:
%        Apr.  7, 2009 - initial version
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Sep.  8, 2009 - fix some error ids
%        Sep. 20, 2009 - minor doc update
%        Sep. 22, 2009 - dropped looped writing of zeros/poles,
%                        confirmation for overwrite with skip option
%        Feb.  5, 2010 - graphical file creation menu
%        Feb. 11, 2011 - mass nargchk fix, fprintf fix
%        Feb.  3, 2012 - doc update
%        Jan. 26, 2014 - abs path exist fix
%        Mar. 10, 2014 - update See also section
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated mAR. 10, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(4,5,nargin));

% default overwrite to false
if(nargin==4 || isempty(o)); o=false; end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:writesacpz:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check z, p, k
if(~isvector(z) || ~isnumeric(z))
    error('seizmo:writesacpz:badZeros',...
        'ZEROS must be a vector of real/complex values!');
elseif(~isvector(p) || ~isnumeric(p))
    error('seizmo:writesacpz:badPoles',...
        'POLES must be a vector of real/complex values!');
elseif(~isreal(k))
    error('seizmo:writesacpz:badConstant',...
        'CONSTANT must be a real scalar!');
end

% directory separator
fs=filesep;

% graphical selection
if(isempty(file))
    [file,path]=uiputfile(...
        {'SAC_PZs_*;sac_pzs_*' 'SAC_PZs Files (SAC_PZs_*,sac_pzs_*)';
        '*.*' 'All Files (*.*)'},...
        'Save SAC PoleZero File as');
    if(isequal(0,file))
        error('seizmo:writesacpz:noFileSelected',...
            'No output file selected!');
    end
    file=[path fs file];
else
    % check file
    if(~isstring(file))
        error('seizmo:writesacpz:fileNotString',...
            'FILE must be a string!');
    end
    if(~isabspath(file)); file=[pwd fs file]; end
    if(exist(file,'file'))
        if(exist(file,'dir'))
            error('seizmo:writesacpz:dirConflict',...
                'SAC PoleZero File: %s\nIs A Directory!',file);
        end
        if(~o)
            fprintf('SAC PoleZero File: %s\nFile Exists!\n',file);
            reply=input('Overwrite? Y/N [N]: ','s');
            if(isempty(reply) || ~strncmpi(reply,'y',1))
                disp('Not overwriting!');
                return;
            end
            disp('Overwriting!');
        end
    end
end

% get total number of poles/zeros
nz=numel(z);
np=numel(p);

% remove poles/zeros at origin
z(abs(z)==0)=[];
p(abs(p)==0)=[];

% get new number of poles/zeros
nnz=numel(z);
nnp=numel(p);

% open file for writing as ascii
fid=fopen(file,'wt');

% check if file is openable
if(fid<0)
    error('seizmo:writesacpz:cannotOpenFile',...
        'SAC PoleZero File: %s\nNot Openable!',file);
end

% write to file
fprintf(fid,'ZEROS %d\n',nz);
if(nnz); fprintf(fid,'%g  %g\n',[real(z(:).'); imag(z(:).')]); end
fprintf(fid,'POLES %d\n',np);
if(nnp); fprintf(fid,'%g  %g\n',[real(p(:).'); imag(p(:).')]); end
fprintf(fid,'CONSTANT %e\n',k);

% close file
fclose(fid);

end
