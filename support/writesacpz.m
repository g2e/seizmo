function []=writesacpz(file,z,p,k)
%WRITESACPZ    Writes out a SAC PoleZero file
%
%    Description: WRITESACPZ(FILE,Z,P,K) writes the SAC PoleZero file FILE
%     using the zeros in Z, the poles in P, and the constant in K.  See the
%     Notes section for file format details.  Z and P must be numeric
%     vectors composed of reals and/or complex conjugate pairs.  K must be
%     a real scalar.
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
%    Tested on: Matlab r2007b
%
%    Usages:    writesacpz(file,z,p,k)
%
%    Examples:
%     Read in a SAC PoleZero file, alter the constant, and write out:
%      [z,p,k]=readsacpz('SAC_PZs_XB_CM32_BHZ_02');
%      k=k*correction_factor;
%      writesacpz('SAC_PZs_XB_CM32_BHZ_02',z,p,k);
%
%    See also: readsacpz, findsacpz, readresp, writeresp, findresp

%     Version History:
%        Apr.  7, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  7, 2009 at 02:45 GMT

% todo:

% check nargin
error(nargchk(4,4,nargin));

% check file
if(~ischar(file))
    error('seizmo:readsacpz:fileNotString',...
        'FILE must be a string!');
end
if(exist(file,'file'))
    if(exist(file,'dir'))
        error('seizmo:writesacpz:dirConflict',...
            'SAC PoleZero File: %s\nIs A Directory!',file);
    else
        warning('seizmo:writesacpz:fileClobber',...
            'SAC PoleZero File: %s\nAlready Exists! Overwriting!',...
            file);
    end
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
if(fid<0); error('SAC PoleZero File: %s\nNot Openable!',file); end

% write to file
fprintf(fid,'ZEROS %d\n',nz);
for i=1:nnz
    fprintf(fid,'%g  %g\n',real(z(i)),imag(z(i)));
end
fprintf(fid,'POLES %d\n',np);
for i=1:nnp
    fprintf(fid,'%g  %g\n',real(p(i)),imag(p(i)));
end
fprintf(fid,'CONSTANT %g\n',k);

% close file
fclose(fid);

end