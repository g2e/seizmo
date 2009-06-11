function [z,p,k]=readsacpz(file)
%READSACPZ    Reads in a SAC PoleZero file
%
%    Usage:    [z,p,k]=readsacpz(file)
%
%    Description: [Z,P,K]=READSACPZ(FILE) reads the SAC PoleZero file FILE,
%     returning the zeros in separate rows of Z, the poles in separate rows
%     of P, and the constant in K.  See the Notes section for file format
%     details.  The output should be compatible with the ZPK output from
%     Signal Processing Toolbox commands like BUTTER.
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
%       component (thus there are 2 complex conjugate pairs).  The last
%       line gives the multiplicative factor.
%     - Comment lines may be added to the SAC PoleZero file by starting the
%       line with a '*' (asterisk).
%     - READSACPZ will read all valid sections -- if there are multiple
%       zeros, poles or constant sections then only the last one will be
%       kept.  Please be aware that this may be important if using programs
%       that append to SAC PoleZero files rather than overwriting them.
%
%    Examples:
%     Read in a SAC PoleZero file, and apply it to some records:
%      [z,p,k]=readsacpz('SAC_PZs_XB_CM32_BHZ_02');
%      [sos,g]=zp2sos(z,p,k);       % converting to a
%      fs=dfilt.df2tsos(sos,g);     % filter structure
%      data=seizmofilter(data,fs);
%
%     Now take a look at the details of the PoleZero filter:
%      fvtool(fs)
%
%    See also: writesacpz, findsacpz, readresp, writeresp, findresp

%     Version History:
%        Apr.  7, 2009 - initial version
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 11, 2009 - add asterisk comment support, octave support
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2009 at 07:55 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end;

% check file
if(~ischar(file))
    error('seizmo:readsacpz:fileNotString',...
        'FILE must be a string!');
end
if(~exist(file,'file'))
    error('seizmo:readsacpz:fileDoesNotExist',...
        'SAC PoleZero File: %s\nDoes Not Exist!',file);
elseif(exist(file,'dir'))
    error('seizmo:writesacpz:dirConflict',...
        'SAC PoleZero File: %s\nIs A Directory!',file);
end

% open file for reading as ascii
fid=fopen(file,'rt');

% check if file is openable
if(fid<0); error('SAC PoleZero File: %s\nNot Openable!',file); end

% set defaults
z=[]; p=[]; k=1;

% read in all lines
c=0; a=cell(1,0);
while 1
    tmp=fgetl(fid);
    if(~ischar(tmp)); break; end
    c=c+1; a{c}=tmp;
end

% close file
fclose(fid);

% keep processing until through all lines
line=1; nlines=numel(a);
while(line<=nlines)
    % skip line if blank
    if(isempty(a{line}))
        line=line+1;
        continue;
    end
    
    % process line
    words=getwords(a{line});
    
    % skip line if blank or comment
    if(isempty(words) || strcmp(words{1}(1),'*'))
        line=line+1;
        continue;
    end
    
    % check number of words
    if(numel(words)~=2)
        error('seizmo:readsacpz:notSACPZ',...
            ['SAC PoleZero File: %s\nLine %d: %s\n'...
            'Does Not Conform To SAC PoleZero Format!'],...
            file,line,a{line});
    end
    
    % separate words
    [field,value]=deal(words{:});
    value=str2double(value);
    
    % check value
    if(isnan(value))
        error('seizmo:readsacpz:notSACPZ',...
            ['SAC PoleZero File: %s\nLine %d: %s\n'...
            'Does Not Conform To SAC PoleZero Format!'],...
            file,line,a{line});
    end
    
    % act based on line
    switch lower(field)
        case 'zeros'
            % check value is fixed
            if(value~=fix(value))
                error('seizmo:readsacpz:notSACPZ',...
                    ['SAC PoleZero File: %s\nLine %d: %s\n'...
                    'Does Not Conform To SAC PoleZero Format!'],...
                    file,line,a{line});
            end
            
            % preallocate zeros
            z=zeros(value,1);
            
            % read in zeros (not at origin)
            line=line+1; zero=1;
            while(line<=nlines && zero<=value)
                % skip line if blank
                if(isempty(a{line}))
                    line=line+1;
                    continue;
                end
                
                % process line
                words=getwords(a{line});
                
                % skip line if blank or comment
                if(isempty(words) || strcmp(words{1}(1),'*'))
                    line=line+1;
                    continue;
                end
                
                % check number of words
                if(numel(words)~=2)
                    error('seizmo:readsacpz:notSACPZ',...
                        ['SAC PoleZero File: %s\nLine %d: %s\n'...
                        'Does Not Conform To SAC PoleZero Format!'],...
                        file,line,a{line});
                end
                
                % separate words
                [value1,value2]=deal(words{:});
                value1=str2double(value1);
                value2=str2double(value2);
                
                % check values
                if(isnan(value1) || isnan(value2)); break; end
                
                % assign values
                if(value2) % complex
                    z(zero)=value1+value2*i;
                else % real
                    z(zero)=value1;
                end
                
                % increment
                zero=zero+1;
                line=line+1;
            end
        case 'poles'
            % check value is fixed
            if(value~=fix(value))
                error('seizmo:readsacpz:notSACPZ',...
                    ['SAC PoleZero File: %s\nLine: %s\n'...
                    'Does Not Conform To SAC PoleZero Format!'],...
                    file,line,a{line});
            end
            
            % preallocate poles
            p=zeros(value,1);
            
            % read in poles (not at origin)
            line=line+1; pole=1;
            while(line<=nlines && pole<=value)
                % skip line if blank
                if(isempty(a{line}))
                    line=line+1;
                    continue;
                end
                
                % process line
                words=getwords(a{line});
                
                % skip line if blank or comment
                if(isempty(words) || strcmp(words{1}(1),'*'))
                    line=line+1;
                    continue;
                end
                
                % check number of words
                if(numel(words)~=2)
                    error('seizmo:readsacpz:notSACPZ',...
                        ['SAC PoleZero File: %s\nLine %d: %s\n'...
                        'Does Not Conform To SAC PoleZero Format!'],...
                        file,line,a{line});
                end
                
                % separate words
                [value1,value2]=deal(words{:});
                value1=str2double(value1);
                value2=str2double(value2);
                
                % check values
                if(isnan(value1) || isnan(value2)); break; end
                
                % assign values
                if(value2) % complex
                    p(pole)=value1+value2*i;
                else % real
                    p(pole)=value1;
                end
                
                % increment
                pole=pole+1;
                line=line+1;
            end
        case 'constant'
            k=value;
            line=line+1;
        otherwise
            error('seizmo:readsacpz:notSACPZ',...
                ['SAC PoleZero File: %s\nLine %d: %s\n'...
                'Does Not Conform To SAC PoleZero Format!'],...
                file,line,a{line});
    end
end

end
