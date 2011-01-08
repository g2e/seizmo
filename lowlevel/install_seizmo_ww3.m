function [ok]=install_seizmo_ww3(mypath)
%INSTALL_SEIZMO_WW3    Compiles WaveWatch III mex files
%
%    Usage:    ok=install_seizmo_ww3(mypath)
%
%    Description:
%     OK=INSTALL_SEIZMO_WW3(MYPATH) compiles the mex files in the 'ww3'
%     directory within the MYPATH directory.  OK indicates if the
%     compilation succeeded or not.
%
%    Notes:
%     - Currently only BDS_unpack_mex5.c is compiled
%
%    Examples:
%
%    See also: INSTALL_SEIZMO, INSTALL_SEIZMO_MMAP, INSTALL_SEIZMO_MATTAUP,
%              INSTALL_SEIZMO_CORE, INSTALL_SEIZMO_OPTIONAL

%     Version History:
%        Dec. 30, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 30, 2010 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% default path input
if(nargin<1 || isempty(mypath)); mypath='.'; end

% check path
fs=filesep;
ww3path=[mypath fs 'ww3'];
if(~exist(ww3path,'dir'))
    error('seizmo:install_seizmo_ww3:badPath',...
        ['Expected path to WW3 directory (' ww3path ') does not exist!']);
end

% get to ww3 directory
cwd=pwd;
cd([mypath fs 'ww3']);

% compile mex file(s)
ok=false;
try
    mex('BDS_unpack_mex5.c');
    ok=true;
catch
    file=[mypath fs 'ww3' fs 'BDS_unpack_mex5.c'];
    warning('seizmo:install_seizmo_ww3:compilationFailed',...
        ['Failed to compile WaveWatch III mex file:\n' ...
        strrep(file,'\','\\') '\n' ...
        'WaveWatch III & SEIZMO interaction will fail until resolved!']);
end

% return to starting directory
cd(cwd);

end
