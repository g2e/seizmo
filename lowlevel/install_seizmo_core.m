function [ok]=install_seizmo_core(mypath,varargin)
%INSTALL_SEIZMO_CORE    Minimal SEIZMO install
%
%    Usage:    ok=install_seizmo_core(mypath)
%              ok=install_seizmo_core(mypath,type)
%
%    Description:
%     OK=INSTALL_SEIZMO_CORE(MYPATH) installs the bare minimum SEIZMO
%     directories to the path.  The path is saved to the pathdef.m file
%     that the path was loaded from at Matlab startup.  The "bare minimum"
%     means read/write support, header editing, struct editing, and
%     common utilities.
%
%     OK=INSTALL_SEIZMO_CORE(MYPATH,TYPE) allows editing the pathdef.m
%     save preference.  See SAVEPATH_SEIZMO for details.  The default is
%     no input to SAVEPATH_SEIZMO.
%
%    Notes:
%
%    Examples:
%
%    See also: INSTALL_SEIZMO, INSTALL_SEIZMO_MATTAUP, INSTALL_SEIZMO_WW3,
%              INSTALL_SEIZMO_OPTIONAL, INSTALL_SEIZMO_MMAP

%     Version History:
%        Dec. 30, 2010 - initial version
%        Jan.  2, 2011 - changed name from minimal to core, update for new
%                        directory layout
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  2, 2011 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% default path input
if(nargin<1 || isempty(mypath)); mypath='.'; end

% check path
fs=filesep;
if(~exist(mypath,'dir'))
    error('seizmo:install_seizmo_core:badPath',...
        ['SEIZMO directory (' mypath ') does not exist!']);
end

% add required seizmo components to path
addpath(...
    mypath,...
    [mypath fs 'lowlevel'],...
    [mypath fs 'rw'],...
    [mypath fs 'hdr'],...
    [mypath fs 'sz'],...
    [mypath fs 'misc'],...
    [mypath fs 'position'],...
    [mypath fs 'time']);
bad=savepath_seizmo(varargin{:});
if(~bad); ok=true; end

end
