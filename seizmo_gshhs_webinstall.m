function [ok]=seizmo_gshhs_webinstall(mypath,varargin)
%SEIZMO_GSHHS_WEBINSTALL    Install GSHHS components for SEIZMO
%
%    Usage:    ok=seizmo_gshhs_webinstall
%              ok=seizmo_gshhs_webinstall(path)
%              ok=seizmo_gshhs_webinstall(path,type)
%
%    Description:
%     OK=SEIZMO_GSHHS_WEBINSTALL downloads the GSHHS zip file to the
%     current directory, extracts its contents to the 'gshhs' directory and
%     installs them on the Matlab/Octave path.  THE DOWNLOAD IS LARGE:
%     OVER 100 MEGABYTES & CAN TAKE HOURS ON A SLOW CONNECTION!  THERE IS
%     NO RESUME EITHER SO IF IT FAILS, YOU WILL LIKELY LOSE EVERYTHING YOU
%     DOWNLOADED!  See the notes below on how to install the GSHHS binaries
%     manually.
%
%     OK=SEIZMO_GSHHS_WEBINSTALL(PATH) install the GSHHS binary files under
%     the directory given by PATH.
%
%     OK=SEIZMO_GSHHS_WEBINSTALL(PATH,TYPE) allows editing the pathdef.m
%     save preference.  See SAVEPATH_SEIZMO for details.  The default is
%     no input to SAVEPATH_SEIZMO.
%
%    Notes:
%     - Get the GSHHS binary files by downloading the archive from here:
%        http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/oldversions/
%       The file is gshhs_1.10.zip under the 1.10 version.  Extract the
%       binaries from the archive using your unzip utility of choice.
%     - If you have the GSHHS binary files extracted, just add the
%       directory containing them to your Matlab or Octave path using:
%        addpath directory/of/gshhs/
%       where the directory/of/gshhs needs to be changed to where the files
%       are on your system.
%
%    Examples:
%     % GSHHS is big: this download can take _hours_.  This function is
%     % only useful if you have a fast connection!
%
%    See also: INSTALL_SEIZMO_MMAP, M_MAP, M_GSHHS

%     Version History:
%        Feb.  5, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  5, 2011 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% default path input
if(nargin<1 || isempty(mypath)); mypath='.'; end

% check path
if(~exist(mypath,'dir'))
    error('seizmo:install_seizmo_mmap:badPath',...
        ['SEIZMO directory (' mypath ') does not exist!']);
end

% attempt gshhs install
try
    cwd=pwd;
    cd(mypath);
    url='http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/';
    url=[url 'oldversions/version1.10/gshhs_1.10.zip'];
    urlwrite(url,'gshhs_1.10.zip');
    unzip('gshhs_1.10.zip');
    addpath('gshhs');
    savepath_seizmo(varargin{:});
    cd(cwd);
    ok=true;
catch
    ok=false;
    cd(cwd);
end

end
