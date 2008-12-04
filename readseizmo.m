function [data]=readseizmo(varargin)
%READSEIZMO    Read datafiles into SEIZMO data structure
%
%    Description: READSEIZMO(FILELIST1,...,FILELISTN) reads compatible 
%     datafiles into a SEIZMO data structure.  Accepts character arrays of 
%     filenames (one filename per row) and/or cell arrays of filenames (one
%     filename per cell).  Wildcards are valid.
%
%     SEIZMO data structure setup:
%
%     Fields for all files:
%      path - path to file
%      name - file name
%      filetype - type of datafile
%      version - version of filetype
%      byteorder - byte-order of file (ieee-le or ieee-be)
%      hasdata - logical indicating if data is read in
%      head - contains header data
%
%     Fields for timeseries files:
%      dep(:,1) - amplitudes
%      ind(:,1) - times (if uneven spacing)
%
%     Fields for spectral amp/phase files:
%      dep(:,1) - spectral amplitudes
%      dep(:,2) - spectral phase
%
%     Fields for spectral real/imag files:
%      dep(:,1) - spectral real
%      dep(:,2) - spectral imaginary
%
%     Fields for general xy files:
%      dep(:,1) - dependent component
%      ind(:,1) - independent component (if uneven spacing)
%
%     Fields for xyz grid files:
%      dep(:,1) - matrix data (x and y evenly spaced; z steps l2r,b2t)
%     
%    Notes:
%     - Multi-component files will replicate the number of columns by the
%       number of components.  So a three component spectral file will have
%       six columns of data.  Dependent components ('dep') share the same 
%       independent component ('ind').
%
%    Tested on: Matlab r2007b
%
%    Header changes: see CHECKHEADER
%
%    Usage:    data=readseizmo(filelist1,...,filelistN)
%
%    Examples:
%     Some simple read statements:
%      data=readseizmo('KATH.R');
%      data=readseizmo('SQRL.R','AAK.R');
%
%     Read in every SEIZMO compatible datafile from the current directory:
%      data=readseizmo('*');
%
%    See also: readheader, readdata, readdatawindow, writeheader, 
%              writeseizmo, bseizmo, changeheader, getheader, listheader

%     Version History:
%        May  30, 2007 - auto-determine byte-order, multiple file support
%        Oct. 29, 2007 - complete rewrite, data now stored as struct
%        Nov.  7, 2007 - doc update
%        Jan. 27, 2008 - complete rewrite, SACHP support
%        Feb. 11, 2008 - code now in RH and RDATA (now just a wrapper)
%        Feb. 29, 2008 - renamed from RSAC to RSEIS
%        Mar.  4, 2008 - doc update
%        Apr. 23, 2008 - wildcard filename support (doc update)
%        Sep. 14, 2008 - doc update
%        Sep. 15, 2008 - history fix
%        Oct. 15, 2008 - doc update for hasdata field
%        Oct. 16, 2008 - doc update for dir and filetype fields
%        Oct. 27, 2008 - minor doc update for struct change
%        Nov. 15, 2008 - update for new name schema (now READSEIZMO)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 15, 2008 at 23:45 GMT

% todo:

% call specific routines
data=readheader(varargin{:});
data=readdata(data);

end
