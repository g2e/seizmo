function [data]=rseis(varargin)
%RSEIS    Read datafiles into SAClab data structure
%
%    Description: RSEIS(FILELIST1,...,FILELISTN) reads datafiles into a 
%     SAClab data structure.  Accepts character arrays of filenames (one 
%     filename per row) and/or cell arrays of filenames (one filename per 
%     cell).  Wildcards are valid.
%
%     SAClab data structure setup:
%
%     Fields for all files:
%      head - contains header data
%      name - filename (may include path)
%      endian - byte-order of file (ieee-le or ieee-be)
%      version - version of datafile
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
%       six columns of data.  Dependent omponents share the same 
%       independent component.
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: arguments must be char or cellstr arrays
%
%    Header changes: NONE
%
%    Usage:    data=rseis(filelist1,...,filelistN)
%
%    Examples:
%     Some simple read statements:
%      data=rseis('KATH.R');
%      data=rseis('SQRL.R','AAK.R');
%
%     Read in every SAClab readible datafile from the current directory:
%      data=rseis('*');
%
%    See also: rh, rdata, rpdw, wh, wseis, bseis

%     Version History:
%        Sometime 2006 - 0.1.0 - auto-determine byte-order
%        Early    2007 - 0.2.0 - multiple file support
%        Nov.  7, 2007 - 0.3.0 - introduced datafile definition
%        Jan. 27, 2008 - 0.4.0 - slicker datafile definition structure
%                        0.5.0 - 'best' datafile definition structure
%        Feb. 29, 2008 - renamed from rsac to rseis
%        Mar.  4, 2008 - doc update
%        Apr. 23, 2008 - wildcard filename support
%        Sep. 14, 2008 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2008 at 22:00 GMT

% todo:

% call specific routines
data=rh(varargin{:});
data=rdata(data);

end
