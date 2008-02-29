function [data]=rseis(varargin)
%RSEIS    Read seismic binary datafiles
%
%    Description: Reads in binary seismic datafiles into a data structure.  
%     Accepts character arrays of filenames (one filename per row) and/or 
%     cell arrays of filenames (one filename per cell).
%
%     Fields for all files:
%      data.head - header
%      data.name - filename (may include path)
%      data.endian - byte-order of file
%      data.version - header version of file
%
%     Fields for timeseries files:
%      data.x(:,1) - amplitudes
%      data.t - times (if uneven spacing)
%
%     Fields for spectral amp/phase files:
%      data.x(:,1) - spectral amplitudes
%      data.x(:,2) - spectral phase
%
%     Fields for spectral real/imag files:
%      data.x(:,1) - spectral real
%      data.x(:,2) - spectral imaginary
%
%     Fields for general xy files:
%      data.x(:,1) - dependent component
%      data.t - independent component (if uneven spacing)
%
%     Fields for xyz grid files:
%      data.x(:,1) - matrix data (x and y evenly spaced; z steps l2r,b2t)
%     
%    Notes:
%     - multi-component files will replicate the number of columns by the
%       number of components.  So a three component spectral file will have
%       six columns of data.
%
%    Usage:    seislab_struct=rseis(['seisfile1'; 'seisfile2'; ...],...
%                                 {'seisfile3' 'seisfile4' ...},...
%                                 'seisfile5','seisfile6',...)
%
%    Examples:
%     data=rseis('KATH.R');
%     data=rseis('SQRL.R','AAK.R');
%     data=rseis(cellarray1,chrarray1,chrarray2,cellarray2);
%
%     read in files from current directory
%      files=dir();
%      data=rseis(files.name);
%
%    See also: rh, rdata, rpdw, wh, wseis, bseis

% call specific routines
data=rh(varargin{:});
data=rdata(data);

end
