function [data]=rseis(varargin)
%RSEIS    Read binary seismic datafiles into SAClab data structure
%
%    Description: Reads binary seismic datafiles into a SAClab data 
%     structure.  Accepts character arrays of filenames (one filename per 
%     row) and/or cell arrays of filenames (one filename per cell).
%     Wildcards are valid.
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
%      x(:,1) - amplitudes
%      t(:,1) - times (if uneven spacing)
%
%     Fields for spectral amp/phase files:
%      x(:,1) - spectral amplitudes
%      x(:,2) - spectral phase
%
%     Fields for spectral real/imag files:
%      x(:,1) - spectral real
%      x(:,2) - spectral imaginary
%
%     Fields for general xy files:
%      x(:,1) - dependent component
%      t(:,1) - independent component (if uneven spacing)
%
%     Fields for xyz grid files:
%      x(:,1) - matrix data (x and y evenly spaced; z steps l2r,b2t)
%     
%    Notes:
%     - Multi-component files will replicate the number of columns by the
%       number of components.  So a three component spectral file will have
%       six columns of data.  Components share the same timing.
%
%    Usage:    SAClab_struct=rseis(['seisfile1'; 'seisfile2'; ...],...
%                                 {'seisfile3' 'seisfile4' ...},...
%                                 'seisfile5','seisfile6',...)
%
%    Examples:
%     data=rseis('KATH.R');
%     data=rseis('SQRL.R','AAK.R');
%     data=rseis(cellarray1,chrarray1,chrarray2,cellarray2);
%
%     Read in SAClab datafiles from current directory:
%      data=rseis('*');
%
%    See also: rh, rdata, rpdw, wh, wseis, bseis

% call specific routines
data=rh(varargin{:});
data=rdata(data);

end
