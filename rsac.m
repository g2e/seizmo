function [data]=rsac(varargin)
%RSAC    Read SAC binary files
%
%    Description: Reads in SAC (seismic analysis code) binary format files 
%     into a data structure.  Accepts character arrays of filenames (one 
%     filename per row) and/or cell arrays of filenames (one filename per 
%     cell).
%
%     Fields for all SAC files:
%      data.head - header
%      data.name - filename (may include path)
%      data.endian - byte-order of file
%      data.version - header version of file
%
%     Fields for SAC timeseries files:
%      data.x(:,1) - amplitudes
%      data.t - times (if uneven spacing)
%
%     Fields for SAC spectral amp/phase files:
%      data.x(:,1) - spectral amplitudes
%      data.x(:,2) - spectral phase
%
%     Fields for SAC spectral real/imag files:
%      data.x(:,1) - spectral real
%      data.x(:,2) - spectral imaginary
%
%     Fields for SAC general xy files:
%      data.x(:,1) - dependent component
%      data.t - independent component (if uneven spacing)
%
%     Fields for SAC xyz grid files:
%      data.x(:,1) - matrix data (x and y evenly spaced; z steps l2r,b2t)
%     
%     Fields for SAC n-component timeseries files:
%      data.x(:,1) - 1st component
%      data.x(:,2) - 2nd component
%      data.x(:,n) - nth component
%      data.t - times (if uneven spacing)
%
%    Usage:    saclab_struct=rsac(['sacfile1'; 'sacfile2'; ...],...
%                                 {'sacfile3' 'sacfile4' ...},...
%                                 'sacfile5','sacfile6',...)
%
%    Examples:
%     data=rsac('KATH.R');
%     data=rsac('SQRL.R','AAK.R');
%     data=rsac(cellarray1,chrarray1,chrarray2,cellarray2);
%
%     read in files from current directory
%      files=dir();
%      data=rsac(files.name);
%
%    By: Garrett Euler (2/2008)   ggeuler@wustl.edu
%
%    See also: rh, rdata, wsac, bsac, sachp, gv, rpdw, wh, gh, lh, ch

% call specific routines
data=rh(varargin{:});
data=rdata(data);

end
