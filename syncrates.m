function [data]=syncsr(data,sr)
%SYNCSR    Resamples SAClab data records to a common sample rate
%
%    Description: SYNCSR(DATA,SR) syncronizes the sample rates of data 
%     records in DATA to a new sample rate SR.  A fir filter is used to
%     avoid aliasing issues but can cause edge effects if the records
%     deviate from zero strongly at the start/end of the record.  Typically
%     using RTREND and TAPER helps to limit these effects.  Uses the Matlab
%     function resample (Signal Processing Toolbox).
%
%    Header Changes: DELTA, NPTS, DEPMEN, DEPMIN, DEPMAX, E
%
%    Notes:
%     - Requires evenly sampled data (use INTERPOLATE for uneven data)
%
%    Usage:  data=syncsr(data,sr)
%
%    Examples:
%     change all records to 5sps
%      data=syncsr(data,5)
%
%    See also: interpolate, iirfilter, deci, stretch

%     Version History:
%        ????????????? - Initial Version
%        June 16, 2008 - Updated documentation, rejects spectral records
%                        fixed no header update bug, and changed the way
%                        records are resampled when Matlab's resample fails
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 16, 2008 at 03:15 GMT

% check nargin
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data,'x'))

% check rate
if(~isnumeric(sr) || ~isscalar(sr) || sr<=0)
    error('SAClab:syncsr:badInput',...
        'sample rate must be a positive numeric scalar')
end

% require evenly sampled records only
if(any(~strcmp(glgc(data,'leven'),'true')))
    error('SAClab:syncsr:evenlySpacedOnly',...
        'illegal operation on unevenly spaced data');
end

% require non-spectral records
iftype=genumdesc(data,'iftype');
if(any(strcmp(iftype,'Spectral File-Real/Imag') | ...
        strcmp(iftype,'Spectral File-Ampl/Phase')))
    error('SAClab:syncsr:illegalOperation',...
        'Illegal operation of spectral file')
end

% make delta groups
delta=gh(data,'delta');
gdt=unique(delta); ng=length(gdt); gi=cell(ng,1);
for i=1:ng; gi{i}=find(delta==gdt(i)).'; end

% find fraction numerator/denominator of sampling
% rate ratio expressed as integers
[n,d]=rat(gdt*sr);

% work on each group
for i=1:ng
    % try resample
    try
        data(gi{i})=chsr(data(gi{i}),n(i),d(i),sr);
    catch
        disp('interpolating some records because resample failed...')
        % interpolate to a sample rate at least 4 times current AND desired
        % sample rates to preserve the data and allow for decimation.
        % here we find the lowest multiple of the desired rate that
        % satisfies that criteria...hopefully that integer is fairly
        % factorable...if not the function crashes
        lowest_multiple=max([4 ceil(4/(gdt*sr))]);
        
        % interpolate temporarily to high rate
        data(gi{i})=interpolate(data(gi{i}),lowest_multiple*sr,'spline');
        
        % now resample
        data(gi{i})=chsr(data(gi{i}),1,lowest_multiple,sr);
    end
end

end

function [data]=chsr(data,n,d,sr)
%CHSR  Resample records to new rate

% combine records
[recs,ind,store,npts]=combo(data);

% resample
recs=resample(recs,n,d);

% new length
nnpts=floor((npts-1)*n/d)+1;

% redistribute records
data=ch(data,'delta',1/sr);
data=distro(data,recs,ind,store,nnpts);

% update header
data=chkhdr(data);

end
