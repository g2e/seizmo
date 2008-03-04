function [data]=syncsr(data,sr)
%SYNCSR    Resamples SAClab data records
%
%    Description: Resamples data records using Matlab's resample routine if
%     possible, otherwise the problem is passed on to intrpol8.  Note that
%     resample uses a cascade of fir filtering interpolations and 
%     decimatations to get to the new rate which will prevent aliasing but
%     may introduce edge effects if the record is not demeaned/detrended/
%     tapered.  If resample fails, the signal is passed to a simpler series
%     of commands that pre-decimate and then interpolate (using intrpol8).
%
%    Notes:
%     - Requires evenly sampled data
%
%    Usage:  [data]=syncsr(data,sr)
%
%    Examples:
%     change all records to 5sps
%      data=syncsr(data,5)
%
%    See also: intrpol8, iirfilter, deci, stretch

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
        disp('Interpolating because resample failed...')
        if(n(i)/d(i)<1)
            % change sampling rate to a lower rate if downsampling to avoid
            % aliasing in the output records
            deci=ceil(d(i)/n(i));
            data(gi{i})=chsr(data(gi{i}),1,deci,1/deci/gdt(i));
        end
        % interpolate to desired rate
        data(gi{i})=intrpol8(data(gi{i}),sr);
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

end
