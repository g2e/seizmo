function [data]=syncsr(data,sr)
%SYNCSR    Resamples SAClab data records
%
%    Description: Resamples data records using Matlab's resample routine if
%     possible, otherwise the problem is passed on to intrpol8.  Note that
%     resample uses a cascade of fir filtering interpolations and 
%     decimatations to get to the new rate.  This means that edge effects
%     are an issue if the record is not demeaned/detrended/tapered.  This
%     also means that aliasing is avoided.  In the case where resample 
%     fails, the signal is pre-decimated and then interpolated to avoid
%     aliasing.  Syncsr does not work correctly with unevenly sampled
%     records (but will not check for them).
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
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% check rate
if(~isnumeric(sr) || ~isscalar(sr))
    error('sample rate must be numeric scalar')
elseif(sr<=0)
    error('negative sample rate not allowed')
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
data=distro(data,recs,ind,store,nnpts);

% update header
data=ch(data,'delta',1/sr);
data=checkup(data);

end
