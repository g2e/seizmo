function [data]=noise_stack_delaz(data,dstep,azstep,ztrans)
%NOISE_STACK_DELAZ    Stacks NCFs based on distance/azimuth
%
%    Usage:    ncfs=noise_stack_delaz(ncfs,dstep,azstep)
%              ncfs=noise_stack_delaz(ncfs,dstep,azstep,ztransform)
%
%    Description:
%     NCFS=NOISE_STACK_DELAZ(NCFS,DSTEP,AZSTEP) stacks the seismic noise
%     correlations in the SEIZMO dataset NCFS by distance and azimuth bins
%     DSTEP and AZSTEP wide in degrees.  Bin intervals start at 0 and only
%     those that contain records are returned.  By default DSTEP is 0.01
%     degrees and AZSTEP is 360 degrees so these options may be omitted.
%
%     NCFS=NOISE_STACK_DELAZ(NCFS,DSTEP,AZSTEP,ZTRANSFORM) sets if Fisher's
%     transform is applied to the correlations for stacking.  If ZTRANSFORM
%     is TRUE Fisher's transform is applied (the default) while ZTRANSFORM
%     set to FALSE will do stacking the correlation values as is.
%
%    Notes:
%     - This differs from NOISE_STACK* in that it focuses on combining NCFs
%       from many different station pairings by distance & azimuth ranges
%       rather than only combining the same station pairing across a range
%       of time.  This technique is useful for improving the estimates of
%       the empirical Green's function for a region in the presence of an
%       anisotropic seismic noise distribution by averaging station pairs
%       of similar distance but different azimuths.
%     - Header info of the output records (noise correlation stacks)
%       related to master and slave station positions & times is
%       representative of only one record going into each stack.  Therefore
%       this info should be ignored.  The header fields RESP6/7 contain the
%       degree distance range and RESP8/9 contain the degree azimuth range
%       for the stack bin.
%
%    Header changes: SCALE (number of records in stack), DEP*, RESP6-9
%
%    Examples:
%     % Typically this is used to average NCFs from station pairs with
%     % similar interstation distance.  Tolerances in the distance range
%     % are best decided by the target frequency range (low frequencies
%     % allow stacking NCFs over a larger distance range):
%     ncfs=readseizmo('stack/zz_all_twoway/*/*');
%     ancfs=noise_stack_delaz(ncfs,0.1);
%
%    See also: NOISE_STACK, NOISE_STACK_ARBITRARY, STACK2STACK

%     Version History:
%        May  29, 2014 - initial version
%        June  4, 2014 - minor doc update
%        June 23, 2014 - minor doc update
%        July 11, 2014 - fd i/o, bugfix: addrecords hack now correct
%        July 17, 2014 - bugfix: solofun needs func handles not strings,
%                        bugfix: convert iftype to string
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 17, 2014 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% defaults
if(nargin<2 || isempty(dstep)); dstep=0.01; end
if(nargin<3 || isempty(azstep)); azstep=360; end
if(nargin<4 || isempty(ztrans)); ztrans=true; end

% require records are correlations
[kuser0,kuser1,iftype]=getheader(data,...
    'kuser0','kuser1','iftype id');
xc=strcmp(kuser0,'MASTER') & strcmp(kuser1,'SLAVE');
if(~all(xc))
    error('seizmo:noise_stack_delaz:badInput',...
        'NCFS does contains non-correlations!');
elseif(numel(unique(iftype))~=1)
    error('seizmo:noise_stack_delaz:badInput',...
        'NCFS contains mixed correlation filetypes!');
elseif(~isnumeric(dstep) || ~isscalar(dstep) || dstep<=0)
    error('seizmo:noise_stack_delaz:badInput',...
        'DSTEP must be a positive scalar!');
elseif(~isnumeric(azstep) || ~isscalar(azstep) || azstep<=0 ...
        || mod(360,azstep)~=0)
    error('seizmo:noise_stack_delaz:badInput',...
        'AZSTEP must be a positive scalar that evenly divides 360!');
elseif(~islogical(ztrans) || ~isscalar(ztrans))
    error('seizmo:noise_stack_delaz:badInput',...
        'ZTRANSFORM must be either TRUE or FALSE!');
end

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);
oldcheckheaderstate=checkheader_state(false);

% attempt stacking
try
    % for debugging
    for i=1:numel(data)
        data(i).misc.stacknames={[data(i).path data(i).name]};
    end
    
    % get some header info
    [gc,az,scale]=getheader(data,'gcarc','az','scale');
    
    % convert fd to cplx
    iftype=iftype{1};
    switch lower(iftype)
        case 'irlim'
            data=solofun(data,@(x)complex(x(:,1),x(:,2)));
        case 'iamph'
            data=solofun(data,@(x)x(:,1).*exp(1j*x(:,2)));
    end
    
    % fisher transform data
    if(ztrans); data=solofun(data,@fisher); end
    
    % multiply by scale
    % - this allows weighted stacks
    scale(isnan(scale))=1;
    data=multiply(data,scale);
    
    % force azimuth to be 0-360
    az=lonmod(az);
    az(az<0)=az(az<0)+360;
    
    % find the minimum & maximum distance intervals
    mingc=fix(min(gc)/dstep)*dstep;
    maxgc=fix(max(gc)/dstep)*dstep;
    
    % number of records and stacks
    nrecs=numel(data);
    stacks=0;
    
    % preallocate bin info for headers
    [dmin,dmax,amin,amax,sscale]=deal(nan(nrecs,1));
    
    % distance bin loop
    for d=mingc:dstep:maxgc
        % azimuth bin loop
        for a=0:azstep:360
            % records within the bin
            in=gc>=d & gc<d+dstep & az>=a & az<a+azstep;
            
            % only stack if there are records within the bin
            if(any(in))
                % found another bin with data
                stacks=stacks+1;
                
                % stack or copy based on number in bin
                if(sum(in)>1)
                    data(nrecs+stacks)=addrecords(data(in),'ref','ignore');
                    sscale(stacks)=sum(scale(in));
                else % only 1 record
                    data(nrecs+stacks)=data(find(in,1));
                    sscale(stacks)=scale(in);
                end
                
                % d range & az range go somewhere too
                dmin(stacks)=d;
                dmax(stacks)=d+dstep;
                amin(stacks)=a;
                amax(stacks)=a+azstep;
            end
        end
    end
    
    % remove original data
    data=data(nrecs+1:end);
    
    % trim excess off header input
    dmin=dmin(1:stacks);
    dmax=dmax(1:stacks);
    amin=amin(1:stacks);
    amax=amax(1:stacks);
    
    % add header info
    data=changeheader(data,'resp6',dmin,'resp7',dmax,...
        'resp8',amin,'resp9',amax);
    
    % divide by scale to get back to an average
    % - updates dep* stats skipped by addrecords hack
    data=divide(data,sscale);
    
    % inverse fisher transform data
    if(ztrans); data=solofun(data,@ifisher); end
    
    % convert cplx to fd
    % - updates dep* to not be complex
    switch lower(iftype)
        case 'irlim'
            data=solofun(data,@(x)[real(x),imag(x)]);
        case 'iamph'
            data=solofun(data,@(x)[abs(x),angle(x)]);
    end
    
    % rename
    data=changename(data,'name',strcat('DEL_',num2str(dmin),'-',...
        num2str(dmax),'_AZ_',num2str(amin),'-',num2str(amax)));
    
    % update headers
    data=changeheader(data,'scale',sscale);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

% toggle checking back
seizmocheck_state(oldseizmocheckstate);
checkheader_state(oldcheckheaderstate);

end


function [d1]=addrecords(d1,varargin)
% simple hack for speed (no dep* update)
try
    d1(1).dep=mean(cat(3,d1.dep),3);
    d1(1).misc.stacknames=getsubfield(d1,'misc','stacknames');
    d1(2:end)=[];
    d1.misc.stacknames=[d1.misc.stacknames{:}]';
catch
    error('seizmo:noise_stack_arbitrary:badNCFs',...
        'NCFs differ in number of points! Cannot stack!');
end
end

