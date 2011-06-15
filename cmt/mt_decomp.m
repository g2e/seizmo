function [varargout]=mt_decomp(mt,option)
%MT_DECOMP    Decompose moment tensor(s)
%
%    Usage:    [iso,dev]=mt_decomp(mt,'iso')
%              [major,minor]=mt_decomp(mt,'majmin')
%              [major,middle]=mt_decomp(mt,'majmid')
%              [middle,minor]=mt_decomp(mt,'midmin')
%              [dblcpl,clvd]=mt_decomp(mt,'maxdc')
%              [clvd,dblcpl]=mt_decomp(mt,'maxclvd')
%              [clvd,dblcpl]=mt_decomp(mt,'clvd')
%              [...,vec]=mt_decomp(...)
%
%    Description:
%     [ISO,DEV]=MT_DECOMP(MT,'ISO') returns the isotropic and deviatoric
%     components of the moment tensor(s) given in MT.  MT must be a Nx6 or
%     3x3xN array where N is the number of moment tensors in MT.  Note that
%     all GlobalCMT & USGS moment tensors are constrained to be purely
%     deviatoric.
%
%     [MAJOR,MINOR]=MT_DECOMP(MT,'MAJMIN') returns major and minor double
%     couples for the moment tensor(s) in MT.  Note that this decomposition
%     is not unique: we use the largest and smallest eigenvalues (in an
%     absolute sense) to define the moments of the double couples but it is
%     just as valid to choose some other pairing of moments (see 'MAJMID' &
%     'MIDMIN' options).
%
%     [MAJOR,MIDDLE]=MT_DECOMP(MT,'MAJMID') returns major and middle double
%     couples for the moment tensor(s) in MT.  Note that this decomposition
%     is not unique: we use the largest and middle eigenvalues (in an
%     absolute sense) to define the moments of the double couples but it is
%     just as valid to choose some other pairing of moments (see 'MAJMIN' &
%     'MIDMIN' options).
%
%     [MIDDLE,MINOR]=MT_DECOMP(MT,'MIDMIN') returns middle and minor double
%     couples for the moment tensor(s) in MT.  Note that this decomposition
%     is not unique: we use the middle and minor eigenvalues (in an
%     absolute sense) to define the moments of the double couples but it is
%     just as valid to choose some other pairing of moments (see 'MAJMIN' &
%     'MAJMID' options).
%
%     [DBLCPL,CLVD]=MT_DECOMP(MT,'MAXDC') returns the maximum double couple
%     with a compensated linear-vector dipole (clvd) component for the
%     moment tensor(s) in MT.  Note that this decomposition is not unique:
%     it is as valid to choose to decompose into a maximum clvd with some
%     double-couple component (see 'MAXCLVD' & 'CLVD' options).  This is
%     the decomposition used by the GlobalCMT & USGS groups for defining
%     the fault planes associated with their moment tensors.
%
%     [CLVD,DBLCPL]=MT_DECOMP(MT,'MAXCLVD') returns a compensated linear
%     vector dipole (clvd) and double couple for the moment tensor(s) in
%     MT.  The largest moment of the clvd is set to the maximum deviatoric
%     moment in this case.  Note that this decomposition is not unique:
%     it is as valid to choose to decompose into a maximum double couple
%     with some clvd component (see 'MAXDC' & 'CLVD' options).
%
%     [CLVD,DBLCPL]=MT_DECOMP(MT,'CLVD') returns a compensated linear
%     vector dipole (CLVD) and double couple for the moment tensor(s) in
%     MT.  The CLVD and double couple share the T & P principal axes.  This
%     is also not a unique decomposition (see 'MAXDC' & 'MAXCLVD' for
%     dc+clvd alternatives).  This is the decomposition suggested by
%     Knopoff & Randall 1970 and is widely used in moment tensor analysis.
%
%     [...,VEC]=MT_DECOMP(...) additionally returns the eigenvectors of the
%     principal axes.  Use with MT_UNDIAG to plot the decomposed tensors.
%
%    Notes:
%     - References:
%        Knopoff & Randall 1970, JGR, V 75, pp. 4957-4963
%        Dziewonski, Chou, & Woodhouse 1981, JGR, V 86, pp. 2825-2852
%        Wallace 1985, JGR, V 90, pp. 11,171-11,176
%        Jost & Hermann 1989, SRL, V 60, pp. 37-57
%        Julian, Miller, & Foulger 1998, Rev Geoph, V 36, No 4, pp. 525-549
%
%    Examples:
%     % Decompose some cmts into a CLVD & best double-couples following
%     % Knopoff & Randall 1970, then undiagonalize and plot:
%     cmts=findcmt('n',10);
%     [clvd,best,vec]=mt_decomp(mt_s2v(cmts),'clvd');
%     clvd=mt_undiag(clvd,vec);
%     best=mt_undiag(best,vec);
%     figure;
%     plotmt(1:10,2,mt_s2v(cmts))
%     hold on;
%     plotmt(1:10,1,clvd)
%     plotmt(1:10,0,best)
%     axis equal tight off
%
%    See also: MT_DIAG, MT_UNDIAG, PLOTMT

%     Version History:
%        June 10, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 10, 2011 at 13:50 GMT

% todo:
% - 3 DC, 3 CLVD (see Jost & Hermann 1989)

% check nargin
error(nargchk(2,2,nargin));

% check tensor format (force 3x3xN)
mtsz=size(mt);
if(~isnumeric(mt) || ~isreal(mt))
    error('seizmo:mt_decomp:badInput',...
        'MT must be a real-valued numeric array!');
elseif(isequal(mtsz(1:2),[3 3]) && any(numel(mtsz)==[2 3]))
    if(numel(mtsz)>2); n=mtsz(3); else n=1; end
elseif(mtsz(2)==6 && numel(mtsz)==2)
    mt=mt_v2g(mt); % convert from Nx6 to 3x3xN
    n=mtsz(1);
else
    error('seizmo:mt_decomp:badInput',...
        'MT must be a harvard moment tensor array as 3x3xN or Nx6!');
end

% check option string
valid={'iso' 'majmin' 'majmid' 'midmin' 'maxdc' 'maxclvd' 'clvd'};
if(~any(strcmpi(option,valid)))
    error('seizmo:mt_decomp:badInput',...
        ['OPTION must be one of the following:\n' ...
        sprintf('''%s'' ',valid{:})]);
end

% first diagonalize (aka principal axis transformation)
% - eigenvalues are sorted from smallest to largest along diagonal
[mt,vec]=mt_diag(mt);

% now get isotropic and deviatoric
tr=nan(n,1); [iso,dev]=deal(nan(3,3,n));
for i=1:n
    tr(i)=trace(mt(:,:,i))/3;
    iso(:,:,i)=diag(tr([i i i]));
    dev(:,:,i)=mt(:,:,i)-iso(:,:,i);
end

% specific decompositions of deviatoric component
switch lower(option)
    case 'iso'
        varargout={iso dev vec};
    case 'majmin' % wallace 1985
        % largest abs eigenvalue = major dc moment
        % smallest abs eigenvalue = minor dc moment
        [major,minor]=deal(nan(3,3,n));
        for i=1:n
            % identify eigenvalues
            egn=diag(dev(:,:,i));
            hi=find(max(abs(egn))==abs(egn),1);
            lo=find(min(abs(egn))==abs(egn),1);
            md=6-hi-lo;
            
            % major
            tmp=zeros(3,1);
            tmp(hi)=egn(hi);
            tmp(md)=-egn(hi);
            major(:,:,i)=diag(tmp);
            
            % minor
            tmp=zeros(3,1);
            tmp(lo)=egn(lo);
            tmp(md)=-egn(lo);
            minor(:,:,i)=diag(tmp);
        end
        varargout={major minor vec};
    case 'majmid' % julian et al 1998
        % largest abs eigenvalue = major dc moment
        % middle abs eigenvalue = minor dc moment
        [major,minor]=deal(nan(3,3,n));
        for i=1:n
            % identify eigenvalues
            egn=diag(dev(:,:,i));
            hi=find(max(abs(egn))==abs(egn),1);
            lo=find(min(abs(egn))==abs(egn),1);
            md=6-hi-lo;
            
            % major
            tmp=zeros(3,1);
            tmp(hi)=egn(hi);
            tmp(lo)=-egn(hi);
            major(:,:,i)=diag(tmp);
            
            % minor
            tmp=zeros(3,1);
            tmp(md)=egn(md);
            tmp(lo)=-egn(md);
            minor(:,:,i)=diag(tmp);
        end
        varargout={major minor vec};
    case 'midmin' % wallace 1985
        % dc share same t axis
        % middle abs eigenvalue = major dc moment
        % smallest abs eigenvalue = minor dc moment
        [major,minor]=deal(nan(3,3,n));
        for i=1:n
            % identify eigenvalues
            egn=diag(dev(:,:,i));
            hi=find(max(abs(egn))==abs(egn),1);
            lo=find(min(abs(egn))==abs(egn),1);
            md=6-hi-lo;
            
            % major
            tmp=zeros(3,1);
            tmp(md)=egn(md);
            tmp(hi)=-egn(md);
            major(:,:,i)=diag(tmp);
            
            % minor
            tmp=zeros(3,1);
            tmp(lo)=egn(lo);
            tmp(hi)=-egn(lo);
            minor(:,:,i)=diag(tmp);
        end
        varargout={major minor vec};
    case '3dc' % Jost & Hermann 1989
        
    case 'maxdc' % dziewonski et al 1981
        % used by globalcmt, usgs catalogs (T, P, B axes shared but mixed)
        % best double couple is 1/2 the difference of largest pos & neg egn
        % clvd is 2 min(abs(egn))/max(abs(egn))
        [clvd,dblcpl]=deal(nan(3,3,n));
        for i=1:n
            % identify eigenvalues
            egn=diag(dev(:,:,i));
            mx=find(max(egn)==egn,1);
            mn=find(min(egn)==egn,1);
            
            % best double couple
            tmp=zeros(3,1);
            tmp(mx)=(egn(mx)-egn(mn))/2;
            tmp(mn)=-(egn(mx)-egn(mn))/2;
            dblcpl(:,:,i)=diag(tmp);
            
            % clvd (remainder)
            clvd(:,:,i)=dev(:,:,i)-dblcpl(:,:,i);
        end
        varargout={dblcpl clvd vec};
    case 'maxclvd' % wallace 1985
        % clvd is largest absolute eigenvalue
        % dc is half difference between two smallest eigenvalues
        [clvd,dblcpl]=deal(nan(3,3,n));
        for i=1:n
            % identify eigenvalues
            egn=diag(dev(:,:,i));
            hi=find(max(abs(egn))==abs(egn),1);
            %lo=find(min(abs(egn))==abs(egn),1);
            
            % clvd
            tmp=-.5*egn(hi)*ones(3,1);
            tmp(hi)=egn(hi);
            clvd(:,:,i)=diag(tmp);
            
            % best double couple (remainder)
            dblcpl(:,:,i)=dev(:,:,i)-clvd(:,:,i);
        end
        varargout={clvd dblcpl vec};
    case 'clvd' % knopoff & Randall 1970
        % dc & clvd have the same t & p axes
        % dc is the difference of 2 abs smallest egn
        % clvd is -2 * abs smallest egn
        [clvd,dblcpl]=deal(nan(3,3,n));
        for i=1:n
            % identify eigenvalues
            egn=diag(dev(:,:,i));
            hi=find(max(abs(egn))==abs(egn),1);
            lo=find(min(abs(egn))==abs(egn),1);
            md=6-hi-lo;
            
            % best double couple
            tmp=zeros(3,1);
            tmp(hi)=(egn(lo)-egn(md));
            tmp(md)=-(egn(lo)-egn(md));
            dblcpl(:,:,i)=diag(tmp);
            
            % clvd (remainder)
            clvd(:,:,i)=dev(:,:,i)-dblcpl(:,:,i);
            %tmp=egn(lo)*ones(3,1); % for
            %tmp(hi)=-2*egn(lo);    % verification
            %clvd(:,:,i)=diag(tmp); % purposes
        end
        varargout={clvd dblcpl vec};
end

end
