function [varargout]=mt_decomp(mt,option)
%MT_DECOMP    Decompose moment tensor(s)
%
%    Usage:    [iso,dev]=mt_decomp(mt,'iso')
%              [major,minor]=mt_decomp(mt,'majmin')
%              [major,middle]=mt_decomp(mt,'majmid')
%              [middle,minor]=mt_decomp(mt,'midmin')
%              [dc1,dc2,dc3]=mt_decomp(mt,'3dc')
%              [clvd1,clvd2,clvd3]=mt_decomp(mt,'3clvd')
%              [dblcpl,clvd]=mt_decomp(mt,'maxdc')
%              [clvd,dblcpl]=mt_decomp(mt,'maxclvd')
%              [clvd,dblcpl]=mt_decomp(mt,'clvd')
%              [...,vec]=mt_decomp(...)
%
%    Description:
%     [ISO,DEV]=MT_DECOMP(MT,'ISO') returns the isotropic and deviatoric
%     components of the moment tensor(s) given in MT.  MT must be either a
%     scalar struct as output by FINDCMT/FINDCMTS, a Nx6 array, or a 3x3xN
%     array where N is the number of moment tensors in MT.  Note that
%     all GlobalCMT & USGS moment tensors are constrained to be purely
%     deviatoric so there will not be any non-zero values in ISO.  THE
%     OUTPUTS ARE DIAGONALIZED 3x3xN matrices!  See the last usage form
%     to get the eigenvectors which may be used with MT_UNDIAG to plot the
%     moment tensors (see the Examples section below too).
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
%     [DC1,DC2,DC3]=MT_DECOMP(MT,'3DC') returns 3 double couples for each
%     moment tensor in MT.  See Jost & Hermann (1989) for algorithm
%     details.
%
%     [CLVD1,CLVD2,CLVD3]=MT_DECOMP(MT,'3CLVD') returns 3 compensated
%     linear-vector dipoles (clvd) for each moment tensor in MT.  See Jost
%     & Hermann (1989) for algorithm details.
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
%     [clvd,best,vec]=mt_decomp(cmts,'clvd');
%     clvd=mt_undiag(clvd,vec);
%     best=mt_undiag(best,vec);
%     figure;
%     plotmt(1:10,2,cmts)
%     hold on;
%     plotmt(1:10,1,clvd)
%     plotmt(1:10,0,best)
%     axis equal tight off
%
%    See also: MT_DIAG, MT_UNDIAG, PLOTMT

%     Version History:
%        June 10, 2011 - initial version
%        Mar. 20, 2013 - 3dc & 3clvd options
%        Mar. 25, 2013 - update for mt_check/mt_change
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 25, 2013 at 13:50 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check tensor format (force 3x3xN)
error(mt_check(mt));
mt=mt_change('g',mt);
n=size(mt,3);

% check option string
valid={'iso' 'majmin' 'majmid' 'midmin' ...
    'maxdc' 'maxclvd' 'clvd' '3dc' '3clvd'};
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
        % this is a bit strange b/c the double couples
        % are calculated from the iso+dev components
        [dc1,dc2,dc3]=deal(zeros(3,3,n));
        for i=1:n
            % dc1 (p & b)
            dc1(1,1,i)=(mt(1,1,i)-mt(2,2,i))/3;
            dc1(2,2,i)=-dc1(1,1,i);
            
            % dc2 (b & t)
            dc2(2,2,i)=(mt(2,2,i)-mt(3,3,i))/3;
            dc2(3,3,i)=-dc2(2,2,i);
            
            % dc3 (t & p)
            dc3(3,3,i)=(mt(3,3,i)-mt(1,1,i))/3;
            dc3(1,1,i)=-dc3(3,3,i);
        end
        varargout={dc1 dc2 dc3 vec};
    case '3clvd' % Jost & Hermann 1989
        % this is a bit strange b/c the clvd components
        % are calculated from the iso+dev components
        [clvd1,clvd2,clvd3]=deal(zeros(3,3,n));
        for i=1:n
            % clvd1 p
            clvd1(1,1,i)=2*mt(1,1,i)/3;
            clvd1(2,2,i)=-mt(1,1,i)/3;
            clvd1(3,3,i)=-mt(1,1,i)/3;
            
            % clvd2 b
            clvd2(1,1,i)=-mt(2,2,i)/3;
            clvd2(2,2,i)=2*mt(2,2,i)/3;
            clvd2(3,3,i)=-mt(2,2,i)/3;
            
            % clvd3 t
            clvd3(1,1,i)=-mt(3,3,i)/3;
            clvd3(2,2,i)=-mt(3,3,i)/3;
            clvd3(3,3,i)=2*mt(3,3,i)/3;
        end
        varargout={clvd1 clvd2 clvd3 vec};
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
