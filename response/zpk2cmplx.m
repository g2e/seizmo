function [h]=zpk2cmplx(f,z,p,k,wpow)
%ZPK2CMPLX    Zeros-Poles-Constant to Complex Spectra
%
%    Usage:    h=zpk2cmplx(f,z,p,k)
%              h=zpk2cmplx(f,z,p,k,wpow)
%
%    Description:
%     H=ZPK2CMPLX(F,Z,P,K) calculates the complex spectra, H, of a
%     zero-pole-constant system at frequencies F.  F is in Hz.  Zeros are
%     given in Z, poles are in P and the constant is given in K.  Z and P
%     are expected to follow the omega convention (see the Notes section
%     for more details).  K should be a real scalar.
%
%     H=ZPK2CMPLX(F,Z,P,K,WPOW) changes the power of the angular frequency
%     divisor in the calculation of the complex response.  By default, WPOW
%     is 0 (division by 1).  If Z,P,K give the displacement response then
%     setting WPOW to 1 will give the velocity response.  Setting WPOW to 2
%     gives the acceleration response in the same case.
%
%    Notes:
%     - Poles and Zeros must be in units of angular frequency.  So a pole
%       at 2Hz on the real axis will be 2*pi*2Hz = 12.566+0i.  ZPK2CMPLX
%       checks that all complex poles and zeros have a conjugate pair.
%     - The formula for getting the complex spectra is:
%
%                 (iW-Z )*(iW-Z )*...*(iW-Z  )
%                      1       2           NZ
%        H(W)=K * ____________________________
%                 (iW-P )*(iW-P )*...*(iW-P  )
%                      1       2           NP
%
%       Where H(W) holds the complex spectra, W represents angular
%       frequencies corresponding to 2*PI*F, and i is the imaginary number.
%
%    Examples:
%     % Read in a SAC PoleZero file, get the complex spectra
%     % and plot the amplitude and phase response in velocity:
%     [z,p,k]=readsacpz('some/pz/file')
%     f=logspace(-3,1,1000);
%     h=zpk2cmplx(f,z,p,k,1);
%     figure; loglog(f,abs(h));
%     figure; plot(f,angle(a)); set(gca,'xscale','log');
%
%    See also: ZPK2AP, READSACPZ, READSACPZ_RDSEED, GETSACPZ, REMOVESACPZ,
%              APPLYSACPZ

%     Version History:
%        Oct. 19, 2009 - initial version
%        Oct. 22, 2009 - fixed NaN for 0Hz in vel/acc response
%        Feb. 11, 2011 - mass nargchk fix
%        Feb.  3, 2012 - doc update
%        Mar. 10, 2014 - update See also section, bugfix: prod needs dim
%                        parameter for single pole or zero, bugfix: zero
%                        out response where freq == 0 not h(1)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 10, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(4,5,nargin));

if(~isreal(f))
    error('seizmo:zpk2cmplx:badHz','F must be a real array!');
end
z=cplxpair(z(:));
p=cplxpair(p(:));
if(~isreal(k) || ~isscalar(k))
    error('seizmo:zpk2cmplx:badConstant','K must be a scalar real!');
end
if(nargin==4 || isempty(wpow)); wpow=0; end
if(~isscalar(wpow) || wpow~=fix(wpow))
    error('seizmo:zpk2cmplx:badPower','WPOW must be a scalar integer!');
end
nf=numel(f); nz=numel(z); np=numel(p);
w=complex(0,2*pi*f(:).');
h=k.*prod(w(ones(nz,1),:)-z(:,ones(nf,1)),1) ...
    ./prod(w(ones(np,1),:)-p(:,ones(nf,1)),1);
h=h./w.^wpow;
if(wpow>0); h(f==0)=0; end

end
