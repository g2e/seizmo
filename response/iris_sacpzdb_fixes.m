% Fixes to IRIS SAC PoleZero Database:
%
% THESE FIXES ARE IMPLEMENTED IN FIX_BAD_SACPZ!
%
% 1. Fixed 30 bad zeros for GE network:
%    Reason: Conjugate pair sign error.
%    Fix: z(4)=conj(z(4));
%    from:
%                     -15.15                         
%                     -176.6                         
%                     -463.1 +                 430.5i
%                     -463.1 +                 430.5i
%                          0                         
%                          0                         
%                          0                         
%    to:
%                     -15.15                         
%                     -176.6                         
%                     -463.1 +                 430.5i
%                     -463.1 -                 430.5i
%                          0                         
%                          0                         
%                          0                         
%
% 2. Fixed 70 bad poles for PN network:
%    Reason: Unpaired complex value turns out to be a typo based on
%            comparison with BK.FARB..BHZ (drop imag part of last pole).
%            Note that I left the other response descrepancies.
%    Fix: p(4)=real(p(4));
%    from:
%                    -0.1486 +                0.1486i
%                    -0.1486 -                0.1486i
%                    -414.69                         
%                   -999.027 -               999.027i
%    to:
%                    -0.1486 +                0.1486i
%                    -0.1486 -                0.1486i
%                    -414.69                         
%                   -999.027                         
%
% 3. Fixed 1 bad pole for CZ network:
%    Reason: Conjugate pair imaginary portion typo. Fixed based on
%            comparison to orthogonal channel (BHE) which is very similar.
%    Fix: p(end)=p(end)-1i
%    from:
%                    -0.0055                         
%                    -0.0173 +                0.0179i
%                    -0.0173 -                0.0179i
%                     -7.212 +                17.415i
%                     -7.212 -                17.415i
%                    -17.415 +                 7.212i
%                    -17.415 -                 7.212i
%                   -23.4459                         
%                   -37.0379                         
%                   -45.4483                         
%                    -8.7307 +               42.4144i
%                    -8.7307 -               41.4144i
%    to:
%                    -0.0055                         
%                    -0.0173 +                0.0179i
%                    -0.0173 -                0.0179i
%                     -7.212 +                17.415i
%                     -7.212 -                17.415i
%                    -17.415 +                 7.212i
%                    -17.415 -                 7.212i
%                   -23.4459                         
%                   -37.0379                         
%                   -45.4483                         
%                    -8.7307 +               42.4144i
%                    -8.7307 -               42.4144i
%
% 4. Fixed 6 bad poles for network AF:
%    Reason: Conjugate pair sign error.
%    Fix: p(2)=conj(p(2))
%    from:
%                    -0.9211 +                  0.94i
%                    -0.9211 +                  0.94i
%    to:
%                    -0.9211 +                  0.94i
%                    -0.9211 -                  0.94i
%
% 5. Fixed 6 bad poles for network AZ:
%    Reason: Conjugate pair sign error.
%    Fix: p(3)=conj(p(3))
%    from:
%                       -981 +                  1009i
%                       -981 -                  1009i
%                      -3290 -                  1263i
%                      -3290 -                  1263i
%    to:
%                       -981 +                  1009i
%                       -981 -                  1009i
%                      -3290 +                  1263i
%                      -3290 -                  1263i
%


% NEED INFO:
% 6. Fixed 4 bad poles for network BK:
%    Reason: Typo in imaginary part of conjugate pair. Which is right?
%    Fix: 
%    from:
%                   -10239.7 +               2725.02i
%                   -10239.7 -               2725.02i
%                   -9512.74 +                 11470i
%                   -9512.74 -                 11470i
%                   -454.526                         
%                   -396.553                         
%    !              -105.139 -               390.123i
%    !              -105.139 +               392.217i
%                   -15.4294                         
%                    -0.0371 +                0.0368i
%                    -0.0371 -                0.0368i
%    to:
%
% 7. Fixed 4 bad poles for network BK:
%    Reason: Typo in imaginary part of conjugate pair. Which is right?
%    Fix: 432 is common
%    from:
%                   -10239.7 +               2725.02i
%                   -10239.7 -               2725.02i
%                   -9512.74 +                 11470i
%                   -9512.74 -                 11470i
%                   -454.526                         
%                   -396.155                         
%    !              -105.023 -               389.432i
%    !              -105.023 +               392.573i
%                   -15.4441                         
%                    -0.0371 +                0.0368i
%                    -0.0371 -                0.0368i
%    to:
%
% 8. Fixed 4 bad poles for network BK:
%    Reason: Typo in imaginary part of conjugate pair. Which is right?
%    Fix: 814 is common
%    from:
%                   -10239.7 +               2725.02i
%                   -10239.7 -               2725.02i
%                   -9512.74 +                 11470i
%                   -9512.74 -                 11470i
%                   -454.526                         
%                   -396.951                         
%    !              -105.254 -               390.814i
%    !              -105.254 +               391.861i
%                   -15.4147                         
%                    -0.0371 +                0.0368i
%                    -0.0371 -                0.0368i
%    to:
%
% 9. Fixed 2 bad poles for network BK:
%    Reason: Typo in imaginary part of conjugate pair. The first is
%            probably right as the second one matches a real value from an
%            adjacent polezero file (I think they did an off-by-one-line
%            thing here - this calls for a complete check on RAMR responses
%            to sort out additional issues for sure).
%    Fix: p(8)=conj(p(7))
%    from:
%                   -10239.7 +               2725.02i
%                   -10239.7 -               2725.02i
%                   -9512.74 +                 11470i
%                   -9512.74 -                 11470i
%                   -454.526                         
%                   -399.485                         
%    !              -103.798 -               391.065i
%    !              -103.798 +               397.349i
%                   -15.5132                         
%                    -0.0371 +                0.0367i
%                    -0.0371 -                0.0367i
%    to:
%                   -10239.7 +               2725.02i
%                   -10239.7 -               2725.02i
%                   -9512.74 +                 11470i
%                   -9512.74 -                 11470i
%                   -454.526                         
%                   -399.485                         
%                   -103.798 -               391.065i
%                   -103.798 +               391.065i
%                   -15.5132                         
%                    -0.0371 +                0.0367i
%                    -0.0371 -                0.0367i
%



% NOT DONE (MAINLY B/C NOT BH*, LH*, OR HH*):
%%%%%%%%%%%%%
% AV is typo city!!! (several different ones)
%   1 p(1)~=conj(p(2))
%   2 p(5)=conj(p(5))
%   3 p(5)=conj(p(5))
%   4 p(5)=conj(p(5))
%   5 p(1)~=conj(p(2))
%   6 p(5)~=conj(p(6))
%   7 p(5)~=conj(p(6))
%   8 p(1)~=conj(p(2))
%   9 p(1)~=conj(p(2))
%  10 p(5)~=conj(p(6))
%  11 p(1)~=conj(p(2))
% BK is HUGE.....
% FA p(1)~=conj(p(2))    (HN[ZEN])
%                     -220.8 +                   265i
%                     -220.8 -                   260i
% H2 z(4)~=conj(z(5))    (BDH, HDH, LDH)
%                      -0.04                         
%                    -0.4708 +                 2.788i
%                    -0.4708 -                 2.788i
%                    -0.0785 +                1.2542i
%                    -0.0785 -               12.5418i
%                          0                         
% SG p(4)~=conj(p(5))    (LGZ)
%                     0.5638 +                2.3355i
%                     0.5638 -                2.3355i
%                    -2.2193                         
%                     0.0992 +                 1.411i
%                     0.0992 -                1.4111i
%                    -0.8904                         
%                     0.4618 +                 0.995i
%                     0.4618 -                 0.995i
%                     0.0872 +                0.7266i
%                     0.0872 -                0.7266i
%                     0.3662                         
%                     0.1415                         
% TA z(2)~=conj(z(3))    (BDE, LDE, VDE, UDE)
%                    -4.4005                         
%                    -0.0026 +                0.0024i
%                          0 -                0.0024i
% XA p(6) has no conj, p(4)~=conj(p(8))   (EL[ZEN])
%                     -401.5 +                 481.8i
%                     -401.5 -                 481.8i
%                     -139.8 +                 612.6i
%                     -391.7 +                 491.2i
%                     -566.1 +                 272.6i
%                     -628.3 +                0.0016i
%                     -566.1 -                 272.6i
%                     -391.8 -                 491.2i
%                     -139.8 -                 612.6i
% XF p(6) has no conj, p(4)~=conj(p(8))   (EH[ZEN])
%                     -8.796 +                 8.974i
%                     -8.796 -                 8.974i
%                     -139.8 +                 612.6i
%                     -391.7 +                 491.2i
%                     -566.1 +                 272.6i
%                     -628.3 +                0.0016i
%                     -566.1 -                 272.6i
%                     -391.8 -                 491.2i
%                     -139.8 -                 612.6i
% XL p(1)*100    (SPE)
%                    -0.1053 +                 7.864i
%                    -10.533 -                 7.864i

%     Version History:
%        Dec.  2, 2009 - initial version
%        Feb.  4, 2010 - fixed some typos
%        May  28, 2010 - AZ fix good, BK fixes pending, note fix_bad_sacpz
%        Feb.  3, 2012 - doc update
%        Mar.  6, 2014 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2014 at 02:25 GMT
