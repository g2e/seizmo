function [n,d] = rrat(x,tol)
%RRAT    Relative rational approximation
%
%    Usage:    [n,d]=rrat(x,tol)
%              s=rrat(x,tol)
%
%    Description:
%     [N,D]=RRAT(X,TOL) returns the numerator and denominator matrices
%     that express the elements in X as a fraction of two small integers.
%     The integer fraction will be within TOL*X of X.  TOL is optional and
%     is by default 1e-6*norm(X(:),1).  This differs from RAT in that TOL
%     actually is relative to X (note that the RAT documentation indicates
%     that it is but in actuality it is not) rather than an absolute
%     tolerance value.
%
%     RRAT(X,TOL) or S=RRAT(X,TOL) returns the continued fraction expansion
%     as a string.
%
%    Notes:
%
%    Examples:
%     % this shows the benefit of relative vs absolute tolerance
%     a=1/300+rand/1e8;
%     rat(a)
%     rrat(a)
%     rat(a,1e-2)
%     rrat(a,1e-2)
%
%    See also: RAT, FORMAT, RATS

%     Version History:
%        Aug. 16, 2010 - gplv3 version derived from GNU Octave's RAT
%
%     Written by Paul Kienzle (RAT in GNU Octave)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 16, 2010 at 11:00 GMT

% todo:

% Copyright (C) 2001, 2007, 2008, 2009 Paul Kienzle
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

  error(nargchk(1,2,nargin));
  error(nargchk(0,2,nargout));

  y = x(:);

  % Replace Inf with 0 while calculating ratios.
  y(isinf(y)) = 0;

  % default norm
  if (nargin < 2)
    tol = 1e-6 * norm(y,1);
  end

  % First step in the approximation is the integer portion

  % First element in the continued fraction.
  n = round(y);
  d = ones(size(y));
  frac = y-n;
  lastn = ones(size(y));
  lastd = zeros(size(y));

  nsz = numel (y);
  steps = zeros(nsz,0);

  % Grab new factors until all continued fractions converge.
  while (1)
    % Determine which fractions have not yet converged.
    % edit by Garrett Euler - makes code match described behavior
    % (tol is in relative terms, not absolute terms)
    idx = find(abs((y-n./d)./y) >= tol);
    if (isempty(idx))
      if (isempty (steps))
        steps = NaN .* ones (nsz, 1);
      end
      break;
    end

    % Grab the next step in the continued fraction.
    flip = 1./frac(idx);
    % Next element in the continued fraction.
    step = round(flip);

    if (nargout < 2)
      tsteps = NaN .* ones (nsz, 1);
      tsteps (idx) = step;
      steps = [steps, tsteps];
    end

    frac(idx) = flip-step;

    % Update the numerator/denominator.
    nextn = n;
    nextd = d;
    n(idx) = n(idx).*step + lastn(idx);
    d(idx) = d(idx).*step + lastd(idx);
    lastn = nextn;
    lastd = nextd;
  end

  if (nargout == 2)
    % Move the minus sign to the top.
    n = n.*sign(d);
    d = abs(d);

    % Return the same shape as you receive.
    n = reshape(n, size(x));
    d = reshape(d, size(x));

    % Use 1/0 for Inf.
    n(isinf(x)) = sign(x(isinf(x)));
    d(isinf(x)) = 0;

    % Reshape the output.
    n = reshape (n, size (x));
    d = reshape (d, size (x));
  else
    n = '';
    nsteps = size(steps, 2);
    for i = 1: nsz
      s = [int2str(y(i)),' '];
      j = 1;

      while (true)
        step = steps(i, j);
        j=j+1;
        if (isnan (step))
          break;
        end
        if (j > nsteps || isnan (steps(i, j)))
          if (step < 0)
            s = [s(1:end-1), ' + 1/(', int2str(step), ')'];
          else
            s = [s(1:end-1), ' + 1/', int2str(step)];
          end
          break;
        else
          s = [s(1:end-1), ' + 1/(', int2str(step), ')'];
        end
      end
      s = [s, repmat(')', 1, j-2)];
      n_nc = size(n,2);
      s_nc = size(s,2);
      if (n_nc > s_nc)
        s(:,s_nc+1:n_nc) = ' ';
      elseif (s_nc > n_nc)
        n(:,n_nc+1:s_nc) = ' ';
      end
      n = cat (1, n, s);
    end
  end

end
