function X = randsphere(m,n,r)
% This function returns an m by n array, X, in which 
% each of the m rows has the n Cartesian coordinates 
% of a random point uniformly-distributed over the 
% interior of an n-dimensional hypersphere with 
% radius r and center at the origin.  The function 
% 'randn' is initially used to generate m sets of n 
% random variables with independent multivariate 
% normal distribution, with mean 0 and variance 1.
% Then the incomplete gamma function, 'gammainc', 
% is used to map these points radially to fit in the 
% hypersphere of finite radius r with a uniform % spatial distribution.
% Roger Stafford - 12/23/05

% Generating uniform distribution in a sphere/ball
% ellieandrogerxyzzy@mindspring.com.invalid (Roger Stafford)
% 21 Apr, 2006 01:26:43
% 
%  In the first line of 'randsphere' we have
%
% X = randn(m,n)
%
%which produces m rows, each consisting of n independent, normally
%distributed random variables, mean 0, variance 1. For the purposes of
%this discussion, assume that m equals 1 so that we would be talking about
%a single set of n values.
%
%  By definition, the sum of the squares, s2, in
%
% s2 = sum(X^2)
%
%gives a value with the chi-squared distribution with n degrees of
%freedom. If you look at Steven Wolfram's MathWorld website at
%
% http://mathworld.wolfram.com/Chi-SquaredDistribution.html,
%
%you can see in its equation (6) that the cumulative probability
%distribution for such a chi-squared distribution is:
%
% Prob(s2 <= c^2) = incomplete_gamma(n/2,c^2/2) / gamma(n/2)
%
%However, matlab's incomplete gamma function, 'gammainc', has been
%regularized by division with gamma(n/2) and its arguments reversed in
%order from this, so we have
%
% Prob(s2 <= c^2) = gammainc(c^2/2,n/2)
%
%with matlab's 'gamminc' function. This is the probability that the vector
%X will fall in a hypersphere of radius less than or equal to c. Note that
%fortunately this depends only on the radius c and not on the direction of
%X.
%
%  On the other hand, the n-dimensional volume of a hypersphere of radius s
%<= r is proportional to s^n so it occupies the fraction, (s/r)^n, of the
%whole desired hypersphere in 'randsphere', r being its specified radius.
%
%  The strategy in 'randsphere' then is to map the original X value so that
%it remains along the same ray but with a radius altered from c to s, in
%such a way that the fraction of volume of the hypersphere with radius s of
%the whole hypersphere will be equal to the probability of having had any
%original X within the hypersphere of radius c. This would then provide
%the hypersphere with a uniform distribution. That is, the probability of
%mapping into the hypersphere of radius s would be proportional to its
%volume, which is what we need.
%
%  This leads to the equality
%
% gammainc(c^2/2,n/2) = (s/r)^n (Matlab's 'gammainc')
%
%or
%
% s = r*(gammainc(c^2/2,n/2))^(1/n)
%
%which provides the necessary mapping relationship between c (=sqrt(s2)) and s.
%
%  The unit vector along the ray is X/sqrt(s2), so we multiply
%
% X/sqrt(s2)*r*(gammainc(s2^2/2,n/2))^(1/n)
%
%to obtain the transformed vector that lies within the desired
%hypersphere. This (hopefully) explains the third and last line in the
%code for 'randsphere'.
%
%  Confession: As it happened, back in December when this problem was being
%considered I had never noticed the 'gammainc' function in matlab. While
%trying to produce a uniform distribution for a hypersphere, after some
%furious analysis, I arrived at the point where I realized the need for an
%incomplete form of the gamma function and wondered how in the world it
%could ever be computed. This led me to rereading the description of
%matlab's gamma function in my manual and lo and behold there was the
%lovely 'gammainc' function in the same section staring me in the face,
%just waiting to be used. My face was a little red at the time!

X = randn(m,n);
s2 = sum(X.^2,2);
X = X.*repmat(r*(gammainc(s2/2,n/2).^(1/n))./sqrt(s2),1,n);
 
