function rho = spharm(t,p,n,m)
theta = p+pi/2;
phi = t;
Pmn = legendre(n,cos(theta(:)));
Pmn = reshape(Pmn(m+1,:),size(theta));
rho = sqrt((2*n+1)/(4*pi)*Gamma(n-m+1)/Gamma(n+m+1))*Pmn.*exp(i*m*phi);


