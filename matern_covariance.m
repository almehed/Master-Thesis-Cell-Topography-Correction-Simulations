function r = matern_covariance(h,sigma,kappa,nu)
% matern_covariance(dist,sigma,kappa,nu) computes the matern covariance
% function C(h) = (sigma^2/(2^(nu-1)*Gamma(nu))) * (kappa*h)^nu *
% K_v(kappa*h) at the distances given in the vector or matrix h.
% sigma^2 is the variance, nu decides the smoothness of C(h) at h=0, and
% kappa decides the practical correlation length r through
% r=sqrt(8*nu)/kappa.
%
% David Bolin (david.bolin@chalmers.se) 2018.

  if nargin < 4; error('Too few input parameters'); end
  if nu<=0; error('nu must be postive'); end
  if kappa<=0; error('kappa must be postive'); end
  if sigma<=0; error('sigma must be postive'); end

  rho = sqrt(8*nu)/kappa;
  hpos = (h>0);
  r = zeros(size(h));
  r(~hpos) = sigma^2;
  [hunique,~,J] = unique(h(hpos));

  B = log(besselk(nu,kappa*hunique));
  B = exp(2*log(sigma)-gammaln(nu)-(nu-1)*log(2) + nu.*log(kappa*hunique) + B);

  if any(isinf(B))
    B(isinf(B)) = sigma^2*exp(-2*(hunique(isinf(B))/rho).^2);
  end
  r(hpos) = B(J);
