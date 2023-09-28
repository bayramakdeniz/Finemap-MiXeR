function [g2x] = g222_new(sigma02,a,glmt,m,q,s,sigma_beta2,delt,pi1,u,k,B_term)

% g2x=(-1/(4*sigma02))*sum(4*s.*a'.*a') -s.*((1-q).*m/delt^2+(q.*m/sigma_beta2) - 1./(s.*s));

qq=1./(1+exp(-k*u));

 %g2x=(-1/(4*sigma02))*sum(4*s.*B_term) -s.*((1-qq).*m/delt^2+(qq.*m/sigma_beta2) - 1./(s.*s));
 
 g2x=(-1/(4*sigma02))*4*s.*B_term-s.*((1-qq).*m/delt^2+(qq.*m/sigma_beta2) - 1./(s.*s));

end