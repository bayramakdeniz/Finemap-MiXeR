function [g3x] = g333_delt(sigma02,a,glmt,m,q,s,sigma_beta2,deltOpt,pi1,u,k)

 %g3x=-(log(sqrt(sigma_beta2)/delt) - (s.^2 + m.^2 )/(2*delt^2)+(s.^2 + m.^2 )/(2*sigma_beta2)+ log(q/pi1)- log((1-q)/(1-pi1)));

delt=deltOpt;
 
qq=1./(1+exp(-k*u));

g3x=-((log(sqrt(sigma_beta2)./delt) - (s.^2 + m.^2 )./(2*delt.^2)+(s.^2 + m.^2 )/(2*sigma_beta2)+ log(qq/pi1)-log((1-qq)/(1-pi1))).*(qq.*(1-qq)*k));
 
%g3x=((log(sqrt(sigma_beta2)/delt) - (s.^2 + m.^2 )/(2*delt^2)+(s.^2 + m.^2 )/(2*sigma_beta2)+ log(qq/pi1)-log((1-qq)/(1-pi1))).*(qq.*(1-qq)*k));
 
  
 
end