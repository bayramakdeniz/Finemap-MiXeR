function [g1x] = g111_new(sigma02,a,glmt,m,q,s,sigma_beta2,delt,pi1,u,k,A_term1,A_term2)

%g1x=((1/sigma02)*sum((a.*glmt'-a.*sum((a'.*m)'))')-((1-q).*m/delt^2)-q.*m/sigma_beta2); 

qq=1./(1+exp(-k*u));

g1x=(((1/sigma02)*(A_term1+A_term2*m'))'-((1-qq).*m/delt^2)-qq.*m/sigma_beta2); 

end