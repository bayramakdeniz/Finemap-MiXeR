function [theta errors  gradients LL non_scaled hh2 u alf smgd] = MyAdamNs2_rep_kf_delt_dum3(sigma02,a,glmt,sigma_beta2,delt,pi1,M,kf)

%%input

% sigma02:variance of error
% a: A matrix
% glmt: z-scores
% sigma_beta2: variance of causal beta values, 'discoverability'
% delta: a small number to rewrite dirac delta function as N(0,delta)
% pi1: prior
% M: numebr of SNP


%% output

%theta: final optimized parameters (beta_i sigma_i pi_i)
% errors: store the increments in each epoch
% upd: store the estimated parameters in each epoch
% grad: gradients in each epoch
% gradients: final gradients when the function is optimized

%%

  %rng(3);
  
  % m=rand(1,M);
  % s=rand(1,M);
  % q=rand(1,M);
   
   m=-0.05*ones(1,M);
   m=0.01*ones(1,M);
   s=0.01*ones(1,M);
   q=0.01*ones(1,M);
   u=-10*ones(1,M);
   u=-8*u/10;
   
   %u([ 110   338   455])=[10 10 10];
   u=3*abs(glmt');
   u(find(abs(u)>15))=15;
   
   u=abs(glmt');
   
   A_term1=a*glmt;
   A_term2=-a*a';
   B_term=sum(a'.*a');
  
   %A_term2(find(abs(A_term2)<0.1*max(max(abs(A_term2)))))=0;
   %%
    %% pca 
%    A3=A_term2;
%    pca_thr=0.9999;
%    mnx=mean(A3);
%    data=A3-mnx;
%    covmatrix = cov(data);
%    [eigvec,eigval] = eig(covmatrix);
%    [sh,ind] = sort((diag(eigval)),'descend');
%    K1=cumsum(sh)/sum(sh);
%    K11=find(K1>pca_thr);
%    K=K11(1)
%    %K=100;
%    reduced = data*eigvec(:,ind(1:K));
%    C_term=reduced;
%    eigs=eigvec(:,ind(1:K))';
   
   %%
   
  
   non_scaled=zeros(1,M);
   
    
    theta=[m s u];
    
    g1x= -1*g111_new(sigma02,a,glmt,m,q,s,sigma_beta2,delt,pi1,u,kf,A_term1,A_term2);
    %g1x=g111_new_pca(sigma02,a,glmt,m,q,s,sigma_beta2,delt,pi1,u,kf,A_term1,C_term,eigs,mnx)

    %g1x= -1*g111(sigma02,a,glmt,m,q,s,sigma_beta2,delt,pi1,u,kf);
    g2x= -1*g222(sigma02,a,glmt,m,q,s,sigma_beta2,delt,pi1,u,kf);
    g3x= -1*g333(sigma02,a,glmt,m,q,s,sigma_beta2,delt,pi1,u,kf);
    
    gradients=[g1x g2x g3x];
 
 
 
    mn = zeros(size(gradients));


    vr = zeros(size(gradients));

     alphaB = 1e-3;
     alpha = 1e-3; % from paper
    
    iteration = 1;

    
    inc=1;
    count=0;
    upd=[];
    differr=1;
    increment=1;
    hh2=[];
    
    errors=[];
    grad=[];
    LL=[];
    deltOpt=0.001*ones(1,M);
    %deltOpt=0.0005*ones(1,M);
   %while (sqrt(sum(gradients.^2))>0.01  &&  iteration<0.5*1*10^3 )
   sumgrad=1000;
   
   while ( iteration<10*10^3 &&   sumgrad>2.5*10^-3 )    
   %while (sqrt(sum(gradients.^2))>0.01  &&  iteration<15*10^3 ) % stopping criterion can be eased to have faster convergence
    count=count+1;
    thtold=1./(1+exp(-kf*u));
    
    alfa=alphaB;
    
 %   alfa=2*alphaB/sqrt(count/4);
    if count<800
        
    alfa=2*alphaB;   
    end  
    
    if count<10
        
    alfa=3*alphaB;   
    end  
    
    if count>2000
    alfa=200*alphaB/(count);   
    end    
    
    %alfa=(2*alphaB./(count/8).^0.5);
    %alfa=1*10^-3*exp(-.00005*count);
    alfa=alphaB;
    alf(count)=alfa;
    
    m=theta(1:M);
    s=theta(M+1:2*M);
    u=theta(2*M+1:end);
    
    g1x= -1*g111_new_delt(sigma02,a,glmt,m,q,s,sigma_beta2,deltOpt,pi1,u,kf,A_term1,A_term2);
    
    %g1x= -1*g111_new_delt_sparse(sigma02,a,glmt,m,q,s,sigma_beta2,deltOpt,pi1,u,kf,A_term1,A_term2);
    %g1x=-1*g111_new_delt_pca(sigma02,a,glmt,m,q,s,sigma_beta2,deltOpt,pi1,u,kf,A_term1,A_term2,C_term,eigs,mnx);

    %g1x= -1*g111(sigma02,a,glmt,m,q,s,sigma_beta2,delt,pi1,u,kf);
    g2x= -1*g222_new_delt(sigma02,a,glmt,m,q,s,sigma_beta2,deltOpt,pi1,u,kf,B_term);
    g3x= -1*g333_delt(sigma02,a,glmt,m,q,s,sigma_beta2,deltOpt,pi1,u,kf);
    
    gradients=[g1x g2x g3x];
    beta1 = 0.9; % from paper
    
    
    updateprev=theta;
    
    
   deltOpt=mean(sqrt((theta(M+1:2*M).^2+theta(1:M).^2)));


    beta2 = 0.99;  % from paper
     

    epsilon = 1e-8; %from paper


  
%alfat=alpha*sqrt(1-beta2^iteration)/(1-beta1^iteration);


% first moment estimate update

mn = beta1 * mn + (1 - beta1) * gradients;
    
% second  moment estimate update

vr = beta2 * vr + (1 - beta2) * gradients.^2;
    
% first moment estimate correction

mhat = mn / (1 - beta1^iteration);
    
% second moment estimate correction

vhat = vr / (1 - beta2^iteration);
    
% update parameters
theta =updateprev-alfa.* mhat ./ (sqrt(vhat) + epsilon);

non_scaled=non_scaled-alfa.* mhat(2*M+1:end) ./ (sqrt(vhat(2*M+1:end)) + epsilon);

%hh2=theta(2*M+1:end);
%theta=updateprev-alfat*mn/(sqrt(vr)+epsilon);

% for q
% q should be between 0 and 1. Also gradient function will lead some
% unmeaningful results (logarithm of a negative number) if pi is not
% between 0 and 1 hence we need to re-arrange q values if this case is not
% satisfied. For now I forced them to stay on border if they exceed the
% limits

 hh=non_scaled;
 dc=alfa.* mhat(2*M+1:end) ./ (sqrt(vhat(2*M+1:end)) + epsilon);

 if length(find(hh>1))>0
  %dec=max(hh>1)-0.999999;
   %dec=max(hh(hh>1))-0.9999;
   dec=max(hh)-1;
   
   %dec(dec<0)=max(dec(dec>0))*fliplr(max(dec(dec<0))./dec(dec<0));
   %dec(dec<0)=min(dec(dec>0));
   %dec(find(dec<0))=min(hh(hh>1));
   %dec(find(dec<0))=0.5*dc(find(dec<0));
  %theta(2*M+1:end)=theta(2*M+1:end)-dec;
  hh2=hh-dec;
  hh2(find(hh2<0))=10^-5;
 end


%theta(2*M+find(theta(2*M+1:end)<0))=10^-4;

%theta(2*M+1:end)=1./(1+exp(-10*(theta(2*M+1:end)-0.5)));

%theta(2*M+find(theta(2*M+1:end)>1))=0.999999;

%theta(2*M+1:end)=rescale(theta(2*M+1:end),10^-4,1-10^-4);

%theta(2*M+1:end)=1./(1+exp(-theta(2*M+1:end)));

%theta(2*M+1:end)=tanh(theta(2*M+1:end));

% hh=theta(2*M+1:end);
% 
% if length(find(hh>1))>0 && length(find(hh<0))>0
%   
% theta(2*M+1:end)=rescale(theta(2*M+1:end),10^-4,1-10^-4);
% end
% 
% 
% 
% if  length(find(hh>1))>0 && length(find(hh<0))==0
%     
%     theta(2*M+1:end)=rescale(theta(2*M+1:end),min(hh)-10^-20,1-10^-4);
%     
% end
%    
% 
% if  length(find(hh>1))==0 && length(find(hh<0))>0
%     
%     theta(2*M+1:end)=rescale(theta(2*M+1:end),10^-4,max(hh));
%     
% end


%theta(2*M+1:end)=rescale(theta(2*M+1:end),10^-4,1-10^-4);


% sigma should not be negative and if this is the case then it forced to
% turn to the  previous value

mm=max(10^-7,min(theta(M+1:2*M)));

theta(M+1:2*M)=rescale(theta(M+1:2*M),mm,max(theta(M+1:2*M)));

%theta(M+find(theta(M+1:2*M)<0))=10^-4;



%theta(M+find(theta(M+1:2*M)<0))=10^-3;

% update value

inc=norm(theta-updateprev)/norm(theta);

% update iteration number
iteration = iteration + 1;



%sumgrad=sum(abs(gradients(2*M+1:end)))/M;
%sumgrad=sum(abs(mhat(2*M+1:end) ./ (sqrt(vhat(2*M+1:end)) + epsilon)));
%tht(count,:)=1./(1+exp(-kf*u));

thtnew=1./(1+exp(-kf*u));

if count>10
sumgrad=sum(abs(thtnew-thtold));
end
% if count>2
% sumgrad=sum(abs(tht(count,:)-tht(count-1,:)))/M;
% end


 LL=1; 


% % not so important just to record the path of theta
% if mod(count,10)==0 || count==1
%     if count==1
% %     upd(1,:)=theta;
% %     grad(1,:)=gradients;   
% %     LL(1) = Lfunc(sigma02,a,glmt,m,1./(1+exp(-u)),s,sigma_beta2,delt,pi1);
%     else
% %     upd(1+count/10,:)=theta;
% %     grad(1+count/10,:)=gradients;
%     LL(1+count/10) = Lfunc(sigma02,a,glmt,m,1./(1+exp(-u)),s,sigma_beta2,delt,pi1);
%     smgd(1+count/10)=sumgrad;
%     end
%     
% end
% %

 smgd(count)=sumgrad;
%% the rest does not effect the algorithm
errors(count)=inc;

if count>1
differr=abs(inc-errors(count-1))/errors(count-1);
end




%f(count)=sqrt(sum(gradients.^2))



   end

theta(2*M+1:end)=1./(1+exp(-kf*u));
hh2(hh2<0)=0;

iteration
end
