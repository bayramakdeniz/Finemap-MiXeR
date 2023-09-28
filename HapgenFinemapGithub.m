clear all

close all

 

%rng('default')

 

 

% A script to run Finemap Mixer by producing artificial SNP for a continious trait

% Just run it as it is, make sure that MyAdamNs2_rep_kf.m g111.m g222.m g333.m are available in your pwd

 


% outputs: 

 

%pp posterior p_i

%bt2 posterior beta_i

 

for tt=1:5 

tt; % number of causals 

       

for ii=1:50 % different trials

 

clearvars -except  data_creator ii tt  err_pi err_miss ground_beta ground_pi type1 type11 power1 power11 tim errpi2 varpower11 vartype11  power_beta power_beta1 snps corr_thresh AUC AUC11 AUCc AUC11c pow1 typ1 pow1c typ1c powROC powROCc typROC typROCc Logf pvn pear_corr snpN pear_corr_train ty pv_tr pv_ts ggain

Ncausal=tt;
totalSamples=10000; %  hapgen chr21 10K samples



N=8000; % number of samples

M=1000; % number of SNPs

h2_d=0.5*10^-2; % desired heritability

 

corr_thresh=.99; % eliminate SNPs corr with .99 or higher

maaf_thr=0.05;

 

 

regions=[1:M]+randi(30000);  % pich a random region

 

 

 

G11=PlinkRead_binary2(N,regions,'chr21'); % get G matrix from this region

G2=double(G11);

%G1=sparse(G2);

Gf=G2;

 

 

%% apply basic QC

corr_thresh=.99; % eliminate SNPs corr with .99 or higher
maaf_thr=0.05;
selectedSNP = finemap_qc(Gf,maaf_thr,corr_thresh);

 


%% correlation matrix calculation, r


G2=Gf(:,selectedSNP); 
G=G2; 

M=size(G,2);
N=size(G,1);


GG= G-mean(G); % centralization of each SNP
H1=(1/N)*sum(GG.*GG);
gcov=cov(GG);

% quick calculation of correlation matrix r

rr=gcov(:)'./(N*sqrt(kron(H1,H1)));
rrr=reshape(rr,M,M)*N;
r=rrr/rrr(1,1); 

 

%% choose causal SNPs randomly

beta=zeros(1,M); % beta vector
correlationThreshold=0.99; % causals are chosen for the SNPs with whose all corr values are below than this threshold

 

eligibleSNPs=find(sum(double(abs(r(1:M,:))>correlationThreshold))<2); 
causalSnps=randperm(length(eligibleSNPs),Ncausal);
beta(eligibleSNPs(causalSnps))=1; % assign causal SNPs

beta1=beta;  %record beta's into beta1
pi1=sum(beta)/M; % prior prob of being causal


%% simu linux part . (creating phenotype vector for a given heritability, etc)

 
a1=sqrt(h2_d/var(GG*beta')); % for adjusting  y
env=randn(N,1); %noise
env=env/sqrt(var(env));

 

b=sqrt(1-h2_d); % for adjusting noise
y=a1*GG*beta'+b*env;  %pheno vector

 

var(y);
h2=var(a1*GG*beta')/var(y); % actual heritability, almost same with the desired

 

%% GWAS analysis to get z-scores:

 

rj=corr(GG,y); % correlation btwn G and y

tj=rj*sqrt(N-2)./sqrt(1-rj.^2); % zscores


 

p_values2 = (2*normcdf(-abs(tj)));
p_values=p_values2;
pd = makedist('Uniform');

 

 

glmt=tj; %z scores
hst=hist(G);
hst2=hst(end,:);
freqs=sqrt(hst2/N); %maf from G
freqs=sum(G)/(2*N);
 

mybeta =std(y)'*rj.*(1./std(GG, 'omitnan'))';

sigma=sqrt(var(y)-2*mybeta.^2.*freqs.*(1-freqs));

%my_se=sigma./sqrt(2*freqs.*(1-freqs)*N);

my_se=mybeta./glmt;

%freqs=maf*ones(1,M) % expected G

Hmin=2*freqs.*(1-freqs);
sigma_beta2=h2/(sum(Hmin)*pi1); 

 
% calculate A matrix
 
aa=(sqrt(N*reshape(repelem(H1,M),M,M)).*r');

a=aa';

 

sigma02=1; % sigma_0^2, default values

 

 


 

%% RUN ADAM


kf=1; % value for reparametrization of p_i=1/(1+e^-(kf u))

adj=0.15;  % adjustment parameter for delta
delt=adj*sqrt(sigma_beta2);

%delt=5*10^-3;


tic


%[updates errors  gradients LL non_scaled hh2 u alf smgd gain tm1 tm2] = MyAdamNs2_rep_kf_delt_dum2pca3noparam(sigma02,a,glmt,sigma_beta2,delt,pi1,M,kf);


% [updates errors  gradients LL non_scaled hh2 u alf smgd gain tm1 tm2] = MyAdamNs2_rep_kf_delt_dum2pca3(sigma02,a,glmt,sigma_beta2,delt,pi1,M,kf);


%[updates errors  gradients LL non_scaled hh2 u alf smgd] = MyAdamNs2_rep_kf_delt_dum3(sigma02,a,glmt,sigma_beta2,delt,pi1,M,kf);
% 


%sigma_beta2=(var(glmt)-1)/sum(Hmin);


%[updates errors  gradients LL non_scaled hh2 u alf smgd] = MyAdamNs2_rep_kf_delt_dum3_noparam(sigma02,a,glmt,sigma_beta2,delt,pi1,M,kf);

[updates errors  gradients LL non_scaled hh2 u alf smgd] = FinemapMiXeRv0(a,glmt);


toc

tim(ii)=toc

%ggain(ii)=gain

 

 

 

bt2=updates(1:M); % posterior beta values

pp=updates(2*M+1:end); % posterior p_i values

 

% two different teststats

 

tststat=pp; %just p_i

tststat2=pp.*bt2; %product of p_i and mu_i

 

%% calculation AUC

 

thrAlf=0.01:0.01:0.99; % type 1

 

for mm=1:length(thrAlf) % for different type1 thresholds power and AUC will be evaluated

 alfaThr=thrAlf(mm);

 

 vv=1:M;

 dd=eligibleSNPs(causalSnps);

 vv(dd)=[];

 

 pp2=tststat(vv); % get tstats for non causals

 pp22=tststat2(vv);

 

 

 [ss cc]=sort(pp2,'ascend'); % sort null tstats

 [sss ccc]=sort(pp22,'ascend');

 

 

 pp3=pp2(cc(1:length(cc)-round(alfaThr*length(cc)))); 

 mythr=pp3(end); % find the threshold among non_causals to get desired alfaThr

 

 

%do the same for other teststat=p_i x mu_i

 

  pp33=pp22(ccc(1:length(ccc)-round(alfaThr*length(ccc))));

  mythr2=pp33(end);

  type1(ii)=length(find(tststat(beta1==0)>mythr))/(M-Ncausal);

  power1(ii)=length(find(tststat(beta1==1)>mythr))/(Ncausal);

 

  power1m(mm)=power1(ii); % power for desired type 1

  type1m(mm)=type1(ii); %  obtained type1, almost same with the desired

 

% same for other tsttat

 

  type1c(ii)=length(find(tststat2(beta1==0)>mythr2))/(M-Ncausal);

  power1c(ii)=length(find(tststat2(beta1==1)>mythr2))/(Ncausal);

  power1mc(mm)=power1c(ii);

  type1mc(mm)=type1c(ii);

 

 

end

 

pow1(ii,:)=power1m; % record power for ii th experiments

typ1(ii,:)=type1m;

 

pow1c(ii,:)=power1mc;
typ1c(ii,:)=type1mc;

 

AUC(ii)=trapz(power1m)/length(type1m) % calculate Area under ROC curve

 

AUCc(ii)=trapz(power1mc)/length(type1mc)

 

end

 

powROC(tt,:)=mean(pow1);

typROC(tt,:)=mean(typ1);

 

powROCc(tt,:)=mean(pow1c);

typROCc(tt,:)=mean(typ1c);

 

 

AUC11(tt)=mean(AUC) % AUC if test stat is prob

AUC11c(tt)=mean(AUCc); %AUC if test stat is beta x prob

 

end

%%



