function selectedSNP = finemap_qc(G,maaf_thr,corr_thresh)

N=size(G,1);
M=size(G,2);

hst=hist(G);
hst2=hst(end,:);
freqs=sqrt(hst2/N); %maf from G

ff=find(freqs<maaf_thr | freqs>(1-maaf_thr)); % to eliminate lower maaf



GG= G-mean(G); % centralization of each SNP

H1=(1/N)*sum(GG.*GG);

gcov=cov(GG);

rr=gcov(:)'./(N*sqrt(kron(H1,H1))); % r matrix

rrr=reshape(rr,M,M)*N;

r=rrr/rrr(1,1); % interestingly r(1,1) is always 1.0001

aa1=[ff];
SNPcount=0;
%prunning with corr and maf threshold

for i=1:M
    
    
    %if sum(i==aa1)>0 || length(find(i==freqs))>0
    if sum(i==aa1)>0
        
    else
        a=find(abs(r(i,:))>corr_thresh);
        at=a(find(a~=i));
        aa1=[aa1 at];
        SNPcount=SNPcount+1;
        selectedSNP(SNPcount)=i;
    end
    
    
    
end


