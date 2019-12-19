import NamedConst
sarr=NamedConst.sarr;
tarr=NamedConst.tarr;
alpha=NamedConst.alpha;
gam=NamedConst.gam;
sig=NamedConst.sig;

% Pump Failure

% Number of pumps
N=length(sarr);

% Total steps that Gibbs sampler will run
TT=1000000;

% Failure rate for each of the 10 pumps
lambda=ones(TT,N);

% Observed Failures


beta=ones(TT,1);
%Gibbs sampler

for i=2:TT
    
    if mod(i,10000)==0
        i
    end
    
    for j=1:N
        lambda(i,j)=gamrnd(alpha+sarr(j),1/(tarr(j)+beta(i-1)));
    end
    beta(i)=gamrnd(gam+N*alpha,1/(sig+sum(lambda(i,:))));
    
end

samples=[lambda beta];
M=10;
kk=1*M; % Radius of SA as in number of atoms, integer
SArad=kk/M; % Actual radius
cent0=1.2*ones(1,N+1);

count=0; % number of points in S
for i=1:length(beta)
    if norm(samples(i,:)-cent0,inf)<=SArad
        count=count+1;
    end
end
count/length(beta)
    

FF1=@(lamm, bet) gampdf(bet,gam,1/sig)*poisspdf(1,lamm)*gampdf(lamm,alpha,1/bet);
FF2=@(lamm, bet) bet^(gam-1)*exp(-bet*sig)*lamm*exp(-lamm)*(lamm)^(alpha-1)*exp(-bet*lamm)*(bet)^(alpha);

pdf2=@(lamm, bet) bet^(gam-1)*exp(-bet*sig)*prod(lamm.^sarr)*exp(-sum(lamm))*(prod(lamm))^(alpha-1)*exp(-bet*sum(lamm))*(bet)^(length(lamm)*alpha);
pdf=@(x) pdf2(x(1:end-1),x(end));

% function m=target(lam, bet,gam,sig,sarr,alpha)
%     w=gampdf(bet,gam,1/sig);
%     for i=1:length(lam)
%         w=w*poisspdf(sarr(i),lam(i))*gampdf(lam(i),alpha,1/bet);
%     end
%     m=w;
% end
% 
% function m=target2(lam, bet,gam,sig,sarr,alpha)
%     w=bet^(gam-1)*exp(-bet*sig);
%     for i=1:length(lam)
%         w=w*poisspdf(sarr(i),lam(i))*(lam(i))^(alpha-1)*exp(-bet*lam(i));
%     end
%     m=w;
% end
% 
% a=1;b=3;c=4
% 100*gampdf(a,b,c)/(a^(b-1)*exp(-a/c))



