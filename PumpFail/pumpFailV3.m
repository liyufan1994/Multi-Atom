import NamedConst
sarr=NamedConst.sarr;
tarr=NamedConst.tarr;
alpha=NamedConst.alpha;
gam=NamedConst.gam;
sig=NamedConst.sig;


% The result in the manuscript is from using seeds:
% 8 9 10 11 12 13

rng(8) 

% Total steps
len=100000;

% Dimension
Dims=length(sarr)+1;

% If we use Gibbs step on E or Metropolis step
useGibbs=1;

% If we use Metropolis within Gibbs on A
usemetGibbs=1;

[M, atmproj]=MatmBuild();
Mout=1000*M;
isinSA=@(x) isinSAun(x,atmproj);
SASamp=randSACdf(3*len,atmproj);
isinSAi=@(x,i) isinSAiun(x,i,atmproj);
pistardiv=@(x,y) pistardivun(x,y, atmproj,M,Mout);


% Within-space transition for A
% We basically just use uniform of radius d around current position
% Note that we include all atoms lines on the radius


radius=NamedConst.Radi; % we assume this is integer, actual distance
d=0.*ceil(radius.*M); % refer to d as the proposal we truly use

% x,y should be both on grid, so abs(x-y).*M should be integer
transproppdf=@(x,y) all(round(abs(x-y).*M)<=d);
transproprnd=@(x) x+round(unifrnd(-d-0.5,d+0.5));


% Within-space transition for E
smet=0.00;
proppdf = @(x,y) mvnpdf(y,x,smet*eye(Dims)); % x is previous location
proprnd = @(x) mvnrnd(x,smet*eye(Dims));



yrec=zeros(len,Dims);

% This vector indicates whether the sample is in E or A
sprec=zeros(len,1);

% Starting point
yrec(1,:)=NamedConst.cent0;%ones(1,length(NamedConst.cent0));
sprec(1)=1;





% Cursor on SASamp
SAScur=1;

for tt=1:len
    if mod(tt,10000)==0
        tt
    end
    
    
    % V is current position after within-space transition; W is proposal in
    % perspective space
    
    if (sprec(tt)==-1) % Xt in A
        
        % Within-A transition
        
        if isinSA(yrec(tt,:)) 
            % The shuffling on SA
            propSA=SASamp(SAScur,:);
            SAScur=SAScur+1;
        else
            propSA=yrec(tt,:);
        end
        
        if usemetGibbs==0 % Use metropolis on A
            
            Vprop=transproprnd(propSA);
            
            accept=min(1, pistardiv(Vprop, propSA));
            
            coinmh=unifrnd(0,1);
            
            if coinmh<=accept
                V=Vprop;
            else
                V=propSA;
            end
            
        else % Use metropolis-within-Gibbs on A
            
            V=propSA; % Starting point
            for j=1:(Dims-1)
                
                
                
                % Proposal for this coordinate
                L=gamrnd(alpha+sarr(j),1/(tarr(j)+V(end)));
                Lprop=floor(L.*M(j))./M(j);
                
                Mj1=M(j);
                if ~isinSAi(Lprop,j)
                    Lprop=floor(L.*Mout(j))./Mout(j); 
                    Mj1=Mout(j);
                end
                
                Mj2=M(j);
                if ~isinSAi(V(j),j)
                    Mj2=Mout(j);
                end
                
                
                % The proposal is just gampdf integrated on each segment
                % of length 1/M; 
                Qprop=gamcdf(Lprop+1/Mj1,alpha+sarr(j), 1/(tarr(j)+V(end)))-gamcdf(Lprop,alpha+sarr(j), 1/(tarr(j)+V(end)));
                Qtm1=gamcdf(V(j)+1/Mj2,alpha+sarr(j), 1/(tarr(j)+V(end)))-gamcdf(V(j),alpha+sarr(j), 1/(tarr(j)+V(end)));
                
                LVprop=V;
                LVprop(j)=Lprop;
                
                accept=min(1, pistardiv(LVprop,V)*(Qtm1/Qprop));
             
                
                coinmh=unifrnd(0,1);
                
                if coinmh<=accept % If not accept, then we do not update the coordinate
                    V(j)=Lprop;
                end
                
            end
            
            % The last coordinate
            Vsum=sum(V(1:end-1));
            
            L=gamrnd(gam+(Dims-1)*alpha,1/(sig+Vsum));
            Lprop=floor(L.*M(end))./M(end);
            
            Mend1=M(Dims);
            if ~isinSAi(Lprop,Dims)
                Lprop=floor(L.*Mout(Dims))./Mout(Dims);
                Mend1=Mout(Dims);
            end
            
            Mend2=M(Dims);
            if ~isinSAi(V(end),Dims)
                Mend2=Mout(Dims);
            end
            
            Qprop=gamcdf(Lprop+1/Mend1,gam+(Dims-1)*alpha,1/(sig+Vsum))-gamcdf(Lprop,gam+(Dims-1)*alpha,1/(sig+Vsum));
            Qtm1=gamcdf(V(end)+1/Mend2,gam+(Dims-1)*alpha,1/(sig+Vsum))-gamcdf(V(end),gam+(Dims-1)*alpha,1/(sig+Vsum));
            
            LVprop=V;
            LVprop(end)=Lprop;
                

            accept=min(1, pistardiv(LVprop, V)*(Qtm1/Qprop));
            
            coinmh=unifrnd(0,1);
            
            if coinmh<=accept % If not accept, then we do not update the coordinate
                V(end)=Lprop;
            end

            
        end
        
        % From here is trans-space transition
        if isinSA(V)
            propSA=SASamp(SAScur,:);
            SAScur=SAScur+1;
        else
            propSA=V;
        end
        
        
        Mnow=zeros(1, Dims);
        for i=1:Dims
            bb=isinSAi(propSA(i),i);
            Mnow(i)=bb*M(i)+(1-bb)*Mout(i);
        end
        
        W=propSA+unifrnd(zeros(1, length(Mnow)), 1./Mnow);
        
    else % Xt in E
        
        if useGibbs==0
            Vprop=proprnd(yrec(tt,:));
            accept=min(1, pdfdiv(Vprop, yrec(tt,:) ));
            coinmh=unifrnd(0,1);
            if coinmh<=accept
                V=Vprop;
            else
                V=yrec(tt,:);
            end
        else
            V=yrec(tt,:);
            for j=1:(Dims-1)
                V(j)=gamrnd(alpha+sarr(j),1/(tarr(j)+V(end)));
            end
            V(end)=gamrnd(gam+(Dims-1)*alpha,1/(sig+sum(V(1:end-1))));
        end
        
        Vatom=floor(V.*M)./M;
        
        Mnow=zeros(1, Dims);
        for i=1:Dims
            bb=isinSAi(Vatom(i),i);
            Mnow(i)=bb*M(i)+(1-bb)*Mout(i);
        end
        
        % This formula finds the corresponding atom
        W=floor(V.*Mnow)./Mnow;  
    end
    
    if (sprec(tt)==1) % Xt in E, propose to go to A

        a2VW=min(1, pistardivpdf(W,V));
        
    else % Xt in A, propose to go to E
        
        Watom=floor(W.*M)./M;
        
        Mnow=zeros(1, Dims);
        for i=1:Dims
            bb=isinSAi(Watom(i),i);
            Mnow(i)=bb*M(i)+(1-bb)*Mout(i);
        end
        
        % Locate corresponding atom in SA for W
        Vcor=floor(W.*Mnow)./Mnow;

        a2VW=min(1, pdfdivpistar(W,Vcor));
        %a2VW=min(1, pdf(W)/pistar(Vcor));
    end
    
    coin=unifrnd(0,1);
    if (coin <=a2VW)
        yrec(tt+1,:)=W; 
        sprec(tt+1)=-1*sprec(tt); % Accept trans-space jump
    else
        yrec(tt+1,:)=V;  
        sprec(tt+1)=sprec(tt); % Reject trans-space jump
    end


end



% "Extract" samples on E and A
clyrec=yrec(sprec==1,:); % Samples on E
imyrec=yrec(sprec==-1,:); % Samples on A
mean(clyrec)

St=[];
Nt=[];
tt=1;
while tt<length(yrec) % Loop through the sample path

    if (sprec(tt)~=-1) % Encounter values not on A, start collect path
        Stnow=0;
        Ntnow=0;
        
        % The path stops when and only when the chain hit SA on A
        
        while (  (sprec(tt)~=-1) || ( ~isinSA(yrec(tt,:)) ) ) &&(tt<length(yrec)) 
            
            
            if sprec(tt)~=-1
                Ntnow=Ntnow+1;
                Stnow=Stnow+yrec(tt,:);
            end
            tt=tt+1;
        end
        St=[St;Stnow];
        Nt=[Nt;Ntnow];
    end
    tt=tt+1;
    
end
ii=11
Sti=St(:,ii);

% The mean estimate
u=mean(clyrec(:,ii));

% The error/variance estimate
sigmart=sqrt((sum((Sti-u.*Nt).^2))/(length(Nt)^2*mean(Nt)^2))



function m=pdf(x)
import NamedConst

sarr=NamedConst.sarr;
tarr=NamedConst.tarr;
alpha=NamedConst.alpha;
gam=NamedConst.gam;
sig=NamedConst.sig;

if min(x)<=0
    m=0;
else
    
    lamm=x(1:end-1);
    bet=x(end);
    m=bet^(gam-1)*exp(-bet*sig)*prod((tarr.*lamm).^sarr)*exp(-sum(lamm.*tarr))*(prod(lamm))^(alpha-1)*exp(-bet*sum(lamm))*(bet)^(length(lamm)*alpha);
end
end

function m=pdfdiv(x,y) 
import NamedConst
sarr=NamedConst.sarr;
tarr=NamedConst.tarr;
alpha=NamedConst.alpha;
gam=NamedConst.gam;
sig=NamedConst.sig;

% The target ratio needs to be evaluated via log
% pi(x)/pi(y)=exp(log(pi(x))-log(pi(y)))

% y is current position and should be always positive
if (min(y)<=0)
    error('Last step not positive');
end

% x is proposal so if it is 0 or negative we simply output 0
if (min(x)<=0)
    m=0;
else
    
    lammx=x(1:end-1);
    betx=x(end);
    lammy=y(1:end-1);
    bety=y(end);
    logmx=(gam-1)*log(betx)+(-betx*sig)+sum(log(lammx.*tarr).*sarr)+(-sum(lammx.*tarr))+(alpha-1)*(sum(log(lammx)))+(-betx*sum(lammx))+(length(lammx)*alpha)*log(betx);
    logmy=(gam-1)*log(bety)+(-bety*sig)+sum(log(lammy.*tarr).*sarr)+(-sum(lammy.*tarr))+(alpha-1)*(sum(log(lammy)))+(-bety*sum(lammy))+(length(lammy)*alpha)*log(bety);
    
    m=exp(logmx-logmy);
end

end

function m=pdfdivpistar(x,y) 
import NamedConst


m=pdfdiv(x,y);

end

function m=pistardivpdf(x,y)
import NamedConst



m=pdfdiv(x,y);

end

function m=pistardivun(x,y, atmproj,M,Mout)
import NamedConst

% Find the volume of block for x and y
Dims=length(x);

logxvol=0;
for i=1:Dims
    if ismembertol(x(i),atmproj{i})
        logxvol=logxvol-log(M(i));
    else
        logxvol=logxvol-log(Mout(i));
    end
end

logyvol=0;
for i=1:Dims
    if ismembertol(y(i),atmproj{i})
        logyvol=logyvol-log(M(i));
    else
        logyvol=logyvol-log(Mout(i));
    end
end


m=pdfdiv(x,y)*exp(logxvol-logyvol);

end

function m=isinSAun(x,atmproj)
    m=0;
    for i=1:length(atmproj)
        if ~ismembertol(x(i), atmproj{i})
            return     
        end
    end
    m=1;
end

function m=isinSAiun(x,i,atmproj)
    m=ismembertol(x,atmproj{i});
end

function SASamp=randSACdf(len, atmproj)
% We have input as a cell array where each cell is an array with all
% projections of the atom

Latm=length(atmproj);


% Find dim of atoms
Narr=ones(1, Latm);
for i=1:Latm
    Narr(i)=length(atmproj{i});
end
N=prod(Narr); 

eval=zeros(1,Latm);
a=cell(1, Latm);

% Sample from SA proportional to pistar
% pistar for all points in SA
if N>1
    pistarSA=zeros(1,N);
    for i=1:N
        if mod(i,100000)==0
            i
        end

        [a{:}]=ind2sub(Narr,i);
        
        for j=1:Latm
            eval(j)=atmproj{j}(a{j});
        end
        pistarSA(i)=pdf(eval);
    end
    pistarSAnm=pistarSA/(sum(pistarSA)); % Normalized probability
    SASampind = randsample(1:N,len,true,pistarSAnm); % These are all the presampled indexes;
    SASamp=zeros(length(SASampind),Latm);
    for i=1:length(SASampind)
        if mod(i,10000)==0
            i
        end
        [a{:}]=ind2sub(Narr, SASampind(i));
        for j=1:Latm
            eval(j)=atmproj{j}(a{j});
        end
        SASamp(i,:)=eval; % Map to actual point
    end
elseif N==1
    SASamp=repmat(cell2mat(atmproj),len,1);
end

end


function [M, atmproj]=MatmBuild()
% An elementary build: basically for each coordinate, choose number of
% partition of unit interval as M(i) such that 1/M(i)=std(i); and then
% incorporate the atom whose cell contains mean and then the atom above it
% and the atom below it

% M(i) needs not to be a integer or even bigger than 1. Its purpose is
% just specify length of a partition cell and we assume the grid starts
% from 0 always. 
import NamedConst

M=1./(NamedConst.Radi./1);

x=floor(NamedConst.cent0.*M)./M;
X=[x-1./M;x; x+1./M; x+2./M];
%X=[x-1./M;x;x+1./M];
atmproj=cell(1, size(X,2));
for i=1:size(X,2)
    atmproj{i}=(X(:,i))';
end


end






