
rng(2)
% Target
pdf = @(x) (0.5*normpdf(x,-3,0.5)+ 0.5*normpdf(x,3,0.5));

% Atom Stationary
pistar=@(x) 1.*(pdf(x))^(1);

% # of Partition cells for unit interval
M=5;

% Within-space transition for A
% We basically just use uniform of radius d around current position
d=5*M;
transproppdf=@(x,y) (abs(x-y)<=(d*(1/M)))/(2*d+1);
transproprnd=@(x) x+(1/M)*randi([-d,+d]);

% Within-space transition for E
proppdf = @(x,y) normpdf(x,y,0.5);
proprnd = @(x) normrnd(x,0.5);

% Total steps
len=100000;
yrec=zeros(len,1);

% This vector indicates whether the sample is in E or A
sprec=zeros(len,1);

% Starting point
yrec(1)=-3;
sprec(1)=1;

% All atoms to be collapsed
lb=2;
ub=4;
SA=[-4:1/M:2 2:1/M:4];


% Sample from SA proportional to pistar
% pistar for all points in SA
if length(SA)>1
    pistarSA=zeros(1,length(SA));
    for i=1:length(SA)
        pistarSA(i)=pistar(SA(i));
    end
    pistarSAnm=pistarSA/(sum(pistarSA)); % Normalized probability
    SASamp = randsample(SA,2*len,true,pistarSAnm); % These are all the presampled points;
elseif length(SA)==1
    SASamp=SA(1).*ones(1,2*len);
    pistarSA=pistar(SA(1));
    pistarSAnm=1;
end

% Cursor on SASamp
SAScur=1;

for tt=1:len
    if mod(tt,10000)==0
        tt
    end
    
    
    % V is current position after within-space transition; W is proposal in
    % perspective space
    
    if (sprec(tt)==-1) 
        
        if ismembertol(yrec(tt), SA)
            propSA=SASamp(SAScur);
            SAScur=SAScur+1;
        else
            propSA=yrec(tt);
        end
        Vprop=transproprnd(propSA);

        accept=min(1, pistar(Vprop)/pistar(propSA));

        
        coinmh=unifrnd(0,1);
        
        if coinmh<=accept
            V=Vprop;
        else
            V=propSA;
        end
        
        if ismembertol(V, SA)
            propSA=SASamp(SAScur);
            SAScur=SAScur+1;
        else
            propSA=V;
        end
        
        W=unifrnd(propSA,propSA+1/M);

    else
        
        Vprop=proprnd(yrec(tt));
        accept=min(1, pdf(Vprop)/pdf(yrec(tt)));
        coinmh=unifrnd(0,1);
        if coinmh<=accept
            V=Vprop;
        else
            V=yrec(tt);
        end
        
        % This formula finds the corresponding atom
        W=floor(V*M)/M;  
    end
    
    if (sprec(tt)==1) % Last step is in E

        a2VW=min(1, (pistar(W))/(pdf(V)));
    else
        % Locate corresponding atom in SA for W
        Vcor=floor(W*M)/M;
        a2VW=min(1, (pdf(W))/(pistar(Vcor)));
    end
    
    coin=unifrnd(0,1);
    if (coin <=a2VW)
        yrec(tt+1)=W; 
        sprec(tt+1)=-1*sprec(tt); % Accept trans-space jump
    else
        yrec(tt+1)=V;  
        sprec(tt+1)=sprec(tt); % Reject trans-space jump
    end


end



% With probability one, the chain does not hit 0 on E
clyrec=yrec(sprec==1);
imyrec=yrec(sprec==-1);
mean(clyrec)

St=[];
Nt=[];
tt=1;
while tt<length(yrec)

    if (sprec(tt)~=-1)
        Stnow=0;
        Ntnow=0;
        while (sprec(tt)~=-1 ||(~ismembertol(yrec(tt),SA)) ) &&(tt<length(yrec)) 
            if sprec(tt)~=-1
                Ntnow=Ntnow+1;
                Stnow=Stnow+yrec(tt);
            end
            tt=tt+1;
        end
        St=[St;Stnow];
        Nt=[Nt;Ntnow];
    end
    tt=tt+1;
    
end

% Mean Estimate
u=mean(clyrec);

% Variance Estimate
sigmart=sqrt((sum((St-u.*Nt).^2))/(length(Nt)^2*mean(Nt)^2))

fplot(pdf);hold on; 
histogram(clyrec,100, 'Normalization','pdf','FaceColor', [0 0 0]);
legend('Target Distribution', 'Samples')


