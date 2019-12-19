% Store state of the chain in a N-by-N matrices (N=32 for standard
% application) where the 3rd coordinate denotes time
numseed=10;
uarr=zeros(1,numseed);
sigarr=zeros(1, numseed);
tourarr=zeros(1, numseed);
reaarr=zeros(1, numseed);

for seed=1:numseed
rng(seed) 
import IsingData

if seed==1
    run('IsingV1.m')
end
file=load('prerun.mat','preDATA');
preDATA=file.preDATA;


N=IsingData.N;
theta=IsingData.theta;
len=10000;
yrec=zeros(N,N,len);
yrec(:,:,1)=ones(N,N); % Initial state
% bm=randi([0 1],N,N);
% bm(bm==0)=-1;
% yrec(:,:,1)=bm;


cache=zeros(1,N^2+1); % The third indx is number of -1's
iscache=zeros(1,N^2+1); 

% Pre-sample all coin tosses required
coinarr=unifrnd(0,1,[32,32,len]);

% Pre-sample another set of coin tosses required (for trans-space transition)
coinarr2=unifrnd(0,1,[len,1]);

% This vector indicates whether the sample is in E or A
sprec=zeros(len,1);
sprec(1)=1; % Starting on E

% This vector stores atom steps
atomrec=zeros(len,1);

% A diagnostic of largest possible Matm
% Criteria: the smallest t for Matm from preDATA must be less than max of
% uniform
if 1==1
    judgearr=zeros(1,N^2+1);
    for Matm=0:N^2
        if mod(Matm,200)==0
            Matm
        end
        
            
        if isempty(preDATA(preDATA(:,2)==Matm,:))
            
            continue
        end
        
        judge=0;
        for jj=1:15
            
            numm1s=Matm;
            inds=datasample(1:N^2, numm1s, 'Replace',false);
            [I,J]=ind2sub([N,N],inds);
            yprop=ones(N);
            for i=1:numm1s
                yprop(I(i),J(i))=-1;
            end
            %
            % % 5. Compute natural statistics of yprop
            tstyprop=0;
            for x=1:N
                for y=1:N
                    
                    
                    if (x==1)
                        xup=N;
                    else
                        xup=x-1;
                    end
                    
                    if (x==32)
                        xdown=1;
                    else
                        xdown=x+1;
                    end
                    
                    if (y==1)
                        yleft=N;
                    else
                        yleft=y-1;
                    end
                    
                    if (y==32)
                        yright=1;
                    else
                        yright=y+1;
                    end
                    
                    s=yprop(x,yleft)+yprop(xup,y)+yprop(x,yright)+yprop(xdown,y);
                    
                    tstyprop=tstyprop+yprop(x,y)*(s);
                end
            end
            tstyprop=tstyprop*0.5;
            
            if tstyprop>=min(preDATA(preDATA(:,2)==Matm,3))
                judge=1;
                break;
            end
        end
        
        
        judgearr(Matm+1)=judge;
    end
end


% Atoms till which to include in SA
Matm=[(0:(length(judgearr)-1))' judgearr'];
Matm=Matm(Matm(:,2)==1);


% First step is to generate representatives for all these atoms
Matmrept=zeros(length(Matm),1);
for i=1:length(Matm)
    Matmrept(i)=min(preDATA(preDATA(:,2)==Matm(i),3));
end

% Second step is to generate random draws from these atoms

logpistar=zeros(1, length(Matm));
for j=1:length(Matm)
    
    numm1s=Matm(j);
    if numm1s==0

        logpistar(j)=theta*2048;
    else
        %pistar(numm1s+1)=exp(prod(log((N^2-numm1s+1):(N^2)))-prod(log(1:numm1s))+theta*cache(numm1s)-Base0);
        logpistar(j)=sum(log((N^2-numm1s+1):(N^2)))-sum(log(1:numm1s))+theta*Matmrept(j);
    end
    
end

pistaru=exp(logpistar-max(logpistar));
pistar=pistaru/(sum(pistaru));



for tt=2:len
    
    % Just print the progress
    if mod(tt,2000)==0
        tt
    end
    
    if sprec(tt-1)==1 % Currently on E
        % 1. First take a step to yprop on E by Gibbs sampler
        
        Yt=yrec(:,:,tt-1); % We will update each entry of this matrixt
        for x=1:N
            for y=1:N
                
                if (x==1)
                    xup=N;
                else
                    xup=x-1;
                end
                
                if (x==32)
                    xdown=1;
                else
                    xdown=x+1;
                end
                
                if (y==1)
                    yleft=N;
                else
                    yleft=y-1;
                end
                
                if (y==32)
                    yright=1;
                else
                    yright=y+1;
                end
                
                s=Yt(x,yleft)+Yt(xup,y)+Yt(x,yright)+Yt(xdown,y);
                
                Ppos=1/(1+exp(-2*theta*s));
                coin=coinarr(x,y,tt);
                if coin<Ppos
                    Yt(x,y)=1;
                else
                    Yt(x,y)=-1;
                end
                %Yt(x,y)=randsample([-1 1],1,true,[1-Ppos Ppos]);
            end
        end
        
        yprop=Yt;
        
        % 2. Take yprop and compute natural statistics
        tstyprop=0;
        for x=1:N
            for y=1:N
                
                
                if (x==1)
                    xup=N;
                else
                    xup=x-1;
                end
                
                if (x==32)
                    xdown=1;
                else
                    xdown=x+1;
                end
                
                if (y==1)
                    yleft=N;
                else
                    yleft=y-1;
                end
                
                if (y==32)
                    yright=1;
                else
                    yright=y+1;
                end
                
                s=yprop(x,yleft)+yprop(xup,y)+yprop(x,yright)+yprop(xdown,y);
                
                tstyprop=tstyprop+yprop(x,y)*(s);
            end
        end
        tstyprop=tstyprop*0.5;
        
        % 3. Find representative: defined as the grid with -1 most spread
        % out
        
        numm1s=nnz(yprop==-1); % Count num of -1's
        
        if ismembertol(numm1s, Matm)
            tstrep=Matmrept(Matm==numm1s);
        else
            tstrep=0;
        end
        
        
        % 5. Compute accpt ratio and decide if we make the transition
        a2VW=min(1, exp(theta*(tstrep-tstyprop)));
        if coinarr2(tt)<a2VW
            sprec(tt)=-1;
            atomrec(tt)=numm1s;
            % Transition to atom, simply leave this space in yrec as 0's.
        else
            sprec(tt)=1;
            yrec(:,:,tt)=yprop;
        end
        
        
    else % Currently on A
        % 1. We will not take a within-A step. So just propose from current step.
        
        % 2. Find representative: defined as the grid with -1 most spread
        % out
        
        numm1s=atomrec(tt-1);
        
        % Shuffling
        if ismembertol(numm1s, Matm)
            if length(Matm)==1
                numm1s=Matm;
            else
                numm1s=randsample(Matm,1,true,pistar);
            end
        end
        
        if ismembertol(numm1s, Matm)
            tstrep=Matmrept(Matm==numm1s);
        else
            tstrep=0;
        end
        
        
        % 4. Propose to go to yprop in E uniformly conditioning on number
        % of -1's
        inds=datasample(1:N^2, numm1s, 'Replace',false);
        [I,J]=ind2sub([N,N],inds);
        yprop=ones(N);
        for i=1:numm1s
            yprop(I(i),J(i))=-1;
        end  
      
        
        % 5. Compute natural statistics of yprop
        tstyprop=0;
        for x=1:N
            for y=1:N
                
                
                if (x==1)
                    xup=N;
                else
                    xup=x-1;
                end
                
                if (x==32)
                    xdown=1;
                else
                    xdown=x+1;
                end
                
                if (y==1)
                    yleft=N;
                else
                    yleft=y-1;
                end
                
                if (y==32)
                    yright=1;
                else
                    yright=y+1;
                end
                
                s=yprop(x,yleft)+yprop(xup,y)+yprop(x,yright)+yprop(xdown,y);
                
                tstyprop=tstyprop+yprop(x,y)*(s);
            end
        end
        tstyprop=tstyprop*0.5;
        
        % 6. Compute accpt ratio and decide if we make the transition
        a2VW=min(1, exp(theta*(tstyprop-tstrep)));
        if coinarr2(tt)<a2VW
            sprec(tt)=1;
            yrec(:,:,tt)=yprop;
        else
            sprec(tt)=-1;
            atomrec(tt)=atomrec(tt-1);
        end
        
    end
    
    
end

% Find natural stats
tstarr=zeros(len,1);
numonesarr=zeros(len,1);
for tt=1:len
    
    if mod(tt,2000)==0
        tt
    end
    
    tst=0;
    Yt=yrec(:,:,tt);
    numones=0;
    for x=1:N
        for y=1:N
            if Yt(x,y)==-1
                numones=numones+1;
            end
            
            if (x==1)
                xup=N;
            else
                xup=x-1;
            end
            
            if (x==32)
                xdown=1;
            else
                xdown=x+1;
            end
            
            if (y==1)
                yleft=N;
            else
                yleft=y-1;
            end
            
            if (y==32)
                yright=1;
            else
                yright=y+1;
            end
            
            s=Yt(x,yleft)+Yt(xup,y)+Yt(x,yright)+Yt(xdown,y);
            
            tst=tst+Yt(x,y)*(s);
        end
    end
    
    tstarr(tt)=tst*0.5;
    numonesarr(tt)=numones;
end

% Compute estimation error
St=[];
Nt=[];
tt=1;
while tt<length(yrec) % Loop through the sample path

    if (sprec(tt)~=-1) % Encounter values not on A, start collect path
        Stnow=0;
        Ntnow=0;
        
        % The path stops when and only when the chain hit SA on A
        
        while (  (sprec(tt)~=-1) || ( ~ismembertol(atomrec(tt), Matm) ) ) &&(tt<length(yrec)) 
            
            
            if sprec(tt)~=-1
                Ntnow=Ntnow+1;
                Stnow=Stnow+tstarr(tt);
            end
            tt=tt+1;
            
        end
        St=[St;Stnow];
        Nt=[Nt;Ntnow];
    end
    tt=tt+1;
    
end

u=mean(tstarr(tstarr~=0));
sigmart=sqrt((sum((St-u.*Nt).^2))/(length(Nt)^2*mean(Nt)^2));
uarr(seed)=u;
sigarr(seed)=sigmart;
tourarr(seed)=length(Nt);
reaarr(seed)=length(sprec(sprec==1));

end
stats=[uarr; sigarr; tourarr; reaarr];
% realtst=tstarr(tstarr~=0);
% mean(realtst)
% std(realtst)
% numonesarrind=[(1:length(numonesarr))' numonesarr tstarr reshape(yrec(1,1,:),length(yrec(1,1,:)),1,1)];
% A=numonesarrind(numonesarr==1022,:);
% B=sortrows(A,3);
% 
% numm1s=99;
% inds=datasample(1:N^2, numm1s, 'Replace',false);
% [I,J]=ind2sub([N,N],inds);
% yprop=ones(N);
% for i=1:numm1s
%     yprop(I(i),J(i))=-1;
% end  
% % 
% % % 5. Compute natural statistics of yprop
% tstyprop=0;
% for x=1:N
%     for y=1:N
% 
% 
%         if (x==1)
%             xup=N;
%         else
%             xup=x-1;
%         end
% 
%         if (x==32)
%             xdown=1;
%         else
%             xdown=x+1;
%         end
% 
%         if (y==1)
%             yleft=N;
%         else
%             yleft=y-1;
%         end
% 
%         if (y==32)
%             yright=1;
%         else
%             yright=y+1;
%         end
% 
%         s=yprop(x,yleft)+yprop(xup,y)+yprop(x,yright)+yprop(xdown,y);
% 
%         tstyprop=tstyprop+yprop(x,y)*(s);
%     end
% end
% tstyprop=tstyprop*0.5
% nnz(yprop(:,17:end)==-1)
% atomreccb=[(1:length(atomrec))' atomrec];
