% Store state of the chain in a N-by-N matrices (N=32 for standard
% application) where the 3rd coordinate denotes time
import IsingData
rng(1);
N=IsingData.N;
theta=IsingData.theta;
len=100000;
preDATA=[];

for init=1:1 % We will presample two runs starting from all 1 and all -1
yrec=zeros(N,N,len);
yrec(:,:,1)=(-1)^init.*ones(N,N); % Initial state  


% Pre-sample all coin tosses required
coinarr=unifrnd(0,1,[32,32,len]);

for tt=2:len
    
    if mod(tt,1000)==0
        tt
    end
    
    % The "heat bath" Gibbs sampler
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
            
        end
    end
    
    yrec(:,:,tt)=Yt;
end

% Find natural stats
tstarr=zeros(len,1);
numonesarr=zeros(len,1);
for tt=1:len
    
    if mod(tt,1000)==0
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
newDATA=[(1:length(numonesarr))' numonesarr tstarr];
preDATA=[preDATA; newDATA];
end

save('prerun.mat','preDATA');


