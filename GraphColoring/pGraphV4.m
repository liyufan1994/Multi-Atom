rng(1)

% 0 indicates Peterson graph, 1 indicates Nauru graph
PetersonOrNauru=1;

% pi1 or pi2 (0 or 1)--pi1 is uniform target, pi2 depends on number of red
% vertices--see paper
pi1Orpi2=0;

% Test functions: # of Red vertexes or Total color used (0 or 1)
RedOrTotal=1;

if PetersonOrNauru==0
    % This represents Peterson Graph
    Edges=...
        [0 1 0 0 1 1 0 0 0 0
        1 0 1 0 0 0 1 0 0 0
        0 1 0 1 0 0 0 1 0 0
        0 0 1 0 1 0 0 0 1 0
        1 0 0 1 0 0 0 0 0 1
        1 0 0 0 0 0 0 1 1 0
        0 1 0 0 0 0 0 0 1 1
        0 0 1 0 0 1 0 0 0 1
        0 0 0 1 0 1 1 0 0 0
        0 0 0 0 1 0 1 1 0 0];
else 
    % This is Nauru Graph
    Edges=...
    [0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
     1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0
     0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
     0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
     1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0
     0 0 0 0 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 1 0
     0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0
     0 0 0 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0
     0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0 0 0
     0 1 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 1 0 0 0 0 0 0
     0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1
     0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0
     1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0
     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1
     0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0];
end
    
% These are a few simpler graph in case reader would like to test on them

% Edges=...
%     [0 1 1 0 0
%      1 0 1 0 1
%      1 1 0 1 0
%      0 0 1 0 1
%      0 1 0 1 0]; % 1 means there is edge, 0 means no edge

%  Edges=...
%      [0 1 1 1
%       1 0 1 0
%       1 1 0 1
%       1 0 1 0];
% 




Edges=logical(Edges);


nv=size(Edges,2);
nvarr=1:nv;


p=7; % Number of Colors
parr=1:p;


len=10000; % Total number of steps

maskR=diag(1:p)*ones(p,nv);

if pi1Orpi2==0
    pdf=@(x) 1;
else
    pdf=@(x) x^5;
end


% All the last color numbers
Num=len; % number of samples
SA=[0 1 2 3 4];

spind=ones(1, length(SA));

coloredcol=zeros(Num,nv+1,length(SA));

for jj=1:length(SA)
    
    nball=SA(jj);
    
    %samplesck=zeros(Num,nball);
    colored=zeros(Num,nv+1);
    
    for i=1:Num
        if mod(i,1000)==0
            i
        end
        if i==2438
            1==1;
        end
        
        valicki=randsample(nv,nball, false)';
        if any(any(Edges(valicki,valicki)))==1
            colored(i,:)=[zeros(1,nv) 1/nchoosek(nv,nball)];
            continue;
        end
        
        Yt=zeros(1,nv);
        Yt(valicki)=p;
        prob=1;
        for j=1:nv
            if Yt(j)==0 % Not update color 4 nodes
                ponemask=true(1,p);
                nbcl=Yt(Edges(j,:)); % neighbor colors
                nbcl(nbcl==0)=[]; % remove 0's
                ponemask(nbcl)=false;
                ponemask(end)=false;
                avc=parr(ponemask);
                
                % avc=setdiff(1:(p-1), Yt(Edges(j,:)));
                if isempty(avc)
                    Yt=zeros(1,nv);
                    break;
                end
                prob=prob/length(avc);
                if length(avc)==1
                    Yt(j)=avc;
                else
                    Yt(j)=randsample(avc,1);
                end
            end
        end
        colored(i,:)=[Yt prob/nchoosek(nv,nball)];
        
    end
    coloredcol(:,:,jj)=colored;
end

% After we have the sample, we may define pistar for each SA element
pistar=zeros(1,length(SA));
for i=1:length(SA)
    probi=coloredcol(:,:,i);
    % Average of the probability of last column without invalid samples
    probi=mean(probi(probi(:,1)~=0,end));
    pistar(i)=pdf(SA(i))/probi;
end

% Store if on A or on E
sprec=zeros(len,1);
sprec(1)=-1; % Starting from A

yrec=zeros(len,nv);
atomrec=zeros(len,1);
atomrec(1)=1;


for tt=2:len
    if mod(tt,1000)==0
        tt
    end
    
    
    if sprec(tt-1)==1 % currently on E
        
        
        % Take a step on E using Gibbs sampler
        Yt=yrec(tt-1,:);
        for ii=1:nv
            
            % Find valid color
            % valcol=setdiff(parr,Yt(logical(Edges(ii,:))));
            
            ponemask=true(1,p);
            ponemask((Yt(Edges(ii,:))))=false;
            
            valcol=parr(ponemask);
            
            % Sample uniformly from valid color
            if length(valcol)==1
                col=valcol;
            else
                col=randsample(valcol,1);
            end
            
            Yprop=Yt;
            Yprop(ii)=col;
            
            coin=unifrnd(0,1);
            
            accpt=(pdf(nnz(Yprop==p)))/(pdf(nnz(Yt==p))); % Proposal is uniform
            
            if coin<min(1, accpt)
                Yt=Yprop;
            end
            
        end
        
        if ismember(nnz(Yt==p),SA)
            at=nnz(Yt==p);
            prob=(1./nchoosek(nv,at));
            Zt=Yt;
            Zt(Zt~=p)=0;
            for kk=1:length(Yt)
                if Zt(kk)~=p
                    ponemask=true(1,p);
                    nbcl=Zt(Edges(kk,:)); % neighbor colors
                    nbcl(nbcl==0)=[]; % remove 0's
                    ponemask(nbcl)=false;
                    ponemask(end)=false;
                    avc=parr(ponemask);
                    
                    % avc=setdiff(1:(p-1),Zt(Edges(kk,:)));
                    prob=prob/length(avc);
                    Zt(kk)=Yt(kk);
                end
            end
            mask=1:length(SA);
            atsanum=mask(SA==at);
            accpt=(pistar(atsanum)*prob)/pdf(at);
            
            coin=unifrnd(0,1);
            if coin<accpt
                
                atomrec(tt)=atsanum;
                sprec(tt)=-1;    
            else
                yrec(tt,:)=Yt;
                sprec(tt)=1;
            end
            
            
        else
            % Not in SA, just not transition, record current step
            yrec(tt,:)=Yt;
            sprec(tt)=1;
        end
        
        
    else % Currently on A
        
        %atsanum=atomrec(tt-1); % Comment out for we want to shuffle
        % shuffle
        atsanum=randsample(1:length(SA),1,true, pistar./sum(pistar));
        
        prop=coloredcol(spind(atsanum),:,atsanum); % Take the next sample
        spind(atsanum)=spind(atsanum)+1; % Move index downward one bit
        
        % Check if the sample is valid
        if prop(:,1)==0
            
            % Not valid then record the result and finish this iteration
            atomrec(tt)=atomrec(tt-1);
            sprec(tt)=-1;
            
        else
            
            % Valid, then check the if accept the proposal
            accpt=pdf(SA(atsanum))/(pistar(atsanum)*prop(end));
            
            coin=unifrnd(0,1);
            if coin<accpt
                yrec(tt,:)=prop(1:end-1);
                sprec(tt)=1;
            else
                atomrec(tt)=atomrec(tt-1);
                sprec(tt)=-1;
            end
            
        end
        
    end
    
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
        
        while (  (sprec(tt)~=-1) || ( ~ismembertol(atomrec(tt), SA) ) ) &&(tt<length(yrec)) 
            
            
            if sprec(tt)~=-1
                Ntnow=Ntnow+1;
                
                if RedOrTotal==0
                    Stnow=Stnow+nnz(yrec(jj,:)==p);
                else
                    Stnow=Stnow+length(unique(yrec(tt,:)));
                end
            end
            tt=tt+1;
            
        end
        St=[St;Stnow];
        Nt=[Nt;Ntnow];
    end
    tt=tt+1;
    
end
SS=0;
num=0;
for jj=1:size(yrec,1)
    if yrec(jj,1)~=0
        if RedOrTotal==0
            SS=SS+nnz(yrec(jj,:)==p);
        else
            SS=SS+length(unique(yrec(jj,:)));
        end
        num=num+1;
    end
end

% Mean estimate
u=SS/num;

% Error
sigmart=sqrt((sum((St-u.*Nt).^2))/(length(Nt)^2*mean(Nt)^2));

% Tour length
tourlength=length(Nt);

