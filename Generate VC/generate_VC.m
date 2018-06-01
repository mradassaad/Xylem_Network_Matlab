function [PLC,Pressures,Kmax,Kini,embNb,ncNb] = ...
    generate_VC(sim,jt,stepsize,itmax,inEmb, ICLogical)
K=zeros(jt,itmax+1);
embNb=zeros(jt,itmax);
ncNb=zeros(jt,itmax);
j=1;
% count=1;
% min=itmax;
sim.RepairAll;
[~,~,Kmax,~,~] = compute_hydraulics(sim,0.1,0);

% while j<=jt && count <=100
while j<=jt
%     disp(num2str(j))
    [K(j,1),~]=sim.AirSeed(0.1+stepsize,1:length(sim.Clusters),inEmb); 
    past=0.1+stepsize;
    i=2;
%     flag=0;
%     while i<=itmax && K(j,i-1)/K(j,1)>0.05 && flag<5
    while i<=itmax && K(j,i-1)/K(j,1)>0.1
        past=past+stepsize;
 
        [Kcurrent,~]=sim.AirSeed(past,0,0);
        if Kcurrent~=-1
            K(j,i)=Kcurrent;
            embNb(j,i)=sum([sim.Conduits.Embolized]);
            ncNb(j,i)=sum([sim.Conduits.NonConducting]);
        else
            K(j,i)=K(j,i-1);
            embNb(j,i)=sum([sim.Conduits.Embolized]);
            ncNb(j,i)=sum([sim.Conduits.NonConducting]);
        end
        i=i+1;

%         if K(j,i-1)/K(j,1)<0.05
%             flag=flag+1;
%         end
    end
    embNb(j,find(embNb(j,:)~=0,1,'last'):end)=embNb(j,find(embNb(j,:)~=0,1,'last'));
    ncNb(j,find(ncNb(j,:)~=0,1,'last'):end)=ncNb(j,find(ncNb(j,:)~=0,1,'last'));
    j=j+1;
end
% embNb=embNb(find(sum(transpose(K))),:);%#ok<FNDSB>
% ncNb=ncNb(find(sum(transpose(K))),:);%#ok<FNDSB>
% K=K(find(sum(transpose(K))),:); %#ok<FNDSB>
sim.RepairAll
if size(K,1)>1 && ~ICLogical
    Kini=mean(K(:,1));
    PLC=1-mean(K)./Kini;
    embNb=mean(embNb)./length(sim.Conduits);
    ncNb=mean(ncNb)./length(sim.Conduits);
elseif size(K,1)>1 && ICLogical
    Kini=K(:,1);
    PLC=1-K./Kini;
    embNb=embNb./length(sim.Conduits);
    ncNb=ncNb./length(sim.Conduits);
else
    Kini=mean(K(:,1));
    PLC=1-K./Kini;
    embNb=embNb./length(sim.Conduits);
    ncNb=ncNb./length(sim.Conduits);
end
Pressures=stepsize:stepsize:(stepsize)*itmax;

PLC=PLC(:,1:end-1);
end