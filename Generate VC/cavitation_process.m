function [PLC,RWD,Pressures,Kmax,Kini,voltot,embNb,ncNb,Embolized,Penetrated] = ...
    cavitation_process(sim,jt,stepsize,itmax,inEmb, AvgBool)
%AvgBool specifies whether to average all VC iterations
%   = 1 No averaging
%   = 0 Averaging occurs
voltot = sum([sim.Conduits.Length].*[sim.Conduits.Diameter].^2.*pi./4);
K = zeros(jt,itmax+1);
embNb = zeros(jt,itmax);
ncNb = zeros(jt,itmax);
Embolized = cell(jt,itmax);
Penetrated = cell(jt,itmax);
RWD = zeros(jt,itmax);
j=1;

sim.RepairAll;
[~,~,Kmax,~,~] = compute_hydraulics(sim,0.1,0);

while j<=jt

    [K(j,1),~,Embolized{j,1},Penetrated{j,1}] =...
        sim.AirSeed(0.1+stepsize,1:length(sim.Clusters),inEmb); 
%     figure;
%     sim.plotNet
%     text(20,52,'\DeltaP=0 MPa','FontSize',14,'FontWeight','bold')
%     ylim([0 55])
%     hold on
% %     axis tight
%     set(gca,'nextplot','replacechildren','visible','off')
%     f = getframe;
%     [im,map] = rgb2ind(f.cdata,256,'nodither');
%     im(1,1,1,20) = 0;

    past=0.1+stepsize;
    i=2;

    while i<=itmax && K(j,i-1)/K(j,1)>0.1
        past=past+stepsize;
 
        [Kcurrent,~,Embolized{j,i},Penetrated{j,i}]=sim.AirSeed(past,0,0);
%         sim.plotNet
%         text(20,52,strcat('\DeltaP=',num2str(past-0.1) ,' MPa'),...
%             'FontSize',14,'FontWeight','bold')
%         ylim([0 55])
%         f = getframe;
%         im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
        
        if Kcurrent~=-1
            K(j,i)=Kcurrent;
            embNb(j,i)=sum([sim.Conduits.Embolized]);
            ncNb(j,i)=sum([sim.Conduits.NonConducting]);
            embcond = sim.Conduits([sim.Conduits.Embolized]);
            RWD(j,i) = sum([embcond.Length].*[embcond.Diameter].^2.*pi./4);
        else
            K(j,i)=K(j,i-1);
            embNb(j,i)=sum([sim.Conduits.Embolized]);
            ncNb(j,i)=sum([sim.Conduits.NonConducting]);
            embcond = sim.Conduits([sim.Conduits.Embolized]);
            RWD(j,i) = sum([embcond.Length].*[embcond.Diameter].^2.*pi./4);
        end
        i=i+1;
        
    end
%     im = im(:,:,1,1:i-1);
    %At failure, set all subsequent values to the one reached at failure
    %pressure
    embNb(j,find(embNb(j,:)~=0,1,'last'):end) =...
        embNb(j,find(embNb(j,:)~=0,1,'last'));
    ncNb(j,find(ncNb(j,:)~=0,1,'last'):end) =...
        ncNb(j,find(ncNb(j,:)~=0,1,'last'));
    RWD(j,find(RWD(j,:)~=0,1,'last'):end) =...
        RWD(j,find(RWD(j,:)~=0,1,'last'));
    j=j+1;
end
%Reset simulation for future use
sim.RepairAll
%Save GIF
% imwrite(im,map,'DancingPeaks.gif','DelayTime',0.3,'LoopCount',inf) %g443800
%Process data before outputting
if size(K,1)>1 && ~AvgBool
    Kini=mean(K(:,1));
    PLC=1-mean(K)./Kini;
    embNb=mean(embNb)./length(sim.Conduits);
    ncNb=mean(ncNb)./length(sim.Conduits);
    RWD = mean(RWD);
elseif size(K,1)>1 && AvgBool
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