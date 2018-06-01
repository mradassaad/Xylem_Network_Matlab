function [Emb_hist, Pen_hist, figs] = cav_process_hist(Dias,ASPs,mVC)
    Embolized = mVC.Embolized;
    Penetrated = mVC.Penetrated;
%     Dias = [obj.Conduits.Diameter];
    [Emb_Y,Emb_E] = discretize(Dias,5);
    Emb_values = (Emb_E(1:end-1) + Emb_E(2:end))./2;
    Emb_valNb(1) = sum(Emb_Y == 1);
    Emb_valNb(2) = sum(Emb_Y == 2);
    Emb_valNb(3) = sum(Emb_Y == 3);
    Emb_valNb(4) = sum(Emb_Y == 4);
    Emb_valNb(5) = sum(Emb_Y == 5);
    
%     ASPs = [obj.ICConnections.ASP];
    [Pen_Y,Pen_E] = discretize(ASPs,5);
    Pen_values = (Pen_E(1:end-1) + Pen_E(2:end))./2;
    Pen_valNb(1) = sum(Pen_Y == 1);
    Pen_valNb(2) = sum(Pen_Y == 2);
    Pen_valNb(3) = sum(Pen_Y == 3);
    Pen_valNb(4) = sum(Pen_Y == 4);
    Pen_valNb(5) = sum(Pen_Y == 5);
    
    Emb_hist = zeros(5,size(Embolized,2));
    Pen_hist = zeros(5,size(Embolized,2));
    for i = 1 :size(Embolized,1)
        Emb_temp = zeros(5,size(Embolized,2));
        Pen_temp = zeros(5,size(Embolized,2));
        for j = 1:size(Embolized,2)
            for l = 1: length(Embolized{i,j})
                D = Embolized{i,j}(l);
                binNb = find(Emb_E(1:end-1)<=D & Emb_E(2:end)>D);
                Emb_temp(binNb,j) = Emb_temp(binNb,j) + 1;
                
                ASP = Penetrated{i,j}(l);
                binNb = find(Pen_E(1:end-1)<=ASP & Pen_E(2:end)>ASP);
                Pen_temp(binNb,j) = Pen_temp(binNb,j) + 1;
            end
        end
        Emb_hist = Emb_hist + Emb_temp;
        Pen_hist = Pen_hist + Pen_temp;
    end
    Emb_hist = cumsum(Emb_hist./size(Embolized,1),2);
    Pen_hist = cumsum(Pen_hist./size(Embolized,1),2);
    
    Emb_hist(Emb_valNb ~= 0,:) = Emb_hist(Emb_valNb ~= 0,:)./...
        Emb_valNb(Emb_valNb ~= 0)';
    Pen_hist(Pen_valNb ~= 0,:) = Pen_hist(Pen_valNb ~= 0,:)./...
        Pen_valNb(Emb_valNb ~= 0)';
    
    lim = find(mVC.PLC>0.95,1,'first') - 1;
    figs(1) = figure('Name','Dc',...
        'NumberTitle','off','Color','w','Units','inches',...
        'Position',[0 0 6 2]);

    ax1 = axes;
    ax1.FontSize = 9;
    yyaxis(ax1,'right')
  
    bar(ax1,mVC.Pressures(1:lim),Emb_hist(:,(1:lim))')
    ax1.YLimMode = 'manual';
    ax1.YLim=[0 0.6];
    ylabel('% embolized')
    barLgdStr1=strcat('D_c=',num2str(Emb_values'.*1e6),'\mum');
    
    yyaxis left
    plot(ax1,mVC.Pressures(1:lim),mVC.PLC(1:lim),'b-o')
    ax1.YLimMode = 'manual';
    ylim([0 1])
    xlim([0 5])
    lgd1=legend(ax1,char('VC',barLgdStr1),'Location','northwest');
    lgd1.FontSize=7;
    xlabel('Pressure (MPa)')
    ylabel('PLC')
    
    figs(2) = figure('Name','ASP',...
        'NumberTitle','off','Color','w','Units','inches',...
        'Position',[0 0 6 2]);
    ax2 = axes;
    ax2.FontSize = 9;
    yyaxis right
    bar(mVC.Pressures(1:lim),Pen_hist(:,(1:lim))')
    ax2.YLimMode = 'manual';
    ax2.YLim=[0 0.1];
    ylabel('% penetrated')
    barLgdStr2=strcat('ASP=',num2str(Pen_values'),'MPa');
    
    yyaxis left
    plot(mVC.Pressures(1:lim),mVC.PLC(1:lim),'b-o')
    ax2.YLimMode = 'manual';
    ylim([0 1])
    xlim([0 5])
    lgd2=legend(char('VC',barLgdStr2),'Location','northwest');
    lgd2.FontSize=7;
    xlabel('Pressure (MPa)')
    ylabel('PLC')
end