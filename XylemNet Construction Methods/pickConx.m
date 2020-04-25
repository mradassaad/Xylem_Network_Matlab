function Conx = pickConx(potConx_rad, potConx_tan, Pe_rad, Pe_tan, radDist)
    
    Pe_rad_rad = (radDist*Pe_rad(1) + (1- radDist)*Pe_rad(2));
    Pe_tan_rad = (radDist*Pe_tan(1) + (1- radDist)*Pe_tan(2));
    Prob_rad = rand(length(potConx_rad),1);
    Prob_tan = rand(length(potConx_tan),1);
    Conx=cell(1,100000);
    count=1;
    
    for k=1:length(potConx_rad)
       if Prob_rad(k)>=(1-mean(Pe_rad_rad(potConx_rad{k}(:,2)))) 
           Conx{count}=potConx_rad{k};
           count=count+1;
       end
    end
    
    for k=1:length(potConx_rad)
       if Prob_tan(k)>=(1-mean(Pe_tan_rad(potConx_tan{k}(:,2)))) 
           Conx{count}=potConx_tan{k};
           count=count+1;
       end
    end
    Conx=Conx(1:count-1);
end