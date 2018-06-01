function Conx = pickConx(potConx,Pc)
    
    Prob=rand(length(potConx),1);
    Conx=cell(1,100000);
    count=1;
    
    for k=1:length(potConx)
       if Prob(k)>=(1-Pc) 
           Conx{count}=potConx{k};
           count=count+1;
       end
    end
    Conx=Conx(1:count-1);
end