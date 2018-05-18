function resizeConduits(obj,Dc,Dc_cv,endWallScheme,varargin)
Dstd=Dc_cv*Dc; %m
Dm=log(Dc^2/sqrt(Dstd^2+Dc^2));
Ds=sqrt(log(Dstd^2/(Dc^2)+1));
Dcs=lognrnd(Dm,Ds,1,length(obj.Conduits));
Dcs=sort(Dcs);
tempLen=[obj.Conduits.Length];
[~,ILen]=sort(tempLen);
fn = DcMapping(length(obj.Conduits),...
    floor(length(obj.Conduits)/20));
for k=1:length(obj.Conduits)
    obj.Conduits(k).updateDiameter(...
        Dcs(fn(ILen==k)).*(rand(1).*0.4+0.8));
end

resizePores(obj,endWallScheme,varargin{:})
end