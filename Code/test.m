zita=10000; %第一问为常数
zitarange=[0:500:100000]; %网格搜索步长
zita=zitarange(3);
t=[0:length(zitarange)];
len=length(zitarange);
%powerall=[length(zitarange),length(mirange)];%用于存储所有计算出的最终功率
figure
x=1:0.1:10;
y=1:0.1:10;
[x, y] = meshgrid(x,y);
mesh(zitarange,zitarange,yyc+yyz);
n=1;
for m=1:length(t1)
    if abs(t1(m)-t2(n))<=0.1 & t1(m)>10
        startz=m;
        startc=n;
        break
    elseif t1(m)>t2(n)
        n=n+1;
    else
        m=m+1;
    end
end

[zitanmax,zitamax]=find (powerall == max(max(powerall)))
[zitanmax,zitamax]