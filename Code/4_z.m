m_float=4866; %kg
r_float=1; %m
m_vib=2433; %kg
r_vib=0.5; %m
rou=1025; %kg/m**3
k=80000; %弹簧的k
zita=10000; %第一问为常数
zitarange=[0:1000:100000] ;%网格搜索步长
z20=0.0; %初值
l0=0.5; %m 弹簧原长
h_vib=0.5; %振子高度
h_float=3; %圆柱高
h_cone=0.8; %圆锥高
l_cen=0.0; %const中计算空心圆锥质心
g=9.8;
z_water=2.000;
m_att=1091.099; %垂荡附加质量
f=1760; %入射波震荡频率
omiga=1.9806; %垂荡激励振幅
gama=528.5018; %兴波阻尼系数

kn=250000; %扭转弹簧的刚度
zitan=1000;%旋转阻尼系数
I_f=7142.493; %纵摇附加转动惯量
L=2140; %纵摇激励力矩振幅
alpha=8890.7; %静水恢复力矩系数
gaman=1655.909; %纵摇兴波阻尼系数
mi=0;
period=10;

%一些给粒子群的小尝试
t=zeros(length(zitarange),length(zitarange),1000);%这个1000是估计值


time=2*pi/omiga*40;
dertax0=m_vib*g/k;
z10=h_vib-dertax0;
powerall=zeros(length(zitarange),length(zitarange));%用于存储所有计算出的最终功率
x0=[0 0 0 0];
z0=[z10 0 0 0]; %初始值；
% 定义 x(1)=z1, x(2)=z1', x(3)=z2, x(4)=z2';
for j1=1:length(zitarange)
for j2=1:length(zitarange)
%计算数值解,每次更新zita（和幂）
zitan=zitarange(j1);
zita=zitarange(j2);
dx=@(t1,x)[x(2); (1/I1)*(-kn*(x(1)-x(3))-zitan*(x(2)-x(4))+m_vib*zc1*sin(x(1))*g); x(4);1/(I2+I_f)*(kn*(x(1)-x(3))+zitan*(x(2)-x(4))-m_float*zc2*sin(x(3))*g+L*cos(omiga*t1)-gaman*x(4)-alpha*x(3))];
dz=@(t2,z)[z(2); (-m_vib*g+(l0-z(1)+z(3))*k-zita*(z(2)-z(4)))/m_vib; z(4);(rou*g* 3.14 * r_float^2 * (z_water-z(3)+1/3*h_cone)-(m_float) * g-k*(l0-z(1)+z(3))+zita*(z(2)-z(4))+f * cos(omiga*t2)-gama*z(4))/(m_float+m_att)];
tspan=[0,time];%time为40个波浪周期

powerz=0.0;
powerc=0.0;
[t1,x]=ode15s(dx,tspan,x0);
[t2,z]=ode15s(dz,tspan,z0);
startz=1;
startc=1;
%找到共同的开始值
for startznum=1:length(t1)
    if t1(startznum)>10
        startz=startznum
        break
    end
end
for startcnum=1:length(t2)
    if t2(startcnum)>10
        startc=startcnum
        break
    end
end
%求纵摇
    for i=startz:length(t1)
        if t1(startz)+period<=t1(i)
            break
        end
        dertat=t1(i+1)-t1(i);
        powerz=powerz+zitan*(x(i,2)-x(i,4))^2*dertat;
    end
%求垂荡
    for i=startc:length(t2)
        if t2(startc)+period<=t2(i)
            break
        end
        dertat=t2(i+1)-t2(i);
        powerc=powerc+zita*(z(i,2)-z(i,4))^2*dertat;
    end
%纵摇垂荡功率求和
power=(powerz+powerc)/period;%这一轮的最终功率
powerall(j1,j2)=power;
end
end
figure
[xpower, ypower] = meshgrid(zitarange,zitarange);
mesh(xpower,ypower,powerall);
xlabel('旋转阻尼系数(N·m·s)','FontSize',size);
ylabel( '直线阻尼系数(N·s/m)','FontSize',size );
zlabel('功率(W)','FontSize',size);
T=table(zitarange.',powerall);
writetable(T,"E:\TJ\TJ2022.9\数模校赛\Q4-zpower.csv");
%writematrix(x,"E:\TJ\TJ2022.9\数模校赛\Q1result-2.csv");