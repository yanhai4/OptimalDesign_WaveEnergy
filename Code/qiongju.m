m_float=4866; %kg
r_float=1; %m
m_vib=2433; %kg
r_vib=0.5; %m
rou=1025; %kg/m**3
k=80000; %弹簧的k
zita=10000; %第一问为常数
zitarange=[0:500:100000] ;%网格搜索步长
z20=0.0; %初值
l0=0.5; %m 弹簧原长
h_vib=0.5; %振子高度
h_float=3; %圆柱高
h_cone=0.8; %圆锥高
l_cen=0.0; %const中计算空心圆锥质心
g=9.8;
z_water=2.000;


m_att=1165.992; %垂荡附加质量
f=4890; %入射波震荡频率
omiga=2.2143; %垂荡激励振幅
gama=167.8395; %兴波阻尼系数
timestart=350;
timelength=10;
dertax0=m_vib*g/k;
z10=h_vib-dertax0;
powerall=[1:length(zitarange)];%用于存储所有计算出的最终功率
x0=[z10 0 0 0]; %初始值；
% 定义 x(1)=z1, x(2)=z1', x(3)=z2, x(4)=z2';
for j=1:length(zitarange)
%计算数值解,每次更新zita（和幂）
zita=zitarange(j);
dx=@(t,x)[x(2); (-m_vib*g+(l0-x(1)+x(3))*k-zita*(x(2)-x(4)))/m_vib; x(4);(rou*g* 3.14 * r_float^2 * (z_water-x(3)+1/3*h_cone)-(m_float) * g-k*(l0-x(1)+x(3))+zita*(x(2)-x(4))+f * cos(omiga*t)-gama*x(4))/(m_float+m_att)];
%dx=@(t,x)[x(2); (-m_vib*g-(x(1)-x(3))*k-10000*sqrt(abs(x(2)-x(4)))*(x(2)-x(4)))/m_vib; x(4);(rou*g* 3.14 * r_float^2 * (z_water-x(3)+1/3*h_cone)-(m_float) * g+k*(x(1)-x(3))+10000*sqrt(abs(x(2)-x(4)))*(x(2)-x(4))+f * cos(omiga*t)-gama*x(4))/(m_float+m_att)];
tspan=[0,360];%time为40个波浪周期
%dertat=[1:length(t)];
%power=[1:length(t)-1];
power=0.0;
[t,x]=ode15s(dx,tspan,x0);
    for i=length(t)-100:length(t)-1
        dertat=t(i+1)-t(i);
        power=power+zita*(x(i,2)-x(i,4))^2*dertat;
    end
timespend=t(length(t)-1)-t(length(t)-100);
power=power/timespend;%这一轮的最终功率
powerall(j)=power;
%[t,z1,z2]=dsolve('-m_vib * g - (z1-z2) * k - zita * (Dz1-Dz2) - m_vib * D2z1','rou * g * 3.14 * r_float ^ 2 * (z_water-z2+1/3*h_cone) - (m_att+m_float) * g + (z1-z2) * k + zita * (Dz1-Dz2) + f * cos(omiga*t)-gama*Dz2-m_float*D2z2','z1(0)=0,z2(0)=0,Dz1(0)=0,Dz2(0)=0',t);
end
figure
p = polyfit(zitarange, powerall, 4);
lstr = { 'Power(W)'};
yy = polyval(p, zitarange); 
%subplot(1,2,1);
plot(zitarange, powerall,'-b',zitarange,yy);
xlabel('Damping Coefficient(N·s/m)','FontSize',size);
%subplot(1,2,2);

%plot(zitarange,yy);
ylabel( '功率(W)','FontSize',size);
writematrix(powerall.',"E:\TJ\TJ2022.9\数模校赛\Q2-1power.csv");
%writematrix(x,"E:\TJ\TJ2022.9\数模校赛\Q1result-2.csv");