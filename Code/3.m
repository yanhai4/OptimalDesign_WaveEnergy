m_float=4866; %kg
r_float=1; %m
m_vib=2433; %kg
r_vib=0.5; %m
rou=1025; %kg/m**3
k=80000; %弹簧的k
zita=10000; %第一问为常数
z10=0.0; %初值
z20=0.0; %初值
l0=0.5; %m 弹簧原长
h_vib=0.5; %振子高度
h_float=3; %圆柱高
h_cone=0.8; %圆锥高
l_cen=0.0; %const中计算空心圆锥质心
g=9.8;
z_water=2.0000;
m_att=1335.535; %垂荡附加质量
f=6250; %入射波震荡频率
omiga=1.7152; %垂荡激励振幅
gama=151.4388; %兴波阻尼系数

kn=250000; %扭转弹簧的刚度
zitan=1000;%旋转阻尼系数
I_f=7001.914; %纵摇附加转动惯量
L=1690; %纵摇激励力矩振幅
alpha=8890.7; %静水恢复力矩系数
gaman=654.3383; %纵摇兴波阻尼系数
mi=0;

dertax0=m_vib*g/k;
z10=h_vib-dertax0;
zc2=(44.6944+3.8*pi*(r_float^2))/s-h_cone;
zc1=0.5*h_vib+z10
time=2*pi/omiga*40;
x0=[0 0 0 0]; %初始值；
%x0=[0 0];
% 定义 x(1)=xita1, x(2)=xita1', x(3)=xita2, x(4)=xita2';
I1=I1*100;
%dx=@(t,x)[x(2); (1/I1)*(-kn*(x(1)-x(3))-zitan*(x(2)-x(4))+m_vib*zc1*sin(x(1))*g); x(4);1/(I2+I_f)*(kn*(x(1)-x(3))+zitan*(x(2)-x(4))+L*cos(omiga*t)-gaman*x(4)-alpha*x(3))];
dx=@(t,x)[x(2); (1/I1)*(-kn*(x(1)-x(3))-zitan*(x(2)-x(4))+m_vib*zc1*sin(x(1))*g); x(4);1/(I2+I_f)*(kn*(x(1)-x(3))+zitan*(x(2)-x(4))-m_float*zc2*sin(x(3))*g+L*cos(omiga*t)-gaman*x(4)-alpha*x(3))];

%dx=@(t,x)[x(2); (-m_vib*g+(l0-x(1)+x(3))*k-10000*((abs(x(2)-x(4)))^mi)*(x(2)-x(4)))/m_vib; x(4);(rou*g* 3.14 * r_float^2 * (z_water-x(3)+1/3*h_cone)-(m_float) * g-k*(l0-x(1)+x(3))+10000*((abs(x(2)-x(4)))^mi)*(x(2)-x(4))+f * cos(omiga*t)-gama*x(4))/(m_float+m_att)];

tspan=[0,time];%time为40个波浪周期
[t,x]=ode15s(dx,tspan,x0);
dertaz=x(:,1)-x(:,3);
%[t,z1,z2]=dsolve('-m_vib * g - (z1-z2) * k - zita * (Dz1-Dz2) - m_vib * D2z1','rou * g * 3.14 * r_float ^ 2 * (z_water-z2+1/3*h_cone) - (m_att+m_float) * g + (z1-z2) * k + zita * (Dz1-Dz2) + f * cos(omiga*t)-gama*Dz2-m_float*D2z2','z1(0)=0,z2(0)=0,Dz1(0)=0,Dz2(0)=0',t);
dy=diff(x(:,4))./diff(t);


lstr = { 'Angular displacement of oscillator(rad)', 'Angular velocity of oscillator(rad/s)', 'Float angular displacement(rad)', 'Float angular velocity(rad/s)' ,'Difference in angular displacement of oscillator float(rad)'};

for i=1:length(lstr)

subplot(3,2,i)
if i<5
plot(t, x(:,i));
elseif i<6
plot(t,dertaz)
else
%xzhoulishi=[1:8630];
%plot(xzhoulishi,dy)
end

xlabel('time(s)','FontSize',size);
ylabel( lstr{i} ,'FontSize',size)

end

z1=x(:,1);
z1dao=x(:,2);
z2=x(:,3);
z2dao=x(:,4);
T=table(t,z1,z1dao,z2,z2dao);
%writetable(T,"E:\TJ\TJ2022.9\数模校赛\Q3result-1-z.csv");