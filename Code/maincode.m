clear;clc
narvs = 30; % 变量个数
x_lb = -30*ones(1,30); % x的下界
x_ub = 30*ones(1,30); % x的上界
[x,fval,exitflag,output] = particleswarm(@fun,narvs,x_lb,x_ub);  

%随迭代次数的变化图
options = optimoptions('particleswarm','PlotFcn','pswplotbestf')   
[x,fval] = particleswarm(@fun,narvs,x_lb,x_ub,options)

%迭代过程
options = optimoptions('particleswarm','Display','iter');
[x,fval] = particleswarm(@fun,narvs,x_lb,x_ub,options)

%推荐混合求解
tic
options = optimoptions('particleswarm','FunctionTolerance',1e-12,'MaxStallIterations',50,'MaxIterations',20000,'HybridFcn',@fmincon);
[x,fval] = particleswarm(@fun,narvs,x_lb,x_ub,options)
toc
