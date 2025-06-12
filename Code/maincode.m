clear;clc
narvs = 30; % ��������
x_lb = -30*ones(1,30); % x���½�
x_ub = 30*ones(1,30); % x���Ͻ�
[x,fval,exitflag,output] = particleswarm(@fun,narvs,x_lb,x_ub);  

%����������ı仯ͼ
options = optimoptions('particleswarm','PlotFcn','pswplotbestf')   
[x,fval] = particleswarm(@fun,narvs,x_lb,x_ub,options)

%��������
options = optimoptions('particleswarm','Display','iter');
[x,fval] = particleswarm(@fun,narvs,x_lb,x_ub,options)

%�Ƽ�������
tic
options = optimoptions('particleswarm','FunctionTolerance',1e-12,'MaxStallIterations',50,'MaxIterations',20000,'HybridFcn',@fmincon);
[x,fval] = particleswarm(@fun,narvs,x_lb,x_ub,options)
toc
