clear all;
% 单尺度一维离散小波变换dwt
% 格式 ① [cA,cD]=dwt(X,'wname');  ② [cA,cD]=dwt(X,Lo_D,Hi_D);
% X为被分析的离散信号；'wname'为分解所用到的小波函数，Lo_D和Hi_D为分解滤波器
% cA和cD分别为返回的低频系数和高频系数向量：length(cA)=length(cD)=length(X)/2或(length(X)+1)/2
% 如果令lx=length(X),
% lf=length(Lo_D),则length(cA)=length(cD)=floor((lx+lf-2)/2)
load noissin;  %装载原始一维信号
s = noissin(1:1000);
% 画出原始信号的波形
subplot(311);plot(s);
title('原始信号');
% 下面用haar小波函数进行一维离散小波变换
[cA1,cD1]=dwt(s,'db2');
subplot(323);plot(cA1);
ylabel('haar(cA1)');
subplot(324);plot(cD1);
ylabel('haar(cD1)');
% 给定一个小波db2，计算与之相关的分解滤波器
[Lo_D,Hi_D]=wfilters('db2','d');
% 用分解滤波器Lo_D,Hi_D，计算信号s的离散小波分解系数
[cA2,cD2]=dwt(s,Lo_D,Hi_D);
subplot(325);plot(cA2);
ylabel('db2(cA2)');
subplot(326);plot(cD2);
ylabel('db2(cD2)');