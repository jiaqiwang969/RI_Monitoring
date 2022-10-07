%目的：导入P1_wavelet_DataBase.m提前算好的mat数据，然后作图分析

clc
clear
close all



%% 导入数据
load('DATA_wavelet/16000-wavelet.mat')
RotorSpeed=16000;


%% PI基准(归一化参数)和阀门开度的基准(75%-28%)
level=1e14;
famen_init=75;famen_end=28;famen_SI=29; %单位：百分比
xuhao_init=1;xuhao_end=75;xuhao_SI=74;

%计算阀门开度和序号的线性对应关系(famen=xuhao*a+b;)
a=(famen_end-famen_init)/(xuhao_end-xuhao_init);
b=famen_init-xuhao_init*a;

label=round((round(xuhao_end/10)*[0:9]*a+b)/5)*5;
xuhao=round((label-b)/a);
for k=1:10
    labelName{k}=[num2str(label(k))];
end
labelName{11}='阀门开度';


%% 1. 小波随阀门开度变化的趋势图（不同转速、不同传感器位置）
%a. PI指标
%% 提取不同频带的小波能量
%频带1: RI频带【10-22】
%频带2: 1BPF 【27-31】
band1=[10:22];
band2=[27:31];
for k=1:length(global_ws)
    PI1(k,:)=sum(global_ws{k}(band1,:))/level;
    PI2(k,:)=sum(global_ws{k}(band2,:))/level;
end

h1=figure
% set(gcf,'OuterPosition',get(0,'screensize'));
axes1 = axes('Parent',h1);
plot(PI1(:,10),'LineWidth',2,'Color','k')
hold on
plot(PI1(:,9),'LineWidth',2,'Color',[0.6350 0.0780 0.1840])
plot(PI1(:,8),'LineWidth',2,'Color',[0.4940 0.1840 0.5560])
plot(PI1(:,7),'LineWidth',2,'Color',[0.3010 0.7450 0.9330])
plot(PI1(:,6),'LineWidth',2,'Color',[96 96 96]/255)
plot(PI1(:,5),'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
plot(PI1(:,4),'LineWidth',2,'Color','y')
plot(PI1(:,3),'LineWidth',2,'Color','b')
plot(PI1(:,2),'LineWidth',2,'Color','g')
plot(PI1(:,1),'LineWidth',2,'Color','r')


%添加虚线
line([xuhao_SI,xuhao_SI],[0,max(max(PI1))*1.2],'linestyle','--','Color','k');

legend1 = legend('C1','R8','R7','R6','R5','R4','R3','R2','R1','B1','失速先兆')
set(legend1,'NumColumns',1,'Location','northwest');
set(axes1,'FontSize',20,'XGrid','on','XTick',xuhao,...
    'XTickLabel',labelName);
grid on
% 创建 ylabel
ylabel({'PI（-）'});
% 创建 xlabel
xlabel({'阀门开度（%）'});
% 设置其余坐标区属性
set(axes1,'FontName','Helvetica Neue','FontSize',24);
title(['转速',num2str(RotorSpeed),'rpm-','RI频带'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-','RI频带','.png'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-','RI频带','.fig'])
cleanfigure
matlab2tikz(['转速',num2str(RotorSpeed),'rpm-','RI频带','.tex'],'width','\figurewidth');




