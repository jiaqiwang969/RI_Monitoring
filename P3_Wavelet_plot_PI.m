%目的：导入P1_wavelet_DataBase.m提前算好的mat数据，然后作图分析

clc
clear
close all



%% 导入数据
load('13000-wavelet.mat')
RotorSpeed=13000;
% figure
% plot(global_ws{1})
% legend('B1','R1','R2','R3','R4','R5','R6','R7','R8','C1')



%% 1. 小波随阀门开度变化的趋势图（不同转速、不同传感器位置）
%a. PI指标
%% 提取不同频带的小波能量
%频带1: RI频带【10-22】
%频带2: 1BPF 【27-31】
band1=[10:22];
band2=[27:31];
for k=1:length(global_ws)
    PI1(k,:)=sum(global_ws{k}(band1,:));
    PI2(k,:)=sum(global_ws{k}(band2,:));
end

h1=figure
axes1 = axes('Parent',h1);
plot(PI1(:,10),'LineWidth',4)
hold on
plot(PI1(:,9),'LineWidth',4,'Color','k')
plot(PI1(:,8),'LineWidth',4)
plot(PI1(:,7),'LineWidth',4)
plot(PI1(:,6),'LineWidth',4)
plot(PI1(:,5),'LineWidth',4,'Color','b')
plot(PI1(:,4),'LineWidth',4,'Color','c')
plot(PI1(:,3),'LineWidth',4,'Color','y')
plot(PI1(:,2),'LineWidth',4,'Color','g')
plot(PI1(:,1),'LineWidth',4,'Color','r')


legend1 = legend('C1','R8','R7','R6','R5','R4','R3','R2','R1','B1')
set(legend1,'NumColumns',2,'Location','northwest');
set(axes1,'FontSize',20,'XGrid','on','XTick',[20 30 40 50 60 70 80 90 100],...
    'XTickLabel',{'100%','90%','80%','70%','60%','50%','40%','30%','20%'});
grid on
% 创建 ylabel
ylabel({'PI'});
% 创建 xlabel
xlabel({'阀门开度'});
% 设置其余坐标区属性
set(axes1,'FontName','Helvetica Neue','FontSize',24);
title(['转速',num2str(RotorSpeed),'rpm-','RI频带'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-','RI频带','.png'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-','RI频带','.fig'])
cleanfigure
matlab2tikz(['转速',num2str(RotorSpeed),'rpm-','RI频带','.tex'],'width','\figurewidth');




