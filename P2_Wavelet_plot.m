%目的：导入P1_wavelet_DataBase.m提前算好的mat数据，然后作图分析

clc
clear
close all



%% 导入数据
load('16000-wavelet.mat')
RotorSpeed=16000;
% figure
% plot(global_ws{1})
% legend('B1','R1','R2','R3','R4','R5','R6','R7','R8','C1')



%% 1. 小波随阀门开度变化的趋势图（不同转速、不同传感器位置）
%a. 不同位置-传感器
sensorLabel=[1:10];

for kk=1:10
h1=figure
set(gcf,'OuterPosition',get(0,'screensize'));
jet_color=colormap(jet(length(global_ws)));
axes1 = axes('Parent',h1);
for k=1:length(global_ws)-2
    waveAmplitude=global_ws{k}(:,sensorLabel);
    plot(fk,global_ws{k}(:,kk),'.-','Color',jet_color(k,:));
    hold on
end
% title('非归一化RI指标：有效')
grid on
% 创建 ylabel
ylabel({'小波幅值谱'});
% 创建 xlabel
xlabel({'阶次'});
% 设置其余坐标区属性
set(axes1,'FontName','Helvetica Neue','FontSize',24);
title(['转速',num2str(RotorSpeed),'rpm-',sensorArray{sensorLabel(kk)},'传感器'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-',sensorArray{sensorLabel(kk)},'传感器','.png'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-',sensorArray{sensorLabel(kk)},'传感器','.fig'])
cleanfigure
matlab2tikz(['转速',num2str(RotorSpeed),'rpm-',sensorArray{sensorLabel(kk)},'传感器','.tex'],'width','\figurewidth');


close all

end

% legend('B1','R1','R2','R3','R4','R5','R6','R7','R8','C1')
