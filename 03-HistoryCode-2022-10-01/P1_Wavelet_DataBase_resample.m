%目的：由于数据量巨大，提前计算wavlet特征，并存储
%提前计算10个传感器的数据
%分析转速：9000rpm-16000rpm

clc
clear
close all




%% 主要参数
%采样率 fk 点/圈
resamplePoint=[100];
%转速     通过键相信号计算得到 rpm
RotorSpeed=9000;
% RotorSpeed=6000:500:13000;
for kk=1:length(RotorSpeed)
%小波阶次
%小波幅值
%小波能量
%阀门开度：以数据的序号呈现；
% set(axes1,'FontSize',24,'XGrid','on','XTick',[20 30 40 50 60 70 80 90 100],...
%     'XTickLabel',{'100%','90%','80%','70%','60%','50%','40%','30%','20%'});
%传感器位置： 参数object包含对应斜排阵列的10个点的位置信息，即data数据的序列号；
sensorArray={'B1';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'C1'};
sensorLabel=10;

%% 导入subfunction
addpath(genpath('subfunction'));


%% 手动选择分析数据
% [fname,location]=uigetfile({'*.mat';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
% if isstr(fname)
%     fname=cellstr(fname);
% end

%% 导入分析数据的参数说明
location=['DATA','/','Compressor2Stall-',num2str(RotorSpeed(kk))];
% location=['DATA','/','Compressor2Stall-',num2str(RotorSpeed(kk))];
load([location,'/','参数说明','/','parameter.mat']); %选择文件导入数据
disp(Note);


%% 保存图像至指定文件夹
save_directory = [strrep(location,'Database','小波谱分析'),'Phd-result-',date];  %频谱图存储文件夹
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end


%% 批量导入某个转速的mat数据
loadMat=[];
for k=1:92
    loadMat{k} = ['Compressor2Stall-',num2str(RotorSpeed(kk)),'-',num2str(k),'.mat'];
end

%% 导入数据
for k=1:length(resamplePoint)
   for  i_file=1:length(loadMat)

    DataBase=importdata(fullfile(location,char(loadMat(i_file))));
    DataBase=V2Pa(DataBase,kulite_transform_ab);

    %% 信号预处理
    sst=DataBase(:,1:end-1);

    %% 降采样处理
    sst1=resample(DataBase(:,1:end-1),resamplePoint(k)*200,fs);
    fsResample=resamplePoint(k)*200;
    %% 
    n = length(sst1);
    dt = 1/fsResample ;
    %% 小波参数
    pad = 1;                                   % 2次幂0填充
    dj = 0.125;                                % this will do 4 sub-octaves per octave
    s0 = 15*dt;                                % this says start at a scale of 6 months
    j1 = 200;                                  % this says do 7 powers-of-two with dj sub-octaves each
    lag1 = 0.75;                               % lag-1 autocorrelation for red noise background
    mother = 'morlet';

    %% 小波变换
    %fk 在【0-40】的阶次范围
    fk=[1:1:100];
    period=1./fk/(RotorSpeed(kk)/60);
    % k0=6;fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
    fourier_factor=1.0330; %在'MORLET' this is k0 (wavenumber), default is 6.
    scale=period/fourier_factor;

    [wave,period,coi] = wavelet_10sensor(sst1(:,object(1:10)),dt,pad,dj,mother,scale);
    power = (abs(wave)).^2 ;        % compute wavelet power spectrum
    fk=1./period/(RotorSpeed(kk)/60);
    global_ws{k,i_file}=reshape(sum(power,2),length(period),10);
   end
end

saveMatName = [num2str(RotorSpeed(kk)),'-caiyang','-wavelet-prestall.mat'];
save(saveMatName)
end




%% PI基准(归一化参数)和阀门开度的基准(75%-28%)
level=100;
famen_init=100;famen_end=29;famen_SI=29; %单位：百分比
xuhao_init=1;xuhao_end=92;xuhao_SI=92;

%计算阀门开度和序号的线性对应关系(famen=xuhao*a+b;)
a=(famen_end-famen_init)/(xuhao_end-xuhao_init);
b=famen_init-xuhao_init*a;


label=75:-5:25;
xuhao=round((label-b)/a);
for k=1:length(label)
    labelName{k}=[num2str(label(k))];
end

%% 1. 小波随阀门开度变化的趋势图（不同转速、不同传感器位置）
%a. 不同位置-传感器
sensorLabel=[1:10];

for kk=1:10
h1=figure
set(gcf,'OuterPosition',get(0,'screensize'));
jet_color=colormap(jet(length(global_ws)+floor(xuhao_end/10)));
axes1 = axes('Parent',h1);
for k=1:xuhao_SI
    waveAmplitude=global_ws{k}(:,sensorLabel);
    plot(fk,global_ws{k}(:,kk),'.-','Color',jet_color(k,:));
    hold on
end
% 创建 colorbar
colorbar(axes1,'TickLabels',labelName);
% title('非归一化RI指标：有效')
grid on
% 创建 ylabel
ylabel({'小波幅值谱（-）'});
% 创建 xlabel
xlabel({'阶次（f/frot）'});
% 设置其余坐标区属性
set(axes1,'FontName','Helvetica Neue','FontSize',24);
title(['转速',num2str(RotorSpeed),'rpm-',sensorArray{sensorLabel(kk)},'传感器'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-',sensorArray{sensorLabel(kk)},'传感器','.png'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-',sensorArray{sensorLabel(kk)},'传感器','.fig'])
cleanfigure
matlab2tikz(['转速',num2str(RotorSpeed),'rpm-',sensorArray{sensorLabel(kk)},'传感器','.tex'],'width','\figurewidth');
close all
end


%b. PI指标
%% 提取不同频带的小波能量
%频带1: RI频带【10-22】
%频带2: 1BPF 【27-31】
band1=[10:20];
band2=[27:31];
for i_file=1:length(global_ws)
    for k=1:size(global_ws,1)
    PI1(k,i_file,:)=sum(global_ws{k,i_file}(band1,:));
    PI2(k,i_file,:)=sum(global_ws{k,i_file}(band2,:));
    end
end
% 


h1=figure
% set(gcf,'OuterPosition',get(0,'screensize'));
axes1 = axes('Parent',h1);
plot(PI1(1,:,10),'LineWidth',2,'Color','k')
hold on
plot(PI1(1,:,9),'LineWidth',2,'Color',[0.6350 0.0780 0.1840])
plot(PI1(1,:,8),'LineWidth',2,'Color',[0.4940 0.1840 0.5560])
plot(PI1(1,:,7),'LineWidth',2,'Color',[0.3010 0.7450 0.9330])
plot(PI1(1,:,6),'LineWidth',2,'Color',[96 96 96]/255)
plot(PI1(1,:,5),'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
plot(PI1(1,:,4),'LineWidth',2,'Color','y')
plot(PI1(1,:,3),'LineWidth',2,'Color','b')
plot(PI1(1,:,2),'LineWidth',2,'Color','g')
plot(PI1(1,:,1),'LineWidth',2,'Color','r')


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




