%2022-09-27
%目的：该代码用于处理斜排阵列时域数据

clc
clear
close all


%% 主要参数
%采样率 fk 点/圈
resamplePoint=500;
%转速     通过键相信号计算得到 rpm
RotorSpeed=12000;
%小波阶次
%小波幅值
%小波能量
%阀门开度：以数据的序号呈现；
% set(axes1,'FontSize',24,'XGrid','on','XTick',[20 30 40 50 60 70 80 90 100],...
%     'XTickLabel',{'100%','90%','80%','70%','60%','50%','40%','30%','20%'});
%传感器位置： 参数object包含对应斜排阵列的10个点的位置信息，即data数据的序列号；
sensorArray={'B1';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'C1'};
sensorLabel=2;

%% 导入subfunction
addpath(genpath('subfunction'));


%% 选择分析数据
[fname,location]=uigetfile({'*.mat';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
if isstr(fname)
    fname=cellstr(fname);
end

%% 导入分析数据的参数说明
load([location,'/','参数说明','/','parameter.mat']); %选择文件导入数据
disp(Note);


%% 保存图像至指定文件夹
save_directory = [strrep(location,'Database','小波谱分析'),'Phd-result-',date];  %频谱图存储文件夹
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end


%% 导入数据
    figure;
    jet_color=colormap(jet(length(fname)));

for i_file=1:length(fname)

    DataBase=importdata(fullfile(location,char(fname(i_file))));
    DataBase=V2Pa(DataBase,kulite_transform_ab);

    %% 信号预处理
    sst1=DataBase(:,1:end-1);
    variance = std(sst1).^2;
    sst = (sst1 - mean(sst1))./sqrt(variance);
    n = length(sst);
    dt = 1/fs ;
    round=[0:length(sst)-1]*dt*RotorSpeed/60;  % 圈数
    xlim = [0,floor(round(end))];              % 圈数的范围
    %% 小波参数
    pad = 1;                                   % 2次幂0填充
    dj = 0.125;                                % this will do 4 sub-octaves per octave
    s0 = 15*dt;                                % this says start at a scale of 6 months
    j1 = 200;                                 % this says do 7 powers-of-two with dj sub-octaves each
    lag1 = 0.75;                               % lag-1 autocorrelation for red noise background
    mother = 'morlet';


    %% 小波变换
    %fk 在【0-40】的阶次范围
    fk=[1:1:40];
    period=1./fk/(RotorSpeed/60);
    % k0=6;fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
    fourier_factor=1.0330; %在'MORLET' this is k0 (wavenumber), default is 6.
    scale=period/fourier_factor;

    [wave,period,coi] = wavelet(sst(:,object(sensorLabel)),dt,pad,dj,mother,scale);
    power = (abs(wave)).^2 ;        % compute wavelet power spectrum
    fk=1./period/(RotorSpeed/60);
    global_ws=sum(power.');
    global_ws_norm=global_ws./max(global_ws);


    plot(fk,global_ws_norm,'.-','Color',jet_color(i_file,:));
    hold on


%% 提取不同频带的小波能量
%频带1: RI频带【10-22】
%频带2: 1BPF 【27-31】

band1=[10:22];
band2=[27:31];
PI1_norm(i_file)=sum(global_ws_norm(band1));
PI1(i_file)=sum(global_ws(band1));
PI2(i_file)=sum(global_ws(band2));



end



figure
plot(1:92,PI1,'LineWidth',4);
hold on
plot(1:92,-PI2+42000000,'LineWidth',4);
plot(1:92,PI1_norm*(1594470000/91),'LineWidth',4);

legend('RI频带:[10-22]','1BPF频带[27-31]','RI频带-归一化:[10-22]')






