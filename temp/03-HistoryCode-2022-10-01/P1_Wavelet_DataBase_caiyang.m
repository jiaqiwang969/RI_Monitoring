%目的：由于数据量巨大，提前计算wavlet特征，并存储
%提前计算10个传感器的数据
%分析转速：9000rpm-16000rpm

clc
clear
close all




%% 主要参数
%采样率 fk 点/圈
resamplePoint=[50 60 70 80 90];
%转速     通过键相信号计算得到 rpm
RotorSpeed=16000;
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
for k=1:75
    loadMat{k} = ['Compressor2Stall-',num2str(RotorSpeed(kk)),'-',num2str(k),'.mat'];
end

%% 导入数据
for k=1:length(resamplePoint)
   for  i_file=1:length(loadMat)-2

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

h2=figure
axes1 = axes('Parent',h2);
for k=1
plot(1:length(global_ws),PI1(:,:,2))
hold on
end
legend('50','60','70','80','90')
set(axes1,'FontSize',24,'XGrid','on','XTick',[20 30 40 50 60 70 80 90 100],...
     'XTickLabel',{'100%','90%','80%','70%','60%','50%','40%','30%','20%'});
xlim([30 92])
grid on
% 创建 ylabel
ylabel({'小波幅值谱'});
% 创建 xlabel
xlabel({'阀门开度'});






h1=figure
jet_color=colormap(jet(length(resamplePoint)));
for k=1:length(resamplePoint)
%     subplot(19,1,k)
    plot(global_ws{k}(:,2),'LineWidth',2,'Color',jet_color(k,:))
    hold on
end
legend('20','25','30','35','40','45','50','55','60')

% legend('50','55','60','65','70','75','80','85','90','95','100','105','110','115','120','125','130','135','140','145','150','155','160','165','170','175','180')
ylabel({'小波幅值'});
xlabel({'阶次'});
title(['转速',num2str(RotorSpeed),'rpm'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-R1-','采样率对比','.png'])
saveas(h1,['转速',num2str(RotorSpeed),'rpm-R1-','采样率对比','.fig'])
cleanfigure
matlab2tikz(['转速',num2str(RotorSpeed),'rpm-R1-','采样率对比','.tex'],'width','\figurewidth');



