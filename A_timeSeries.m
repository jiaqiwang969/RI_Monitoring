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
    [wave,period,scale,coi] = wavelet(sst(:,object(sensorLabel)),dt,pad,dj,s0,j1,mother);
    power = (abs(wave)).^2 ;        % compute wavelet power spectrum
    fk=1./period/(RotorSpeed/60);
    global_ws=sum(power.');

figure;
plot(fk,global_ws,'x-')

end













%% 对数据做时域信号分析
N=resamplePoint;             % 数据长度
fk=RotorSpeed/60*resamplePoint;
tt=(0:N-1)/fk;               % 时间刻度
figure
plot(tt, Tdata_resample.surfaces(1).v(:,object(sensorLabel)),'k'); xlim([0 max(tt)]);
xlabel('时间/s'); ylabel('幅值');
title('调频信号波形图')
set(gcf,'color','w');





%% 傅里叶变换
signal=Tdata_resample.surfaces(1).v(:,object(sensorLabel));
freq=[0:length(signal)/2.56-1]*fk/length(signal);
tmpFreq=abs(fft(signal))*2/length(signal)
figure
plot(freq/200,tmpFreq(1:length(signal)/2.56));



%% 小波变换

[y,f,coi] = cwt(DataBase(:,sensorLabel),fs,'wavetype','morlet');
figure
imagesc(tt,flip(f),abs(y).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt\_cmor谱图2'); %ylim([0 70]);


wind=hanning(wlen);          % 窗函数
noverlap=wlen-1;             % 重叠部分长度
[B,freq,time]=spectrogram(DataBase(:,object(sensorLabel)),wind,noverlap,wlen,fs);
figure
imagesc(time,freq,abs(B)); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('STFT谱图-spectrogram'); ylim([50 350]);

figure
[wt,f]=cwt(x,'amor',fs);
imagesc(tt,f,abs(wt).^2); axis xy;
xlabel('时间/s'); ylabel('对数频率/Hz');
title('cwt谱图'); 
ylim([50 350]);

figure
args = {tt,f,abs(wt).^2};
surf(args{:},'edgecolor','none');
view(0,90);
axis tight;
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt谱图'); 
ylim([50 350]);


 



