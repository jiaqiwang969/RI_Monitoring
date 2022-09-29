%目的：由于数据量巨大，提前计算wavlet特征，并存储
%提前计算10个传感器的数据
%分析转速：9000rpm-16000rpm

clc
clear
close all




%% 主要参数
%采样率 fk 点/圈
% resamplePoint=500;
%转速     通过键相信号计算得到 rpm
RotorSpeed=14000;
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
for i_file=1:length(loadMat)

    DataBase=importdata(fullfile(location,char(loadMat(i_file))));
    DataBase=V2Pa(DataBase,kulite_transform_ab);

    %% 信号预处理
    sst=DataBase(:,1:end-1);
    n = length(sst);
    dt = 1/fs ;
    %% 小波参数
    pad = 1;                                   % 2次幂0填充
    dj = 0.125;                                % this will do 4 sub-octaves per octave
    s0 = 15*dt;                                % this says start at a scale of 6 months
    j1 = 200;                                  % this says do 7 powers-of-two with dj sub-octaves each
    lag1 = 0.75;                               % lag-1 autocorrelation for red noise background
    mother = 'morlet';


    %% 小波变换
    %fk 在【0-40】的阶次范围
    fk=[1:1:40];
    period=1./fk/(RotorSpeed(kk)/60);
    % k0=6;fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
    fourier_factor=1.0330; %在'MORLET' this is k0 (wavenumber), default is 6.
    scale=period/fourier_factor;

    [wave,period,coi] = wavelet_10sensor(sst(:,object(1:10)),dt,pad,dj,mother,scale);
    power = (abs(wave)).^2 ;        % compute wavelet power spectrum
    fk=1./period/(RotorSpeed(kk)/60);
    global_ws{i_file}=reshape(sum(power,2),length(period),10);

end

saveMatName = [num2str(RotorSpeed(kk)),'-wavelet.mat'];
save(saveMatName)
end

% figure
% plot(global_ws{i_file})
% legend('B1','R1','R2','R3','R4','R5','R6','R7','R8','C1')



