%%目的：为得到与论文对应的RI频带结构
%主要关注B1传感器



clc
clear
close all

%%
addpath(genpath('subfunction'));
addpath(genpath('Test'));
addpath(genpath('synchrosqueezing'));


%% 保存图像至指定文件夹
save_directory = ['Result-',date];  %频谱图存储文件夹
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end

%DataInfo:包含数据的位置、数据转速、数据序列、对应磁阀开度、传感器位置等信息
%RotorSpeed、xuhao_end;famen_init、famen_end
DataInfo.condition(1,:) =[6000;92-1;29;95];
DataInfo.condition(2,:) =[6500;94-1;29;95];
DataInfo.condition(3,:) =[7000;95-1;29;95];
DataInfo.condition(4,:) =[7500;93-1;29;95];
DataInfo.condition(5,:) =[8000;93-1;29;95];
DataInfo.condition(6,:) =[8500;95-1;29;95];
DataInfo.condition(7,:) =[9000;95-1;29;95];
DataInfo.condition(8,:) =[9500;95-1;29;95];
DataInfo.condition(9,:) =[10000;95-1;29;95];
DataInfo.condition(10,:)=[10500;93-1;29;95];
DataInfo.condition(11,:)=[11000;94-1;29;95];
DataInfo.condition(12,:)=[11500;97-1;29;95];
DataInfo.condition(13,:)=[12000;93-1;29;95];
DataInfo.condition(14,:)=[12500;91-1;29;95];
DataInfo.condition(15,:)=[13000;92-1;29;95];
DataInfo.condition(16,:)=[14000;67-1;29;67];
DataInfo.condition(17,:)=[15000;76-1;29;67];
DataInfo.condition(18,:)=[16000;75-1;29;67];


%目的：统一选择磁阀开度范围（横坐标固定）
%计算磁阀开度和序号的线性对应关系(famen=b-xuhao*a;)
x1=DataInfo.condition(:,4);x2=DataInfo.condition(:,3);
y1=1;y2=DataInfo.condition(:,2);
DataInfo_a=(y2-y1)./(x2-x1);
DataInfo_b=y1-DataInfo_a.*x1;

for k=1:length(DataInfo.condition)
    DataInfo.location{k,1}=['DATA','/','Compressor2Stall-',num2str(DataInfo.condition(k,1))];
    DataInfo.parameterLocation{k,1}=[DataInfo.location{k,1},'/','参数说明','/','parameter.mat'];
end

DataInfo.sensorArray={'B1';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'C1'};
%10个传感器的相对弦长位置
DataInfo.sensorLoc=[-1;5;13;23;35;47;58;71;82;101];

N_caiyang=4;%每个叶道的点数
resamplePoint=N_caiyang*29; %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
RIband=[11:25];    %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
sensorLoc=[1:10];

%转速的影响-FFT
the_freq=[];
amplitude_freq=[];
UI_FFT=[];
Process.sst=cell(length(DataInfo.condition),max(DataInfo.condition(:,2)));

%主要关注12000rpm，因为和仿真对应
for k=13
    Process.fsResample=round(resamplePoint*DataInfo.condition(k)/60/10)*10;
    
    Process.sensorParameter=load([DataInfo.location{k},'/','参数说明','/','parameter.mat']); %选择文件导入数据
    Process.sensorParameter.object(11)=[1]
Process.sensorParameter.object(12)=[2]
Process.sensorParameter.object(13)=[3]
Process.sensorParameter.object(14)=[4]
Process.sensorParameter.object(15)=[5]
Process.sensorParameter.object(16)=[6]
Process.sensorParameter.object(17)=[7]

    for  i_file=90%利用好一点的计算机并行计算
        % %         Process.loadMat=[];
        %         for k1=1:DataInfo.condition(k,2)
        %             Process.loadMat{k1} = ['Compressor2Stall-',num2str(DataInfo.condition(k)),'-',num2str(k1),'.mat'];
        %         end
        DataBase=importdata(fullfile(DataInfo.location{k},['Compressor2Stall-',num2str(DataInfo.condition(k)),'-',num2str(i_file),'.mat']));
        DataBase=V2Pa(DataBase,Process.sensorParameter.kulite_transform_ab);

        % 执行等角度采样操作:完全按照阶次谱来执行
        % n_round_Point=(Pulse(n_round+1)-Pulse(n_round))/29;
        % xuhao=Pulse(n_round)+(n_blade-1)*n_round_Point+(1/(2*N_caiyang)+n_caiyang/N_caiyang)*n_round_Point
        QUAN=112; %圈数按照最小转速(7000rpm)统一固定为112，目的是为了后续统一频率分布
        [Pulse,Rotor_Speed]=keyRotation_RealTime(DataBase(:,end),Process.sensorParameter.fs);
        xuhao=[];
        for n_round=1:QUAN
            n_round_Point=(Pulse(n_round+1)-Pulse(n_round))/29;
            for n_blade=1:29
                for n_caiyang=1:N_caiyang
                    xuhao=[xuhao round(Pulse(n_round)+(n_blade-1)*n_round_Point+(1/(2*N_caiyang)+n_caiyang/N_caiyang)*n_round_Point)];
                end
            end
        end
        sst1All{k,i_file}=DataBase(xuhao,Process.sensorParameter.object);
        %%利用流动失稳的循环平稳特性分离出tonal noise和broadband noise
        %背景流
        data_tonal_rms=permute(mean(reshape(sst1All{k,i_file},QUAN,29*N_caiyang,length(Process.sensorParameter.object)),1),[2,3,1]);
        data_tonal{k,i_file}=kron(ones(QUAN,1),data_tonal_rms);
        data_broadband{k,i_file}=sst1All{k,i_file}-data_tonal{k,i_file};
        %尝试减去方差！！看UI曲线的变化
%         data_broadband{k,i_file}= data_broadband{k,i_file}./std(data_broadband{k,i_file});
           
        %       data_diff{k,i_file}=reshape(diff(reshape(sst1All{k,i_file},QUAN,N_caiyang*29,length(Process.sensorParameter.object(sensorLoc))),1),(QUAN-1)*N_caiyang*29,length(Process.sensorParameter.object(sensorLoc)));

        %% 方法1:不平均
        P2 = 2*abs(fft(data_broadband{k,i_file})/length(data_broadband{k,i_file}));
        the_freq1{k,i_file}= N_caiyang*29*(0:(length(data_broadband{k,i_file})/2.56))/length(data_broadband{k,i_file});
        amplitude_freq1{k,i_file}= P2(1:length(data_broadband{k,i_file})/2.56+1,:);
        nlabel1=find(the_freq1{k,i_file}<RIband(end) & the_freq1{k,i_file}>RIband(1));
        UI_FFT_noAverage(k,i_file,:)=sum(amplitude_freq1{k,i_file}(nlabel1,:).^2)*(the_freq1{k,i_file}(2)-the_freq1{k,i_file}(1));

    end

    %% 方法2:stft：平均+overlap
    
        for iSensor=1:length(Process.sensorParameter.object)
            for  i_file=90%利用好一点的计算机并行计算
                %various options and parameters
                CWTopt=struct('gamma',eps,'type','morlet','mu',6,'s',2,'om',0,'nv',64,'freqscale','linear');
                STFTopt=struct('gamma',eps,'type','gauss','mu',1,'s',2.5,'om',0.5,'winlen',256*4,'squeezing','full');
                %Short-time Fourier transform (STFT)
                dt=1/(29);
                [Sx,fs,dSx] = stft_fw(data_broadband{k,i_file}(:,iSensor), dt, STFTopt);
                the_freq2{k,i_file}= 2*[0:size(Sx,1)-1]/size(Sx,1)*29;
                amplitude_freq2x{k,i_file}=mean(abs(Sx)/2,2);  %这里除以2，才使得PI结果统一？
    
                %         figure
                %         plot(the_freq1{k,i_file},amplitude_freq1{k,i_file}(:,2));
                %         hold on
                %         plot(the_freq2{k,i_file},amplitude_freq2{k,i_file});
                %
                %         % 用于验证信号
                %        xNew = stft_iw(Sx, fs, STFTopt).';
                %        figure;tplot(Sx, 1:length(data_broadband{k,i_file}(:,2)), fs); colorbar; title('STFT','FontSize',14); xlabel('Time (seconds)','FontSize',14); ylabel('Frequency (hz)', 'FontSize',14);
                %        figure(); plot(1:length(data_broadband{k,i_file}(:,2)),[data_broadband{k,i_file}(:,2),xNew(:,1)]); title('Inverse Synchrosqueezing Signal');
    
                % % STFT Synchrosqueezing transform
                %
                % [Tx, fs, Sx, Sfs, Sw, dSx] = synsq_stft_fw(1:length(data_broadband{k,i_file}(:,2)), data_broadband{k,i_file}(:,2), STFTopt);
                % xNew = synsq_stft_iw(Tx, fs, STFTopt).';
                % figure(); tplot(Tx, 1:length(data_broadband{k,i_file}(:,2)), fs); colorbar; title('STFT Synchrosqueezing','FontSize',14); xlabel('Time (seconds)','FontSize',14); ylabel('Frequency (hz)', 'FontSize',14);
                % figure(); plot(1:length(data_broadband{k,i_file}(:,2)),[data_broadband{k,i_file}(:,2),xNew(:,1)]); title('Inverse Synchrosqueezing Signal');
    
    
                %%
                nlabel2=find(the_freq2{k,i_file}<RIband(end) & the_freq2{k,i_file}>RIband(1));
                UI_FFT_stftAverage(k,i_file,iSensor)=sum(amplitude_freq2x{k,i_file}(nlabel2).^2)*(the_freq2{k,i_file}(2)-the_freq2{k,i_file}(1));
            end
            for i_file=90
                amplitude_freq2{k,i_file}(:,iSensor)=amplitude_freq2x{k,i_file};
            end
        end



    %% 方法3:wavelet：平均+overlap

%     for iSensor=1:length(Process.sensorParameter.object(sensorLoc))
%         parfor  i_file=1:DataInfo.condition(k,2)%利用好一点的计算机并行计算
%             fk=[1:0.5:58];
%             pad = 0;                                   % 2次幂0填充
%             dj = 0.125;                                % this will do 4 sub-octaves per octave
%             j1 = 200;                                  % this says do 7 powers-of-two with dj sub-octaves each
%             lag1 = 0.75;                               % lag-1 autocorrelation for red noise background
%             mother = 'morlet';
%             fsResample=DataInfo.condition(k,1)/60*N_caiyang*29;
%             dt = 1/fsResample ;
%             s0 = 15*dt;
% 
%             period=1./fk/(DataInfo.condition(k,1)/60);
%             % k0=6;fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
%             fourier_factor=1.0330; %在'MORLET' this is k0 (wavenumber), default is 6.
%             scale=period/fourier_factor;
%             fk=1./period/(DataInfo.condition(k,1)/60);
% 
%             [wave,period] = wavelet_feature(data_broadband{k,i_file},dt,pad,dj,mother,scale);
%             the_freq3{k,i_file}= fk;
%             amplitude_freq3{k,i_file}= permute(mean(abs(wave),2),[1,3,2]);
%             nlabel3=find(the_freq3{k,i_file}<RIband(end) & the_freq3{k,i_file}>RIband(1));
%             UI_FFT_waveletAverage(k,i_file,:)=sum(amplitude_freq3{k,i_file}(nlabel3,:).^2)*(the_freq3{k,i_file}(2)-the_freq3{k,i_file}(1));
%         end
%     end
% 
end

% figure
% plot(Wavlet.fk,abs(Wavlet.wave(:,1,2)))

figure;
plot(the_freq2{k,i_file},amplitude_freq2{k,i_file}(:,2))

