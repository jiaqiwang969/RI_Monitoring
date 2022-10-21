% %出图代码，Figure：傅里叶和小波变换结果
%12000rpm
%三种不同的方法对比
%采样率116点/圈
%RIband=[8:24]
%仅对比B1和R1

clc
clear
close all


%% 参数
N_caiyang=4;%每个叶道的点数
resamplePoint=N_caiyang*29; %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
RIband=[8:24];    %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
sensorLoc=[1:10];




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
DataInfo.condition(3,:) =[7000;95-1;25;95];
DataInfo.condition(4,:) =[7500;93-1;25+1/3;95];
DataInfo.condition(5,:) =[8000;93-1;25+1/3*2;95];
DataInfo.condition(6,:) =[8500;95-1;25+1/3*3;95];
DataInfo.condition(7,:) =[9000;95-1;25+1/3*4;95];
DataInfo.condition(8,:) =[9500;95-1;25+1/3*5;95];
DataInfo.condition(9,:) =[10000;95-1;25+1/3*6;95];
DataInfo.condition(10,:)=[10500;93-1;25+1/3*7;95];
DataInfo.condition(11,:)=[11000;94-1;25+1/3*8;95];
DataInfo.condition(12,:)=[11500;97-1;25+1/3*9;95];
DataInfo.condition(13,:)=[12000;93-1;25+1/3*10;95];
DataInfo.condition(14,:)=[12500;91-1;25+1/3*11;95];
DataInfo.condition(15,:)=[13000;92-1;29;95];
DataInfo.condition(16,:)=[14000;67-1;29;70];
DataInfo.condition(17,:)=[15000;76-1;29;70];
DataInfo.condition(18,:)=[16000;75-1;29;70];



%目的：统一选择磁阀开度范围（横坐标固定）
%计算磁阀开度和序号的线性对应关系(famen=b-xuhao*a;)
x1=DataInfo.condition(:,4);x2=DataInfo.condition(:,3);
y1=1;y2=DataInfo.condition(:,2);
DataInfo_a=(y2-y1)./(x2-x1);
DataInfo_b=y1-DataInfo_a.*x1;

% 完善工作 利用上述数据信息，构造一个与UI_FFT同等长度的矩阵18*96
XUHAO=repmat([1:max(DataInfo.condition(:,2))],size(DataInfo.condition,1),1);
FAMEN=(XUHAO-DataInfo_b)./DataInfo_a;

%%

for k=1:length(DataInfo.condition)
    DataInfo.location{k,1}=['DATA','/','Compressor2Stall-',num2str(DataInfo.condition(k,1))];
    DataInfo.parameterLocation{k,1}=[DataInfo.location{k,1},'/','参数说明','/','parameter.mat'];
end

DataInfo.sensorArray={'B1';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'C1'};
%10个传感器的相对弦长位置
DataInfo.sensorLoc=[-1;5;13;23;35;47;58;71;82;101];


%转速的影响-FFT
the_freq=[];
amplitude_freq=[];
UI_FFT=[];
Process.sst=cell(length(DataInfo.condition),max(DataInfo.condition(:,2)));


for k=3:15
    Process.fsResample=round(resamplePoint*DataInfo.condition(k)/60/10)*10;
    Process.sensorParameter=load([DataInfo.location{k},'/','参数说明','/','parameter.mat']); %选择文件导入数据
    parfor  i_file=1:DataInfo.condition(k,2)%利用好一点的计算机并行计算
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
        sst1All{k,i_file}=DataBase(xuhao,Process.sensorParameter.object(sensorLoc));
        %%利用流动失稳的循环平稳特性分离出tonal noise和broadband noise
        %背景流
        data_tonal_rms=permute(mean(reshape(sst1All{k,i_file},QUAN,29*N_caiyang,length(Process.sensorParameter.object(sensorLoc))),1),[2,3,1]);
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

%     %% 方法2:stft：平均+overlap
%     %
%         for iSensor=1:length(Process.sensorParameter.object(sensorLoc))
%             parfor  i_file=1:DataInfo.condition(k,2)%利用好一点的计算机并行计算
%                 %various options and parameters
%                 CWTopt=struct('gamma',eps,'type','morlet','mu',6,'s',2,'om',0,'nv',64,'freqscale','linear');
%                 STFTopt=struct('gamma',eps,'type','gauss','mu',1,'s',2.5,'om',0.5,'winlen',256*4,'squeezing','full');
%                 %Short-time Fourier transform (STFT)
%                 dt=1/(N_caiyang*29);
%                 [Sx,fs,dSx] = stft_fw(data_broadband{k,i_file}(:,iSensor), dt, STFTopt);
%                 the_freq2{k,i_file}= 2*[0:size(Sx,1)-1]/size(Sx,1)*29;
%                 amplitude_freq2x{k,i_file}=mean(abs(Sx)/2,2);  %这里除以2，才使得PI结果统一？
%     
%                 %         figure
%                 %         plot(the_freq1{k,i_file},amplitude_freq1{k,i_file}(:,2));
%                 %         hold on
%                 %         plot(the_freq2{k,i_file},amplitude_freq2{k,i_file});
%                 %
%                 %         % 用于验证信号
%                 %        xNew = stft_iw(Sx, fs, STFTopt).';
%                 %        figure;tplot(Sx, 1:length(data_broadband{k,i_file}(:,2)), fs); colorbar; title('STFT','FontSize',20); xlabel('Time (seconds)','FontSize',20); ylabel('Frequency (hz)', 'FontSize',20);
%                 %        figure(); plot(1:length(data_broadband{k,i_file}(:,2)),[data_broadband{k,i_file}(:,2),xNew(:,1)]); title('Inverse Synchrosqueezing Signal');
%     
%                 % % STFT Synchrosqueezing transform
%                 %
%                 % [Tx, fs, Sx, Sfs, Sw, dSx] = synsq_stft_fw(1:length(data_broadband{k,i_file}(:,2)), data_broadband{k,i_file}(:,2), STFTopt);
%                 % xNew = synsq_stft_iw(Tx, fs, STFTopt).';
%                 % figure(); tplot(Tx, 1:length(data_broadband{k,i_file}(:,2)), fs); colorbar; title('STFT Synchrosqueezing','FontSize',20); xlabel('Time (seconds)','FontSize',20); ylabel('Frequency (hz)', 'FontSize',20);
%                 % figure(); plot(1:length(data_broadband{k,i_file}(:,2)),[data_broadband{k,i_file}(:,2),xNew(:,1)]); title('Inverse Synchrosqueezing Signal');
%     
%     
%                 %%
%                 nlabel2=find(the_freq2{k,i_file}<RIband(end) & the_freq2{k,i_file}>RIband(1));
%                 UI_FFT_stftAverage(k,i_file,iSensor)=sum(amplitude_freq2x{k,i_file}(nlabel2).^2)*(the_freq2{k,i_file}(2)-the_freq2{k,i_file}(1));
%             end
%             for i_file=1:DataInfo.condition(k,2)
%                 amplitude_freq2{k,i_file}(:,iSensor)=amplitude_freq2x{k,i_file};
%             end
%         end
% 
% RIband=[8:22];    %!!!!!!!!!!!!!!!!重要参数!!!!!!!!!!!!!!!!!!!!!
% 
% 
%     %% 方法3:wavelet：平均+overlap
% 
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
% % 
end

%% 不同转速的结果拟合


famen_range_fit=[70 35];
famen_range_large=[80 30];

xuhao_7000_large=find(FAMEN(3,:)<famen_range_large(1) & FAMEN(3,:)>famen_range_large(2));
xuhao_7000_fit=find(FAMEN(3,:)<famen_range_fit(1) & FAMEN(3,:)>famen_range_fit(2));

p_7000 = polyfit(FAMEN(3,famen_range_fit),UI_FFT_noAverage(3,famen_range_fit,2),5);
max_7000_35=polyval(p_7000,35);
initial_7000_70=polyval(p_7000,70);
Iri_7000_35=(max_7000_35-initial_7000_70);


%%
for k=3:15
xuhao=find(FAMEN(k,:)<famen_range_large(1) & FAMEN(k,:)>famen_range_large(2));
p(k,:) = polyfit(FAMEN(k,xuhao),UI_FFT_noAverage(k,xuhao,2),5);
end

%% 拟合初始值
initial=[];
for k=3:15
    initial=[initial polyval(p(k,:),70)];
end

figure
plot([7000:500:13000],initial)
hold on
rs=polyfit([7000:500:13000],initial,8)
plot([7000:500:13000],polyval(rs,[7000:500:13000]))

%%
xuhao_13000_fit=find(FAMEN(15,:)<famen_range_fit(1) & FAMEN(15,:)>famen_range_fit(2));
xuhao_13000_large=find(FAMEN(15,:)<famen_range_large(1) & FAMEN(15,:)>famen_range_large(2));

p_13000 = polyfit(FAMEN(15,xuhao_13000_large),UI_FFT_noAverage(15,xuhao_13000_large,2),5);
max_13000_35=polyval(p_13000,35);
max_13000_30=polyval(p_13000,30);
initial_13000_70=polyval(p_13000,70);
Iri_13000_35=1;


Iri_13000_alpha=@(alpha)  (polyval(p_13000,alpha)-initial_13000_70)/max_13000_35;
Iri_7000_alpha=@(alpha)  (polyval(p_7000,alpha)-initial_7000_70)/max_13000_35;

Isegm_alpha=@(alpha) (Iri_13000_alpha(alpha)-Iri_7000_alpha(alpha));
% Isegm_rate=@(alpha) (Iri_13000_alpha(alpha)./Iri_7000_alpha(alpha));

%利用等比数列进行近似
%需要知道a=alpha_13000-alpha_12500;
%r=(E_7500-E_7000)/(E_13000-E_125000)
% Final(13,:)=[13000,326374];
% Final(12,:)=[12500,272559];
% Final(2,:)=[7500,38421];
% Final(1,:)=[7000,30742];
% r=((38421-30742)/(326374-272559))^(1/11);
%已知：Sn=I_segm(Isegm_alpha)和r
%等比数列求和公式：Sn=a(1-r^n)/(1-r);
q=0.93%linspace(0.99,0.8,20);
 for kq=1:length(q)
% a=@(alpha) 12*Isegm_alpha(alpha)*(1-q(kq))/(1-q(kq)^((13000-7000)/500));

Isegm_nihe=@(alpha) (Iri_13000_alpha(alpha)-Iri_7000_alpha(alpha))/(Iri_13000_alpha(35)-Iri_7000_alpha(35))

%若已知两条曲线的所有值
Iri_1 = @(RS, alpha) (alpha<70).*(Iri_13000_alpha(alpha)-Isegm_alpha(alpha)/(1-q(kq)^12)*(1-q(kq)^((13000-RS)/500)));
xishu=(Iri_1(13000,35)/Iri_1(13000,30));
Iri_2 = @(RS, alpha) Iri_1(RS, alpha)*xishu;%将最高点从35移回到30；
% Isegm_alpha(alpha).*(13000-RS)/500;%.*((1-(alpha-famen_range_fit(2))./(famen_range_fit(1)-famen_range_fit(2)))).^(2);
Iri_1_real = @(RS, alpha) Iri_1(RS, alpha).*max_13000_35+polyval(rs,RS);
% %若仅知道13000rpm曲线的所有值，以及7000的一个点
% Iri_2 = @(RS, alpha) Iri_13000_alpha(alpha)-Isegm_alpha(35).*(13000-RS)/500.*((1-(alpha-famen_range_fit(2))./(famen_range_fit(1)-famen_range_fit(2))).^2);
% Iri_2_real = @(RS, alpha) Iri_2(RS, alpha).*max_13000_35+polyval(rs,RS);


% %变化规律并非线性
% figure
% plot(linspace(0,1,36),Isegm_nihe([70:-1:35]))
% %https://cdn.mathpix.com/snip/images/J8nsUEaecZA59onYyQbgkxOmxis5DArczc026apaAR4.original.fullsize.png
% hold on
% plot(linspace(0,1,36),linspace(0,1,36).^3)


pp = polyfit([70:-1:35],Isegm_nihe([70:-1:35]),2)

% (1-(alpha-famen_range_fit(2))./(famen_range_fit(1)-famen_range_fit(2))

% %% 最终拟合曲线的结果
% famen_range=[95 24];
% alpha1=famen_range(1):-1:famen_range(2);
% figure
% for k=3:15
%     RS1=DataInfo.condition(k,1);
%     xuhao=find(FAMEN(k,:)<famen_range(1) & FAMEN(k,:)>famen_range(2));
%     plot(FAMEN(k,xuhao),(UI_FFT_noAverage(k,xuhao,2)-Iri_1_real(RS1,70))/max_13000_35*xishu)
%     hold on
%     plot(alpha1,Iri_2(RS1,alpha1))
% end
% set(gca, 'XDir','reverse')
% xlim([25 90])

% %% 还得将高度统一到30%
% famen_range=[95 24];
% alpha1=famen_range(1):-1:famen_range(2);
% figure
% xishu=(Iri_1(13000,35)/Iri_1(13000,30));
% for k=3:15
%     RS1=DataInfo.condition(k,1);
%     xuhao=find(FAMEN(k,:)<famen_range(1) & FAMEN(k,:)>famen_range(2));
%     plot(FAMEN(k,xuhao),(UI_FFT_noAverage(k,xuhao,2)-Iri_1_real(RS1,70))/max_13000_35*xishu)
%     hold on
%     plot(alpha1,Iri_2(RS1,alpha1))
% end
% set(gca, 'XDir','reverse')
% xlim([25 90])



%% 误差分析
% %每个转速都求一个均方差
% famen_range=[90 27]
% figure;
% 
% for k=3:15
%     RS1=DataInfo.condition(k,1);
%     xuhao=find(FAMEN(k,:)<famen_range(1) & FAMEN(k,:)>famen_range(2));
%     I_RI_exp=(UI_FFT_noAverage(k,xuhao,2)-Iri_1_real(RS1,70))./Iri_1_real(13000,30);
% plot(FAMEN(k,xuhao),I_RI_exp)
% hold on
% plot(FAMEN(k,xuhao),Iri_2(RS1,FAMEN(k,xuhao)))
% %     hold on
% %     plot(alpha1,Iri_1(RS1,alpha1))
% end
% set(gca, 'XDir','reverse')
% 

%% 数据
famen_range=[95 31];
alpha1=famen_range(1):-1:famen_range(2);
% figure
xishu=(Iri_1(13000,35)/Iri_1(13000,30));
Variance=[];
for k=3:15
    RS1=DataInfo.condition(k,1);
    xuhao=find(FAMEN(k,:)<famen_range(1) & FAMEN(k,:)>famen_range(2));
    I_RI_exp=(UI_FFT_noAverage(k,xuhao,2)-Iri_1_real(RS1,70))/max_13000_35*xishu;
    I_RI_base=Iri_2(RS1,FAMEN(k,xuhao));
%     hold on
%     plot(FAMEN(k,xuhao),(I_RI_exp-I_RI_base).^2);
    Variance(k)=sqrt(sum((I_RI_exp-I_RI_base).^2)/length(I_RI_exp));
    Max(k)=I_RI_base(end);
end
set(gca, 'XDir','reverse')
xlim([25 90])


h1=figure('Visible', 'on');
% set desired output size
set(h1, 'Units','centimeters')
height = 18*1.3;
width = 25*1.3;

% the last two parameters of 'Position' define the figure size
set(h1, 'Position',[25 5 width height],...
    'PaperSize',[width height],...
    'PaperPositionMode','auto',...
    'InvertHardcopy', 'off',...
    'Renderer','painters'...     %recommended if there are no alphamaps
    );
set(gcf, 'color', 'white');

subplot(2,1,1)
bar(DataInfo.condition(3:15,1),Variance(3:15)./Max(3:15))
max_variance(kq)=max(Variance(3:15)./Max(3:15));
title(['公比q=',num2str(q(kq)),'-','Max(\sigma)=',num2str(max_variance(kq))],'FontSize',20)
% ylabel('$\sigma=\sqrt{\sum_{i}(I_{RIexp}(i)-I_{RIbase}(i))^2/N}/I_{RIbase}(N)$','interpreter','latex','FontSize',14)
ylabel('方差\sigma','FontSize',20)
x=1;
hold off
ylim([0 0.16])
set(gca,'FontSize',20);


subplot(2,1,2)
figure
famen_range=[95 24];
alpha1=famen_range(1):-1:famen_range(2);

xishu=(Iri_1(13000,35)/Iri_1(13000,30));
for k=3:15
    RS1=DataInfo.condition(k,1);
    xuhao=find(FAMEN(k,:)<famen_range(1) & FAMEN(k,:)>famen_range(2));
    plot(FAMEN(k,xuhao),(UI_FFT_noAverage(k,xuhao,2)-Iri_1_real(RS1,70))/max_13000_35*xishu)
    hold on
    plot(alpha1,Iri_2(RS1,alpha1))
end
set(gca, 'XDir','reverse')
xlim([24 75])
ylim([0 1])
ylabel('I_{RI}','FontSize',20)
xlabel('磁阀开度(%)','FontSize',20);
set(gca,'FontSize',20);




saveas(h1,[save_directory,'/','FIG05-variance-','q',num2str(q(kq)),'-',num2str(max_variance(kq)),'-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'pdf')
% saveas(h1,[save_directory,'/','FIG04-caiyanglv-',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'fig')
% saveas(h1,[save_directory,'/','FIG04-caiyanglv-',num2str(DataInfo.condition(i_rotorspeed,1)),'rpm','-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'pdf')


end



h2=figure
bar(q,max_variance)
xlabel('公比q(Isegm为等比数列)')
ylabel('最大方差值')
saveas(h2,[save_directory,'/','FIG05-maxvariance-','-',num2str(max_variance(kq)),'-采样率',num2str(resamplePoint),'-RIband-',num2str(RIband(1)),'-',num2str(RIband(end))],'pdf')





















