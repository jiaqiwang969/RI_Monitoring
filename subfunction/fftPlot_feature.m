function [the_freq,amplitude_freq]=fftPlot_dB_universal_feature(signal,fs,Scale_data,rotorspeed)
data_fft = fs;
N_overlap = data_fft/2; 
N_seg = round((length(signal)-data_fft)/(data_fft-N_overlap)) - 1;
data_freq = zeros(data_fft,size(signal,2));
for i = 1:N_seg
    data = signal(round((i-1)*(data_fft-N_overlap)+1):round((i-1)*(data_fft-N_overlap)+data_fft),:);
    temp_freq = (abs(fft(data))*2/data_fft).^2; %为了求能量
    data_freq = data_freq + temp_freq;
end
data_freq = data_freq/N_seg;
the_freq = [0:round(data_fft/Scale_data) - 1]*fs/data_fft/(rotorspeed/60);  %数据频域离散刻度
amplitude_freq=data_freq(1:round(data_fft/Scale_data),:);






% figure
% freq_dB =20*log10(data_freq./2e-5);     % 计算声压级
% freq_dB=freq_dB(1:data_fft/Scale_data,:);
% freq_dB_scale=freq_dB(1:data_fft/Scale_data,:);
% [Pulse,rotor_speed]=keyRotation(signal(:,end),fs);
%横坐标为阶次
% plot(the_freq/(rotor_speed/60), freq_dB_scale(:,sensorLoc)); hold on
% legend1 = legend('B1','R1','R2','R3','R4','R5','R6','R7','R8','C1')
% set(legend1,'NumColumns',1,'Location','northwest');
% %     ik = floor(fs/(Scale_data*rotor_speed/60*29));
% %     for nk=1:ik
% %       plot(nk*rotor_speed/60*29,min(min(freq_dB))-1,'^');
% %       hold on
% %     end
% grid on
% for k=1:length(band1)
%     line([band1(k),band1(k)],[110,180],'linestyle','--','Color','k');
% end
% set(gca,'XLim',[0 40])
% set(gca,'YLim',[110 180])
% 
% %      title({[testTime,'-FFT频谱分析'];[fname,'-转速: ',num2str(rotor_speed),'-采样率：',num2str(fs)]},'FontSize',14)
% xlabel('频率/Hz','FontSize',16);ylabel('幅值/dB','FontSize',16);
% end

