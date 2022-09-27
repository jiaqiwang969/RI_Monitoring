clear all; clc; close all;

N=1024;                      % 数据长度
fs=1000;                     % 采样频率
tt=(0:N-1)/fs;               % 时间刻度
x=chirp(tt,100,1,250);       % Chirp 信号x
figure
plot(tt,x,'k'); xlim([0 max(tt)]);
xlabel('时间/s'); ylabel('幅值');
title('调频信号波形图')
set(gcf,'color','w');

wlen = 256;                         % 帧长
% inc = wlen/4;                     % 帧移
inc = 1; 
nfft = wlen;                        % nfft长
[Spec_stft, f, t] = mystftfun(x, wlen, inc, nfft, fs);
figure
imagesc(t,f,abs(Spec_stft).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('STFT谱图-mystftfun'); ylim([50 350]);

f = linspace(1, fs/2, N/2);
Spec_cwt = cwt_cmor(x,4,2,f,fs);
figure
imagesc(tt,f,abs(Spec_cwt).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt\_cmor谱图'); ylim([50 350]);

Spec_cwt = cwt_cmor_2(x,4,2,f,fs);
figure
imagesc(tt,f,abs(Spec_cwt).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt\_cmor谱图2'); ylim([50 350]);


wind=hanning(wlen);          % 窗函数
noverlap=wlen-1;             % 重叠部分长度
[B,freq,time]=spectrogram(x,wind,noverlap,wlen,fs);
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

