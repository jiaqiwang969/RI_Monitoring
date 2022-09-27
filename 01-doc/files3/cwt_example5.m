clc
clear
%ECG
load wecg
fs = 180;
N = numel(wecg);
tm = 0:1/fs:N/fs-1/fs;
plot(tm,wecg)
grid on
axis tight
title('Human ECG')
xlabel('Seconds')

figure
[cfs_cwt,f_cwt] = cwt(wecg,fs);
surf(tm,f_cwt,abs(cfs_cwt),'edgecolor','none'); view(0,90);
xlabel('Seconds'); ylabel('Hz');
title('cwt')

figure;
f = linspace(1, fs/2, N/2);
Spec_cwt = cwt_cmor(wecg,1,2,f,fs);
surf(tm,f,abs(Spec_cwt),'edgecolor','none'); view(0,90);
xlabel('Seconds'); ylabel('Hz');
title('cwt\_cmor')

Spec_cwt = cwt_cmor_2(wecg,1,2,f,fs);
surf(tm,f,abs(Spec_cwt),'edgecolor','none'); view(0,90);
xlabel('Seconds'); ylabel('Hz');
title('cwt\_cmor2')


figure
wlen = 64;                         % 帧长
wind=hanning(wlen);          % 窗函数
noverlap=wlen-1;             % 重叠部分长度
nfft = wlen;                        % nfft长
[S,F,T] = spectrogram(wecg,wind,noverlap,nfft,fs);
surf(T,F,abs(S),'edgecolor','none'); view(0,90);
xlabel('Seconds'); ylabel('Hz');
title('stft-spectrogram')

inc = 1; 
[Spec_stft, f_stft, t_stft] = mystftfun(wecg, wlen, inc, nfft, fs);
figure
surf(t_stft,f_stft,abs(Spec_stft),'edgecolor','none'); view(0,90);
axis tight
xlabel('Seconds'); ylabel('Hz');
title('stft-mystftfun')


