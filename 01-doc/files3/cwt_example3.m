clc
clear
load batsignal
N = numel(batsignal);
t = 0:DT:(N*DT)-DT;
fs = 1/DT;
f = linspace(1, fs/2, N/2);
fft_sig = fft(batsignal);
figure
subplot(211), plot(t,batsignal); xlabel('t/Sec'); 
subplot(212),plot(f,abs(fft_sig(1:N/2)));  xlabel('Freq/Hz')

Spec_cwt = cwt_cmor(batsignal,1,4,f,fs);
figure
imagesc(t,f,abs(Spec_cwt).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt\_cmor谱图');

Spec_cwt = cwt_cmor_2(batsignal,1,4,f,fs);
figure
imagesc(t,f,abs(Spec_cwt).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt\_cmor谱图2');

[cfs,f] = cwt(batsignal,'bump',1/DT,'VoicesPerOctave',32);
args = {t,f,abs(cfs).^2};
figure
surf(args{:},'edgecolor','none');
view(0,90);
axis tight;
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt谱图');

wlen = 64;                         % 帧长
% inc = wlen/4;                       % 帧移
inc = 1; 
nfft = wlen*2;                        % nfft长
[Spec_stft, f, t_stft] = mystftfun(batsignal, wlen, inc, nfft, fs);
figure
imagesc(t_stft,f,abs(Spec_stft).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('STFT谱图-mystftfun'); 

wind=hanning(wlen);          % 窗函数
noverlap=wlen-1;             % 重叠部分长度
[B,freq,time]=spectrogram(batsignal,wind,noverlap,nfft,fs);
figure
imagesc(time,freq,abs(B)); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('STFT谱图-spectrogram'); 
