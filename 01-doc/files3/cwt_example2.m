clear all; clc; close all;

fs = 200;
t = 0 : 1/fs : 15;

c1 = 2 * pi * 10;            % initial frequency of the chirp excitation
c2 = 2 * pi * 5/2;           % set the speed of frequency change be 1 Hz/second
c3 = 2 * pi * 1/3;
c4 = 2 * pi * -1/40;

Sig = sin(c1 * t + c2 * t.^2 / 2 + c3 * t.^3 /3 + c4 * t.^4 /4);  
N = length(Sig);
IF = (c1 + c2*t + c3 *t.^2 + c4*t.^3)/2/pi;

figure
subplot(311), plot(t,Sig); xlabel('t/Sec'); 
subplot(312), plot(t,IF);  xlabel('t/Sec'); ylabel('Freq/Hz');
fft_Sig = fft(Sig);
subplot(313),plot([0:N/fs:N/fs*(N/2-1)],abs(fft_Sig(1:N/2)));  xlabel('Freq/Hz')

wlen = 256;                         % 帧长
% inc = wlen/4;                       % 帧移
inc = 1; 
nfft = wlen;                        % nfft长
[Spec_stft, f, t_stft] = mystftfun(Sig, wlen, inc, nfft, fs);
figure
imagesc(t_stft,f,abs(Spec_stft).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('STFT谱图-mystftfun'); 
ylim([10 50]);

f = linspace(1, fs/2, N/2);
Spec_cwt = cwt_cmor(Sig,1,4,f,fs);
figure
imagesc(t,f,abs(Spec_cwt).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt\_cmor谱图');
ylim([10 50]);

Spec_cwt = cwt_cmor_2(Sig,1,4,f,fs);
figure
imagesc(t,f,abs(Spec_cwt).^2); axis xy;
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt\_cmor谱图2');
ylim([10 50]);


wind=hanning(wlen);          % 窗函数
noverlap=wlen-1;             % 重叠部分长度
[B,freq,time]=spectrogram(Sig,wind,noverlap,wlen,fs);
figure
imagesc(time,freq,abs(B)); axis xy;
ylim([10 50]);
xlabel('时间/s'); ylabel('频率/Hz');
title('STFT谱图-spectrogram'); 


figure
[wt,f]=cwt(Sig,'amor',fs);
imagesc(t,f,abs(wt)); axis xy;
xlabel('时间/s'); ylabel('对数频率/Hz');
title('cwt谱图'); 


figure
args = {t,f,abs(wt).^2};
surf(args{:},'edgecolor','none');
view(0,90);
ylim([10 60]);
xlabel('时间/s'); ylabel('频率/Hz');
title('cwt谱图'); 