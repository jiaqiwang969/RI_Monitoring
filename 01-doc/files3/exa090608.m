% -------------------------------------------------------------------
% exa060608.m, for fig 9.6.7
% 产生 meyer 小波的时域波形，同时给出了尺度函数的时域波形； 
%--------------------------------------------------------------------
clear;

lowb=-8;uppb=8;n=128;
[scaling,wavelet,t]=meyer(lowb,uppb,n);

subplot(2,2,1)
plot(t,scaling);grid on
title('meyer:  Phi')
subplot(2,2,2)
plot(t,wavelet);grid on
title('meyer: Psi')
%
a=[1];
   fs=1;
[h,w]=freqz(scaling,a,256,'whole',fs);
hr=abs(h);
subplot(2,2,3)
plot(w(1:128),hr(1:128));grid on
title(' The FT of Phi')

[h,w]=freqz(wavelet,a,256,'whole',fs);
hr=abs(h);
subplot(2,2,4)
plot(w(1:128),hr(1:128));grid on
title(' The FT of Psi')