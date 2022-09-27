% -------------------------------------------------------------------
% exa060606.m, for fig 9.6.6
% ���� sym4 С����ʱ���Σ�ͬʱ�����˳߶Ⱥ�����ʱ���Σ� 
%--------------------------------------------------------------------
clear;

iter=3
[phi4,psi4,xval4] = wavefun('sym4',iter);
subplot(2,2,1)
plot(xval4,phi4);grid on
title('Sym4:  Phi')
subplot(2,2,2)
plot(xval4,psi4);grid on
title('Sym4:  Psi')
%
a=[1];
   fs=1;
[h,w]=freqz(phi4,a,256,'whole',fs);
hr=abs(h);
subplot(2,2,3)
plot(w(1:128),hr(1:128));grid on
title(' The FT of Phi')

[h,w]=freqz(psi4,a,256,'whole',fs);
hr=abs(h);
subplot(2,2,4)
plot(w(1:128),hr(1:128));grid on
title(' The FT of Psi')
