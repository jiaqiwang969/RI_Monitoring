function coefs = cwt_cmor_2(sig,Fb,Fc,f,fs)
%利用复morlet小波进行连续小波变换
%输入：coefs(length(f),length(sig))：连续小波变换系数
%输出：一维信号sig；
%输出：Fb,Fc: 分别指morlet小波的时域衰减常数和中心频率(Hz)具体请参考cmorwavf；
%输出：f：一维频率向量(Hz)
%输出：fs：采样频率(Hz)
%

    sig = sig(:)';%强行变成y向量，避免前面出错
    L = length(sig);
    
%计算尺度
    scal = fs*Fc./f;

%计算小波
    shuaijian = 0.001;  %取小波衰减长度为0.1%
    tlow2low = sqrt(Fb*log(1/shuaijian));       %单边cmor衰减至0.1%时的时间长度，参照cmor的表达式

%小波的积分函数
    iter = 10;%小波函数的区间划分精度
    t = linspace(-tlow2low,tlow2low,2^iter);
    dt = t(2) - t(1);
    Kesi_m = cmorwavf(-tlow2low,tlow2low,2^iter,Fb,Fc);
    %卷积前准备
    t = t-t(1);
    tMax = t(end);   
    coefs = zeros(length(scal),L);      %预初设coefs

%小波与信号的卷积
    for k = 1:length(scal)    %一个scal一行
        a = scal(k);    %a是这一行的尺度函数
        j = 1+floor((0:a*tMax)/(a*dt));   %j的最大值为是确定的，尺度越大，划分的越密。相当于把一个小波拉伸的越长。
        if length(j)==1 , j = [1 1]; end
        waveinscal = Kesi_m(j);     %把积分值扩展到j区间，然后左右颠倒。f为当下尺度的积分小波函数
        coefs(k,:) = 1./sqrt(a)*wkeep1(conv(sig,waveinscal),L);      %wkeep1取diff(wconv1(ySIG,f))里长度为lenSIG的中间一段
        %
    end
 end