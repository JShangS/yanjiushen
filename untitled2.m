clear
clc
T =10e-6;  %LFM信号脉宽10us
fc = 0;  %中心频率
TBP = 150;  %时宽带宽积
K = TBP/T^2;  %线性调频率
B = K*T;  %带宽
oversamplerate=5;  %过采样率
fs = oversamplerate*B;  
Ts=1/fs;
t=(-T/2:Ts:T/2);  %采样点位置
L = length(t);
t0 = 0;  %接收时延为0


sr = rectpuls((t-t0)/T).*exp(1i*pi*K*(t-t0).^2);  %接收目标回波
h = rectpuls(t/T).*exp(-1i*pi*K*(t).^2);  %匹配滤波器
sout = conv(sr,h,'same');  %时域卷积  same是取中间部分
figure()
plot(abs(sout))

Sout = 20*log10(abs(sout)/max(abs(sout)));  %归一化并转化成dB 转dB时幅度取20 功率取10
% IRW = 0.886/(K*T) %计算IRW
% Figure 1: 基带线性调频信号的匹配滤波
% figure(1)
% subplot(2,2,1);
% plot(t/1e-6,real(sr));
% xlabel('t/us');ylabel('幅度 ');title('(a)原始信号实部 ');
% 
% subplot(2,2,2);
% % plot(t/1e-6,Sout);
% plot(abs(sout));
% % axis([-0.4 0.4 -40 1.25]);
% xlabel('t/us');ylabel('幅度 /dB');title('(b) 脉压(经扩展) ');
% 
% subplot(2,2,3);
% plot(t/1e-6,abs(sout));
% % axis([-1 1 -5 600]);
% xlabel('t/us');ylabel('幅度');title('(c) 压缩后的信号');
% 
% subplot(2,2,4);
% sout=real(sout);
% plot(t/1e-6,angle(sout))
% % axis([-1 1 -5 5]);
% xlabel('t/us');ylabel('幅度');title('(d) 压缩后信号的相位（经扩展）')

% figure(2)
% 
% sssout=T*sinc(K*T*(t-t