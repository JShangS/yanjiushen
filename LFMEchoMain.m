% Author: Taylor-W, PHD, Electronic Engineering
%         National University Of Defence Technology
% Creation Date  : 2020-11-3
% Modification Date:
% Copyright 2020 Taylor-W.

%% Introduction
% This program is used for the second homework of synthetic aperture rada,
% including Achieve one dimensional range imaging，complete pulse 
% compression
%% 参数设定
clear 
clc
fc = 1e9;  %中心频率1GHz
Tp = 30e-6;  %脉宽30us
B = 30e6;  %带宽30MHz
aos = 1.2; %过采样因子
fs = aos*B;  %过采样率
ts = 1/fs;  %采样间隔
K = B/Tp;  %调频率
R = [10e3,11e3,11e3+3,11e3+50];  %目标位置
t_short = 2*R(1)/(3e8);
t_long = 2*R(4)/(3e8);
%根据目标位置可计算最小和最大时延分别是0.067ms和0.07367ms
%因此接收信号时间为0.067ms~Tp+0.07367ms

tr1 = t_short;  %接收波门起始时间
tr2 = Tp+t_long+1e-6;  %接收波们终止时间，为保证完全接收，加上1us
tr = (tr1:ts:tr2);  %接收时间
%% 发射LFM信号
st=rectpuls((tr-tr1-Tp/2)/Tp).*exp(1i*2*pi*fc*(tr-tr1-...
              Tp/2)+1i*pi*K*(tr-tr1-Tp/2).^2);
%Figure 1：发射信号波形
% figure(1);
% plot((tr-tr1-Tp/2)/1e-6,abs(st));xlim([-Tp/1e-6 Tp/1e-6]);
% xlabel('t/us');ylabel('幅度');
% title('发射信号');
%% 回波信号
c = 3e8;  %光速
tao = 2*R/c;  %时延
s = zeros(length(tao),length(tr));  %分配内存
for i=1:length(tao)
    s(i,:) = rectpuls((tr-tao(i)-Tp/2)/Tp).*exp(1i*2*pi*fc*(tr-tao(i)-...
             Tp/2)+1i*pi*K*(tr-tao(i)-Tp/2).^2);  %单个点目标的回波信号
end
S = sum(s);  %多个点目标总回波信号
%Figure 2：回波信号波形
% figure(2);
% subplot(211);
% plot(tr/1e-6,abs(s(4,:)));
% xlabel('t/us');ylabel('幅度');
% title('单个点目标回波信号');
% subplot(212);
% plot(tr/1e-6,abs(S));
% xlabel('t/us');ylabel('幅度');
% title('总回波信号');
%% 正交解调
sout_Dm=S.*exp(-1i*2*pi*fc*tr);%正交解调即直接乘上exp(-1i*2*pi*fc*tr)
%% 时域匹配滤波/脉冲压缩
ht=exp(-1i*pi*K*(-Tp/2:ts:Tp/2).^2);  %匹配滤波器
sout_Dm_c=conv(sout_Dm,ht,'same');  %脉压
Sout_Dm_c=abs(sout_Dm_c)./max(abs((sout_Dm_c)));  %归一化
%Figure 3： 时域脉冲压缩和一维距离像
figure(3);
subplot(2,1,1)
plot(tr*1e6,Sout_Dm_c);
xlabel('时间/us');ylabel('归一化幅度');title('Demodulation时域脉压');
subplot(2,1,2)
plot((tr-Tp/2).*c/2/1e3,Sout_Dm_c);
xlabel('距离/Km');ylabel('归一化幅度');title('Demodulation一维距离像');
%% 频域脉冲压缩
N = length(tr);
htt = rectpuls((tr-tr1-Tp/2)/Tp).*exp(-1i*pi*K*((tr-tr1-Tp/2)).^2);
HF = fft(htt,N);
Sout_Dm_F =( fft(sout_Dm));
sout_Dm_fc = (abs(ifft(Sout_Dm_F.*HF)));
Sout_Dm_fc = (sout_Dm_fc)./max(sout_Dm_fc);%归一化
%Figure 4： 频域脉压一维距离像
figure(4);
plot((tr-Tp).*c/2/1e3,Sout_Dm_fc);
xlabel('距离/Km');ylabel('归一化幅度');
title('Demodulation频域脉压一维距离像');
%% Dechirp接收
Rref = 11e3;%参考距离
tao_ref = 2*Rref/c;%参考距离时延
tdelta=tr-tao_ref;
s_ref = rectpuls((tr-tao_ref-Tp/2)/Tp).*exp(1i*2*pi*fc*(tr-tao_ref-...
        Tp/2)+1i*pi*K*(tr-tao_ref-Tp/2).^2);
%参考信号
s_Dc = S.*conj(s_ref);%Dechirp处理的输出信号
S_Dc = ifftshift(ifft(s_Dc));%直接进行IFFT可得脉冲压缩结果
sout_Dc = abs(S_Dc.*exp(-1i*pi*K*tdelta.^2));%引入补偿相位消除RVP
sout_Dc = sout_Dc./max(sout_Dc);%归一化
f=linspace(-fs/2,fs/2,N);
%Figure 5：Dechirp接收一维距离像
figure(5);
plot(f*c/2/K/1e3+11,sout_Dc);%前一项为与参考距离的相对距离
xlabel('距离/Km');ylabel('归一化幅度');title('Dechirp输出一维距离像');