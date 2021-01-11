close all
clear all
clc
% f0=1e9;
fc=100e6;%载频
B=2e6;%带宽
Tao=100e-6;%脉宽
Fs=2*B;%采样频率
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
mu=B/Tao;%条频率R_start=0.5e4;
C=3e8;
lamda=C/fc;
L=length(t);
Hz=linspace(0,Fs,L);
sig = exp(-1j*2*pi*(mu/2*(t).^2));
sig_fft = fftshift(fft(sig));
% % 脉压系数
ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
ht=conj(fliplr(ht_t));
ht_fft=fftshift(fft(ht));%
pc=ifftshift(ifft(ht_fft.*sig_fft));
plot(abs(pc))