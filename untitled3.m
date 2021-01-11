%%%%%%%%%%%一维距离像%%%%%%%%%%%%%%
clear all
close all
clc
%% radar parameter setting
fc=1e9;
Tp=30e-6;
B=30e6;
K=B/Tp;
C=3e8;
PRI=200e-6;  %%PRI为200us 对应距离30km
T_start=2e-6; %%波门开启时间，对应距离300m
T_end=2e-6;   %%波门关闭时间，对应距离300m
oversamplingrate=5;  
fs=B*oversamplingrate; 





% t0=2*R0/c;
t=(-Tp/2:1/fs:Tp/2-1/fs);
%% LFM signal (正交解调方式)
sr=exp(1i*pi*K*(t).^2);   %发射信号
% sr0=exp(1i*pi*K*(t-t0).^2);  %接受回波信号
n_receive=PRI*fs-Tp*fs;
sig=[sr,zeros(1,n_receive)];

plot(imag(sr))