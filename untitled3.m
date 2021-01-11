%%%%%%%%%%%һά������%%%%%%%%%%%%%%
clear all
close all
clc
%% radar parameter setting
fc=1e9;
Tp=30e-6;
B=30e6;
K=B/Tp;
C=3e8;
PRI=200e-6;  %%PRIΪ200us ��Ӧ����30km
T_start=2e-6; %%���ſ���ʱ�䣬��Ӧ����300m
T_end=2e-6;   %%���Źر�ʱ�䣬��Ӧ����300m
oversamplingrate=5;  
fs=B*oversamplingrate; 





% t0=2*R0/c;
t=(-Tp/2:1/fs:Tp/2-1/fs);
%% LFM signal (���������ʽ)
sr=exp(1i*pi*K*(t).^2);   %�����ź�
% sr0=exp(1i*pi*K*(t-t0).^2);  %���ܻز��ź�
n_receive=PRI*fs-Tp*fs;
sig=[sr,zeros(1,n_receive)];

plot(imag(sr))