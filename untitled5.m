%%%%%%%%%%%一维距离像%%%%%%%%%%%%%%
clear all
close all
clc
%% radar parameter setting
fc=1e9;
Tp=10e-6;
B=10e6;
K=B/Tp;
% C=3e8;
PRI=200e-6;  %%PRI为200us 对应距离30km
% T_start=2e-6; %%波门开启时间，对应距离300m
% T_end=2e-6;   %%波门关闭时间，对应距离300m
aos=1000;  
fs=B*aos; 


t=0:1/fs:PRI-1/fs;

st=exp(1i*pi*K*(t-Tp/2).^2);
plot((t-Tp/2)/1e-6,imag(st))
% axis([40,60,0,1])
figure
plot(t,K*t)