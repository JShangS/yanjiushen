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

R0=10e3; %R0�ھ���10km��


t0=2*R0/C;
t=(-Tp/2:1/fs:Tp/2-1/fs);
%% LFM signal (���������ʽ)
sr=exp(1i*pi*K*(t).^2);   %�����ź�


sr0=exp(1i*pi*K*(t-t0).^2);  %���ܻز��ź�
n_receive=PRI*fs-Tp*fs;
sig0=[sr0,zeros(1,n_receive)];

pulse=repmat(sig0,1,5);   %5��PRI
plot(real(pulse))

%% pulse compression
hmt= fliplr(sr);  %ʱ�䷴��
hmt=conj(hmt);    %ȡ������

sig0_f=fft(sig0);

n_bu=length(sig0_f-length())
hmt=[hmt,zeros(1,n_bu)]
Hmf=fft(hmt);

sout=ifft(sig0_f.*Hmf);
sout=fftshift(sout);

plot(abs(sout))