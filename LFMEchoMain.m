% Author: Taylor-W, PHD, Electronic Engineering
%         National University Of Defence Technology
% Creation Date  : 2020-11-3
% Modification Date:
% Copyright 2020 Taylor-W.

%% Introduction
% This program is used for the second homework of synthetic aperture rada,
% including Achieve one dimensional range imaging��complete pulse 
% compression
%% �����趨
clear 
clc
fc = 1e9;  %����Ƶ��1GHz
Tp = 30e-6;  %����30us
B = 30e6;  %����30MHz
aos = 1.2; %����������
fs = aos*B;  %��������
ts = 1/fs;  %�������
K = B/Tp;  %��Ƶ��
R = [10e3,11e3,11e3+3,11e3+50];  %Ŀ��λ��
t_short = 2*R(1)/(3e8);
t_long = 2*R(4)/(3e8);
%����Ŀ��λ�ÿɼ�����С�����ʱ�ӷֱ���0.067ms��0.07367ms
%��˽����ź�ʱ��Ϊ0.067ms~Tp+0.07367ms

tr1 = t_short;  %���ղ�����ʼʱ��
tr2 = Tp+t_long+1e-6;  %���ղ�����ֹʱ�䣬Ϊ��֤��ȫ���գ�����1us
tr = (tr1:ts:tr2);  %����ʱ��
%% ����LFM�ź�
st=rectpuls((tr-tr1-Tp/2)/Tp).*exp(1i*2*pi*fc*(tr-tr1-...
              Tp/2)+1i*pi*K*(tr-tr1-Tp/2).^2);
%Figure 1�������źŲ���
% figure(1);
% plot((tr-tr1-Tp/2)/1e-6,abs(st));xlim([-Tp/1e-6 Tp/1e-6]);
% xlabel('t/us');ylabel('����');
% title('�����ź�');
%% �ز��ź�
c = 3e8;  %����
tao = 2*R/c;  %ʱ��
s = zeros(length(tao),length(tr));  %�����ڴ�
for i=1:length(tao)
    s(i,:) = rectpuls((tr-tao(i)-Tp/2)/Tp).*exp(1i*2*pi*fc*(tr-tao(i)-...
             Tp/2)+1i*pi*K*(tr-tao(i)-Tp/2).^2);  %������Ŀ��Ļز��ź�
end
S = sum(s);  %�����Ŀ���ܻز��ź�
%Figure 2���ز��źŲ���
% figure(2);
% subplot(211);
% plot(tr/1e-6,abs(s(4,:)));
% xlabel('t/us');ylabel('����');
% title('������Ŀ��ز��ź�');
% subplot(212);
% plot(tr/1e-6,abs(S));
% xlabel('t/us');ylabel('����');
% title('�ܻز��ź�');
%% �������
sout_Dm=S.*exp(-1i*2*pi*fc*tr);%���������ֱ�ӳ���exp(-1i*2*pi*fc*tr)
%% ʱ��ƥ���˲�/����ѹ��
ht=exp(-1i*pi*K*(-Tp/2:ts:Tp/2).^2);  %ƥ���˲���
sout_Dm_c=conv(sout_Dm,ht,'same');  %��ѹ
Sout_Dm_c=abs(sout_Dm_c)./max(abs((sout_Dm_c)));  %��һ��
%Figure 3�� ʱ������ѹ����һά������
figure(3);
subplot(2,1,1)
plot(tr*1e6,Sout_Dm_c);
xlabel('ʱ��/us');ylabel('��һ������');title('Demodulationʱ����ѹ');
subplot(2,1,2)
plot((tr-Tp/2).*c/2/1e3,Sout_Dm_c);
xlabel('����/Km');ylabel('��һ������');title('Demodulationһά������');
%% Ƶ������ѹ��
N = length(tr);
htt = rectpuls((tr-tr1-Tp/2)/Tp).*exp(-1i*pi*K*((tr-tr1-Tp/2)).^2);
HF = fft(htt,N);
Sout_Dm_F =( fft(sout_Dm));
sout_Dm_fc = (abs(ifft(Sout_Dm_F.*HF)));
Sout_Dm_fc = (sout_Dm_fc)./max(sout_Dm_fc);%��һ��
%Figure 4�� Ƶ����ѹһά������
figure(4);
plot((tr-Tp).*c/2/1e3,Sout_Dm_fc);
xlabel('����/Km');ylabel('��һ������');
title('DemodulationƵ����ѹһά������');
%% Dechirp����
Rref = 11e3;%�ο�����
tao_ref = 2*Rref/c;%�ο�����ʱ��
tdelta=tr-tao_ref;
s_ref = rectpuls((tr-tao_ref-Tp/2)/Tp).*exp(1i*2*pi*fc*(tr-tao_ref-...
        Tp/2)+1i*pi*K*(tr-tao_ref-Tp/2).^2);
%�ο��ź�
s_Dc = S.*conj(s_ref);%Dechirp���������ź�
S_Dc = ifftshift(ifft(s_Dc));%ֱ�ӽ���IFFT�ɵ�����ѹ�����
sout_Dc = abs(S_Dc.*exp(-1i*pi*K*tdelta.^2));%���벹����λ����RVP
sout_Dc = sout_Dc./max(sout_Dc);%��һ��
f=linspace(-fs/2,fs/2,N);
%Figure 5��Dechirp����һά������
figure(5);
plot(f*c/2/K/1e3+11,sout_Dc);%ǰһ��Ϊ��ο��������Ծ���
xlabel('����/Km');ylabel('��һ������');title('Dechirp���һά������');