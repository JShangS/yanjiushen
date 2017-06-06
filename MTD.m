close all
clear all
clc
% f0=1e9;
fc=100e6;%��Ƶ
B=2e6;%����
Tao=100e-6;%����
Fs=2*B;%����Ƶ��
t=-Tao/2:1/Fs:Tao/2-1/Fs;%����ʱ��
mu=B/Tao;%��Ƶ��
R_start=0.5e4;
C=3e8;
lamda=C/fc;
delt_R=C/(2*Fs);%%�������뵥Ԫ
PRF=1000;
Tr=1/PRF;
Vr_start=100;%��ʼ�ٶ�
a=0;
pusle_num=128;%������
PV=PRF*lamda/4;
L=length(t);
Hz=linspace(0,Fs,L);
M=pusle_num;
    %%%��������
    echo=zeros(M,L);%�ز�
    echo_bu=zeros(M,L);%keystone�任�ûز�
    echo_fft=zeros(M,L);%Ƶ��ز�
    echo_bu_fft=zeros(M,L);%keystone�任��Ƶ��ز�
    pc_result_fft=zeros(M,L);%��ѹƵ���ź�
    pc_result=zeros(M,L);%��ѹʱ���ź�
    pc_result_bu=zeros(M,L);%keystone�任��Ƶʱ����ѹ
    pc_result_bu_t=zeros(M,L);%keystone�任��Ƶʱ����ѹ��ʱ����
    MTD=zeros(M,L);%δ��keystone��MTD
    MTD_bu_fft=zeros(M,L);%��keystone��MTD

for i=1:pusle_num
    Vr(i)=Vr_start+a*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%�ز��ӳ�
end
% % ��ѹϵ��
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj(fliplr(ht_t));
    ht_fft=fftshift(fft(ht));%
figure()
plot(Hz,abs(ht_fft))
%%��ѹ
tic
for i=1:pusle_num
   echo(i,:)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));
   echo_bu(i,:)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));
   echo_fft(i,:)=fftshift(fft(echo(i,:)));
   echo_bu_fft(i,:)=fftshift(fft(echo_bu(i,:)));
   pc_result(i,:)=ifft(ifftshift(echo_fft(i,:).*ht_fft));%ifftshift
   pc_result_fft(i,:)=fftshift(fft(pc_result(i,:)));%fftshift
end

% figure()
% plot(Hz,abs(echo_fft(1,:)));

figure()
surf(real(echo));
% view(0,90)

figure(), 
imagesc(abs(pc_result));
xlabel('��ʱ��')
ylabel('��ʱ��')
title('ԭʼ�ز���ѹ')
grid on
% axis([400,550,0,140,0,800])
view(0,90)
MTD=fftshift(fft(pc_result,[],1));
toc
figure()
mesh((abs(MTD)));
title('δ������λ���')
xlabel('����')
ylabel('�����յ�Ԫ')
figure
plot((abs(MTD(:,66))))
% axis([500,800,0,64])
