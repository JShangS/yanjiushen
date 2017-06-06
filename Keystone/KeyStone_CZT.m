clc
clear all
close all
%%ʹ��CZT��%%%%%%%%%%%%%%%%%%%
%%%matlab�Դ�czt�����������þ��ǲ���������DFT�������ض��ǶȽ��в���
%%%%exp(-1j*2pi*theta*/N),�����ڿ�ʱ��Ƶ��ά������ͬ���뵥ԪΪһƵ�ʵ�Ԫ
%%%ע��Ƶ�ʵ�Ķ���fftshift
fc=100e6;%��Ƶ
B=1e6;%����
Tao=20e-6;%����
Fs=5*B;%����Ƶ��
t=-Tao:1/Fs:Tao-1/Fs;
mu=B/Tao;%��Ƶб��
R_start=1.1e4;%��ʼ����
C=3e8;%����
lamda=C/fc;%����
delt_R=C/(2*Fs);%���뵥Ԫ
PRF=1000;%�ظ�Ƶ��
Tr=1/PRF;%�ظ�����
V_max=PRF*lamda/4;%���ģ���ٶ�
Vr_start=4000;%��ʼ�ٶ�
a=0;%���ٶ�
mi=7;
M=2^mi;%������
L=length(t);%%��ʱ��Ƶ����
Vr=zeros(1,M);
fd=zeros(1,M);
delt_t=zeros(1,M);
for i=1:M
    Vr(i)=Vr_start+a*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%�ز��ӳ�
end
%%��ѹ
echo=zeros(L,M);
echo_bu=zeros(L,M);
echo_fft=zeros(L,M);
echo_bu_fft=zeros(L,M);
pc_result_fft=zeros(L,M);
pc_result=zeros(L,M);
% pc_result_bu=zeros(L,M);
%���ݺ���
ht_t=exp(1j*2*pi*(mu/2*(t).^2))';
ht=conj(fliplr(ht_t));
ht_fft=(fft(ht));
%����ز�����ѹ
for i=1:M
   echo(:,i)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));%fd(i)%t+Tr*(i-1)+
   echo_bu(:,i)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));%fd(i),t+Tr*(i-1)+
   echo_fft(:,i)=fftshift(fft(echo(:,i)));%fftshift
   echo_bu_fft(:,i)=fftshift(fft(echo_bu(:,i)));%fftshift
   pc_result(:,i)=ifftshift(ifft((echo_fft(:,i).*ht_fft)));%ifftshift
   pc_result_fft(:,i)=fftshift(fft(pc_result(:,i)));%fftshift
end
figure()
mesh(abs(pc_result))
view(-90,90)
%%%Chirp_Z%%%%%%%%%%%%%%%%%
LL=2^(mi+1);%LL��Chirp��Z�任
df=B/L;
ft=(0:L-1)'*df;%-Fs/2+
alpha=(fc+ft)./fc;
n=(1:L);
W=(exp(-1j.*alpha.*2*pi/M));
eta=B/(fc*L);
m=(-LL/2:LL/2-1)';
k=6;%ģ��ϵ��
pc_result_bu_t1=zeros(L,LL);
pc_result_bu_t2=zeros(L,LL);
for i=1:L
    pc_result_bu_t1(i,:)=czt(pc_result_fft(i,:),LL,W(i));
end
for i=1:L
    for j=1:LL
       pc_result_bu_t2(i,j)=pc_result_bu_t1(i,:)*exp(1i*2*pi/LL*m*(j))*exp(1j*2*pi*k*j*(1/(1+eta*i))); %������ifft
    end
end
pc_result_bu=ifft(pc_result_bu_t2,[],1);
figure()
mesh(abs(pc_result_bu))
view(-90,90)
%%%MTD
MTD_fft=(fft(pc_result,[],2));
MTD_bu_fft=(fft(pc_result_bu,[],2));%fftshift
figure()
mesh((abs(MTD_fft)));
title('δ������λ���')
figure()
%%%ȡǰM��Ϊ������������ΪCZT�ӳ�����
mesh((abs(MTD_bu_fft(:,1:M))));%fftshift
title('���������λ���')
