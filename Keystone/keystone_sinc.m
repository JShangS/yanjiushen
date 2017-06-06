clc
clear all
close all
%��keystone�������С������������ȶԻز�����������任
%%
fc=100e6;%��Ƶ
B=5e6;%����
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
Vr_start=1000;%��ʼ�ٶ�
a=0;%���ٶ�
pusle_num=128;%������
Vr=zeros(1,pusle_num);
fd=zeros(1,pusle_num);
delt_t=zeros(1,pusle_num);
for i=1:pusle_num
    Vr(i)=Vr_start+a*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%�ز��ӳ�
end

%  echo=zeros(pusle_num,length())
%%��ѹ
echo=zeros(length(t),pusle_num);
echo_bu=zeros(length(t),pusle_num);
echo_fft=zeros(length(t),pusle_num);
echo_bu_fft=zeros(length(t),pusle_num);
pc_result_fft=zeros(length(t),pusle_num);
pc_result=zeros(length(t),pusle_num);
pc_result_bu=zeros(length(t),pusle_num);
%���ݺ���
ht_t=exp(1j*2*pi*(mu/2*(t).^2))';
ht=conj(fliplr(ht_t));
ht_fft=(fft(ht));
%����ز�����ѹ
for i=1:pusle_num
   echo(:,i)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));%fd(i)%t+Tr*(i-1)+
   echo_bu(:,i)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));%fd(i),t+Tr*(i-1)+
   echo_fft(:,i)=fftshift(fft(echo(:,i)));
   echo_bu_fft(:,i)=fftshift(fft(echo_bu(:,i)));
   pc_result(:,i)=ifft(fftshift(echo_fft(:,i).*ht_fft));
   pc_result_fft(:,i)=fftshift(fft(pc_result(:,i)));%fftshift
end
figure();
mesh(abs(pc_result))
xlabel('��ʱ��')
ylabel('��ʱ��')
view(90,90)
title('ԭʼ�ز���ѹ')
grid on
%%%keystone,ʹ��sinc��ֵ��
df=B/length(t);
L=length(t);
M=pusle_num;
%��ʱ��Ƶ��
ft=(0:length(t)-1)'*df;%%��ʱ��任��Ķ�ӦƵ��
%����ϵ��
alpha=fc./(ft+fc);
Tm=(0:M-1);
Taom=alpha*Tm;
k=0;%round(round(fd(1)/(PRF)))-1;%ģ��ϵ��
%%%sinc��ֵ
for i=1:L
    for j=1:M%
        t1=sinc(Taom(i,:)-Tm).*kaiser(128).';
        pc_result_bu(i,j)=sum(pc_result_fft(i,:).*t1)*exp(1i*2*pi*k*j*alpha(i));%
    end
end
pc_result_bu=ifft(pc_result_bu,[],1);
figure()
mesh(abs(pc_result_bu))
% contour(abs(pc_result_bu));
xlabel('��ʱ��')
ylabel('��ʱ��')
view(90,90)
grid on
title('keystone��ز���ѹ')
%%%MTD
MTD_fft=(fft(pc_result,[],2));
MTD_bu_fft=(fft(pc_result_bu,[],2));%fftshift

figure()
mesh((abs(MTD_fft)));
title('δ������λ���')

figure()
mesh((abs(MTD_bu_fft)));%fftshift
title('���������λ���')
