clc
clear all
close all
%%ʹ��Ƶ�򲹳���%%%%%%%%%%%%%%%%%%%
fc=1000e6;%��Ƶ
B=5e6;%����
Tao=20e-6;%����
Fs=5*B;%����Ƶ��
t=-Tao:1/Fs:Tao-1/Fs;
mu=B/Tao;%��Ƶб��
R_start=0.65e4;%��ʼ����
C=3e8;%����
lamda=C/fc;%����
delt_R=C/(2*Fs);%���뵥Ԫ
PRF=1000;%�ظ�Ƶ��
Tr=1/PRF;%�ظ�����
V_max=PRF*lamda/4;%���ģ���ٶ�
Vr_start=4200;%��ʼ�ٶ�
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
% mesh(abs(pc_result_fft))
tic
%%%����λ�Ʒ� 
d_v=lamda/(2*M*Tr);
Max_result=0;
df=Fs/length(t);
ft=(0:length(t)-1)'*df;%-B/2+
for v=4000:d_v:4500
    pc_result_t1=pc_result_fft;
    for i=1:M
        Cn=exp(-1j*2*pi*ft.*(i-1).*Tr*2*v/C);%%��������
        pc_result_t2(:,i)=pc_result_t1(:,i).*Cn;
    end
    pc_result_t3=ifft(pc_result_t2,[],1);
    MTD_t=(fft(pc_result_t3,[],2));
    Max_t=max(max(abs(MTD_t)));
    if Max_t>Max_result
        V_result=v;
        Max_result=Max_t;
        MTD_bu=MTD_t;
        pc_result_bu=pc_result_t3;
%     else
%         break;
    end
end
toc
MTD_fft=(fft(pc_result,[],2));
figure(), 
mesh(abs(pc_result));
xlabel('��ʱ��')
ylabel('��ʱ��')
view(90,90)
title('ԭʼ�ز���ѹ')
figure()
mesh((abs(MTD_fft)));
title('δ������λ���')
xlabel('����')
ylabel('�����յ�Ԫ')
figure(), 
mesh(abs(pc_result_bu));
xlabel('��ʱ��')
ylabel('���뵥Ԫ')
view(90,90)
title('Ƶ�򲹳����ز���ѹ')
figure()
mesh((abs(MTD_bu)));
title('������λ���')
xlabel('�����յ�Ԫ')
ylabel('����')