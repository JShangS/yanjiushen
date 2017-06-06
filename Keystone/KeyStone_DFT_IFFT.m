clc
clear all
close all
%%ʹ��DFT-IFFT��%%%%%%%%%%%%%%%%%%%
%%%1��ע��Ƶ�ʵ�Ķ���fftshift
%%%2��ģ��ϵ�����趨����Ҫ,keystone�任�ƺ������Ƶ���й�
%%%3����PC���˵ȼ�����������һ����һ�Σ�
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
Vr_start=200;%��ʼ�ٶ�
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
% pc_result_fft=zeros(length(t),pusle_num);
pc_result=zeros(length(t),pusle_num);
% pc_result_bu=zeros(length(t),pusle_num);
%���ݺ���
ht_t=exp(1j*2*pi*(mu/2*(t).^2))';
ht=conj(fliplr(ht_t));
ht_fft=(fft(ht));
% figure
% plot(real(ht_t))
%����ز�����ѹ
for i=1:pusle_num
   echo(:,i)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));%fd(i)%t+Tr*(i-1)+
   echo_bu(:,i)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));%fd(i),t+Tr*(i-1)+
   echo_fft(:,i)=(fft(echo(:,i)));%fftshift
   echo_bu_fft(:,i)=(fft(echo_bu(:,i)));%fftshift
   pc_result(:,i)=ifft((echo_fft(:,i).*ht_fft));%ifftshift
%    pc_result_fft(:,i)=(fft(pc_result(:,i)));
   pc_result_sample(:,i)=downsample(pc_result(:,i),2);%%����
   pc_result_fft(:,i)=(fft(pc_result_sample(:,i)));%fftshift
end
% figure();
% mesh(abs(pc_result_fft))

figure();
mesh(abs(pc_result))
xlabel('��ʱ��')
ylabel('��ʱ��')
view(90,90)
title('ԭʼ�ز���ѹ')
grid on

%%%%ʹ��DFT-IFFT�����б任%%%
M=pusle_num;%������
L=length(pc_result_fft(:,1));%��ʱ�������
m=(-M/2:M/2-1)';
% l=(-L/2:L/2-1)';
% df=Fs/length(t);
% ft=-Fs/2+(0:length(t)-1)'*df;%-fs/2+
k=round(2.2*round(fd(1)/(PRF/2)));%ģ��ϵ��
% k=17;
echo_bu_fft_t2=zeros(L,M);
echo_bu_fft_t1=zeros(L,M);
eta=B/((fc*length(t)));


% % % ��S(l,k')
for i=1:L
    for j=1:M%
       echo_bu_fft_t1(i,j)=pc_result_fft(i,:)*exp(-1i*2*pi*(1+eta*i)/M*(j-M/2-1)*m); %������KeyStone
    end
end
% % % ��S(l,m')
% echo_bu_fft_t2 = ifft(echo_bu_fft_t1,[],2);
for i=1:L;
    for j=1:M
       echo_bu_fft_t2(i,j)=echo_bu_fft_t1(i,:)*exp(1i*2*pi/M*m*(j-M/2-1))*exp(1j*2*pi*k*j*(1/(1+eta*i))); %������ifft
    end
end

% figure
% mesh((abs(echo_bu_fft_t2)))
pc_result_bu=ifft(echo_bu_fft_t2,[],1);%��ʱ��ifft�任

% for i=1:pusle_num
%     pc_result_bu(:,i)=(ifft(echo_bu_fft_t2(:,i).*ht_fft));
% end


figure()
mesh((abs(pc_result_bu)))
view(90,90)
xlabel('��ʱ��')
ylabel('��ʱ��')
grid on
title('keystone��ز���ѹ')

%%%MTD
MTD_fft=zeros(length(t),pusle_num);
MTD_bu_fft=zeros(L,pusle_num);
MTD_fft=(fft(pc_result,[],2));
MTD_bu_fft=(fft(pc_result_bu,[],2));%fftshift



figure()
mesh((abs(MTD_fft)));
title('δ������λ���')
figure()
mesh((abs(MTD_bu_fft)));%fftshift
title('���������λ���')