clc
clear all
close all
%%使用CZT法%%%%%%%%%%%%%%%%%%%
%%%matlab自带czt函数，其作用就是不再是正常DFT而是以特定角度进行采样
%%%%exp(-1j*2pi*theta*/N),还是在快时间频率维进行相同距离单元为一频率单元
%%%注意频率点的对齐fftshift
fc=100e6;%载频
B=1e6;%带宽
Tao=20e-6;%脉宽
Fs=5*B;%采样频率
t=-Tao:1/Fs:Tao-1/Fs;
mu=B/Tao;%调频斜率
R_start=1.1e4;%初始距离
C=3e8;%光速
lamda=C/fc;%波长
delt_R=C/(2*Fs);%距离单元
PRF=1000;%重复频率
Tr=1/PRF;%重复周期
V_max=PRF*lamda/4;%最大不模糊速度
Vr_start=4000;%初始速度
a=0;%加速度
mi=7;
M=2^mi;%脉冲数
L=length(t);%%快时间频点数
Vr=zeros(1,M);
fd=zeros(1,M);
delt_t=zeros(1,M);
for i=1:M
    Vr(i)=Vr_start+a*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%回拨延迟
end
%%脉压
echo=zeros(L,M);
echo_bu=zeros(L,M);
echo_fft=zeros(L,M);
echo_bu_fft=zeros(L,M);
pc_result_fft=zeros(L,M);
pc_result=zeros(L,M);
% pc_result_bu=zeros(L,M);
%传递函数
ht_t=exp(1j*2*pi*(mu/2*(t).^2))';
ht=conj(fliplr(ht_t));
ht_fft=(fft(ht));
%构造回波并脉压
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
LL=2^(mi+1);%LL点Chirp—Z变换
df=B/L;
ft=(0:L-1)'*df;%-Fs/2+
alpha=(fc+ft)./fc;
n=(1:L);
W=(exp(-1j.*alpha.*2*pi/M));
eta=B/(fc*L);
m=(-LL/2:LL/2-1)';
k=6;%模糊系数
pc_result_bu_t1=zeros(L,LL);
pc_result_bu_t2=zeros(L,LL);
for i=1:L
    pc_result_bu_t1(i,:)=czt(pc_result_fft(i,:),LL,W(i));
end
for i=1:L
    for j=1:LL
       pc_result_bu_t2(i,j)=pc_result_bu_t1(i,:)*exp(1i*2*pi/LL*m*(j))*exp(1j*2*pi*k*j*(1/(1+eta*i))); %脉冲向ifft
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
title('未补偿相参积累')
figure()
%%%取前M列为脉冲数，后面为CZT加出来的
mesh((abs(MTD_bu_fft(:,1:M))));%fftshift
title('补偿后回相参积累')
