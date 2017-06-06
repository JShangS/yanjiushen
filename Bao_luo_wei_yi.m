clc
clear all
close all
%%使用包络移位法%%%%%%%%%%%%%%%%%%%
fc=100e6;%载频
B=5e6;%带宽
Tao=20e-6;%脉宽
Fs=5*B;%采样频率
t=-Tao:1/Fs:Tao-1/Fs;
mu=B/Tao;%调频斜率
R_start=0.65e4;%初始距离
C=3e8;%光速
lamda=C/fc;%波长
delt_R=C/(2*Fs);%距离单元
PRF=1000;%重复频率
Tr=1/PRF;%重复周期
V_max=PRF*lamda/4;%最大不模糊速度
Vr_start=4200;%初始速度
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
view(90,90)
tic
%%%包络位移法 
d_v=lamda/(2*M*Tr);
Max_result=0;
for v=1000:d_v:6000
    pc_result_t1=pc_result;
    Move_block=round(v.*(0:M-1)*Tr./delt_R);%移动单元数
    l_t=length(Move_block);
    for i=1:M
        pc_result_t2(:,i)=circshift(pc_result_t1(:,i),Move_block(i));
    end
    MTD_t=(fft(pc_result_t2,[],2));
    Max_t=max(max(abs(MTD_t)));
    if Max_t>Max_result
        V_result=v;
        Max_result=Max_t;
        MTD_bu=MTD_t;
        pc_result_bu=pc_result_t2;
%     else
%         break;
    end
end
toc
MTD_fft=(fft(pc_result,[],2));
figure(), 
mesh(abs(pc_result));
xlabel('慢时间')
ylabel('距离单元')
view(90,90)
title('原始回波脉压')
figure()
mesh((abs(MTD_fft)));
title('未补偿相参积累')
xlabel('多普勒单元')
ylabel('距离单元')
figure(), 
mesh(abs(pc_result_bu));
xlabel('慢时间')
ylabel('距离单元')
view(90,90)
title('包络插值位移补偿回波脉压')
figure()
mesh((abs(MTD_bu)));
title('补偿相参积累')
xlabel('多普勒单元')
ylabel('距离单元')