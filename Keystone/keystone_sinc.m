clc
clear all
close all
%用keystone方法进行“补偿”，即先对回波进行了坐标变换
%%
fc=100e6;%载频
B=5e6;%带宽
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
Vr_start=1000;%初始速度
a=0;%加速度
pusle_num=128;%脉冲数
Vr=zeros(1,pusle_num);
fd=zeros(1,pusle_num);
delt_t=zeros(1,pusle_num);
for i=1:pusle_num
    Vr(i)=Vr_start+a*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%回拨延迟
end

%  echo=zeros(pusle_num,length())
%%脉压
echo=zeros(length(t),pusle_num);
echo_bu=zeros(length(t),pusle_num);
echo_fft=zeros(length(t),pusle_num);
echo_bu_fft=zeros(length(t),pusle_num);
pc_result_fft=zeros(length(t),pusle_num);
pc_result=zeros(length(t),pusle_num);
pc_result_bu=zeros(length(t),pusle_num);
%传递函数
ht_t=exp(1j*2*pi*(mu/2*(t).^2))';
ht=conj(fliplr(ht_t));
ht_fft=(fft(ht));
%构造回波并脉压
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
xlabel('慢时间')
ylabel('快时间')
view(90,90)
title('原始回波脉压')
grid on
%%%keystone,使用sinc插值做
df=B/length(t);
L=length(t);
M=pusle_num;
%快时间频域
ft=(0:length(t)-1)'*df;%%快时间变换后的对应频率
%伸缩系数
alpha=fc./(ft+fc);
Tm=(0:M-1);
Taom=alpha*Tm;
k=0;%round(round(fd(1)/(PRF)))-1;%模糊系数
%%%sinc插值
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
xlabel('慢时间')
ylabel('快时间')
view(90,90)
grid on
title('keystone后回波脉压')
%%%MTD
MTD_fft=(fft(pc_result,[],2));
MTD_bu_fft=(fft(pc_result_bu,[],2));%fftshift

figure()
mesh((abs(MTD_fft)));
title('未补偿相参积累')

figure()
mesh((abs(MTD_bu_fft)));%fftshift
title('补偿后回相参积累')
