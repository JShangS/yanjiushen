clc
clear all
close all
%%使用DFT-IFFT法%%%%%%%%%%%%%%%%%%%
%%%1、注意频率点的对齐fftshift
%%%2、模糊系数的设定很重要,keystone变换似乎与采样频率有关
%%%3、对PC做了等间隔抽样（间隔一个抽一次）
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
Vr_start=200;%初始速度
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
% pc_result_fft=zeros(length(t),pusle_num);
pc_result=zeros(length(t),pusle_num);
% pc_result_bu=zeros(length(t),pusle_num);
%传递函数
ht_t=exp(1j*2*pi*(mu/2*(t).^2))';
ht=conj(fliplr(ht_t));
ht_fft=(fft(ht));
% figure
% plot(real(ht_t))
%构造回波并脉压
for i=1:pusle_num
   echo(:,i)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));%fd(i)%t+Tr*(i-1)+
   echo_bu(:,i)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));%fd(i),t+Tr*(i-1)+
   echo_fft(:,i)=(fft(echo(:,i)));%fftshift
   echo_bu_fft(:,i)=(fft(echo_bu(:,i)));%fftshift
   pc_result(:,i)=ifft((echo_fft(:,i).*ht_fft));%ifftshift
%    pc_result_fft(:,i)=(fft(pc_result(:,i)));
   pc_result_sample(:,i)=downsample(pc_result(:,i),2);%%抽样
   pc_result_fft(:,i)=(fft(pc_result_sample(:,i)));%fftshift
end
% figure();
% mesh(abs(pc_result_fft))

figure();
mesh(abs(pc_result))
xlabel('慢时间')
ylabel('快时间')
view(90,90)
title('原始回波脉压')
grid on

%%%%使用DFT-IFFT法进行变换%%%
M=pusle_num;%脉冲数
L=length(pc_result_fft(:,1));%快时间采样数
m=(-M/2:M/2-1)';
% l=(-L/2:L/2-1)';
% df=Fs/length(t);
% ft=-Fs/2+(0:length(t)-1)'*df;%-fs/2+
k=round(2.2*round(fd(1)/(PRF/2)));%模糊系数
% k=17;
echo_bu_fft_t2=zeros(L,M);
echo_bu_fft_t1=zeros(L,M);
eta=B/((fc*length(t)));


% % % 求S(l,k')
for i=1:L
    for j=1:M%
       echo_bu_fft_t1(i,j)=pc_result_fft(i,:)*exp(-1i*2*pi*(1+eta*i)/M*(j-M/2-1)*m); %脉冲向KeyStone
    end
end
% % % 求S(l,m')
% echo_bu_fft_t2 = ifft(echo_bu_fft_t1,[],2);
for i=1:L;
    for j=1:M
       echo_bu_fft_t2(i,j)=echo_bu_fft_t1(i,:)*exp(1i*2*pi/M*m*(j-M/2-1))*exp(1j*2*pi*k*j*(1/(1+eta*i))); %脉冲向ifft
    end
end

% figure
% mesh((abs(echo_bu_fft_t2)))
pc_result_bu=ifft(echo_bu_fft_t2,[],1);%快时间ifft变换

% for i=1:pusle_num
%     pc_result_bu(:,i)=(ifft(echo_bu_fft_t2(:,i).*ht_fft));
% end


figure()
mesh((abs(pc_result_bu)))
view(90,90)
xlabel('慢时间')
ylabel('快时间')
grid on
title('keystone后回波脉压')

%%%MTD
MTD_fft=zeros(length(t),pusle_num);
MTD_bu_fft=zeros(L,pusle_num);
MTD_fft=(fft(pc_result,[],2));
MTD_bu_fft=(fft(pc_result_bu,[],2));%fftshift



figure()
mesh((abs(MTD_fft)));
title('未补偿相参积累')
figure()
mesh((abs(MTD_bu_fft)));%fftshift
title('补偿后回相参积累')