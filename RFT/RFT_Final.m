close all
clear all
clc
%用RFT方法,类似于MTD，用最佳多普勒去匹配，但是会带来量化误差和旁瓣
%%实验参数为：
%JIA XU, JI YU, YINGNING P, et al. 
% Radon-Fourier Transform for Radar Target Detection, 
% I-Generalized Doppler Filter bank[J]. 
% IEEE Transactions on Aerospace and Electronic. 2011(47(2)): 1183-1202.
close all
clear all
clc
% % %用RFT方法,类似于MTD，用最佳多普勒去匹配，但是会带来量化误差和旁瓣
% f0=1e9;
%%参数计算
fc=150e6;%载频
B=15e6;%带宽%64
Tao=50e-6;%脉宽%16
Fs=2*B;%采样频率%1*B
Ts=1/Fs;
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
mu=B/Tao;%条频率
C=3e8;
delt_R=C/(2*Fs)%%采样距离单元
R_start=750*delt_R;%相对初始距离
R_real=150e3;%真实初始距离
lamda=C/fc;
delt_R=C/(2*Fs);%%采样距离单元
PRF=200;
Tr=1/PRF;
Vr_start=200;%初始速度
TCI=0.5;%相参积累时间s
a=0;
pusle_num=round(TCI*PRF);%脉冲数
PV=PRF*lamda/4;
L=length(t);
Range=R_real-delt_R.*(-L/2:L/2-1);%真实距离
M=pusle_num;
PSLR1=M*lamda/2/delt_R
    %变量声明
    echo=zeros(L,M);%回波
    echo_fft=zeros(L,M);%频域回波
    pc_result_fft=zeros(L,M);%脉压频域信号
    pc_result=zeros(L,M);%脉压时域信号
for i=1:pusle_num
    Vr(i)=Vr_start+a*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%回拨延迟
end
tic
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2)).';
    ht=conj((ht_t));%fliplr
    ht_fft=fft(ht);
% % 脉压
for i=1:pusle_num
   echo(:,i)=2*exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));
   echo_fft(:,i)=(fft(echo(:,i)));
   pc_result(:,i)=ifft((echo_fft(:,i).*ht_fft));
   pc_result_fft(:,i)=(fft(pc_result(:,i)));
% %    归一化
%    pc_result(i,:)=pc_result(i,:)./max(abs(pc_result(i,:)));
end
toc
tic
MTD_fft=(fft(pc_result,[],1));
toc
% Pc=importdata('d:\Pc.txt');
% Gv_cpu=importdata('d:\Gv_cpu.txt');
% figure
% mesh(abs(Pc))
aaaa=abs(pc_result)';
figure(), 
mesh(abs(pc_result))
xlabel('慢时间/脉冲')
ylabel('快时间/距离单元')
view(0,90)
title('原始回波脉压')
grid on


% RFT法

delta_V=lamda/(2*M*Tr);
vb=lamda/(2*Tr);%第一盲速
P=3;%round(1024*4/pusle_num);%%最大盲速因子
V=-3*vb:delta_V:3*vb-delta_V;%%速度搜索
num_sou=length(V);
Gv=zeros(num_sou,L);
s_t=0;%%求积分临时量
% R00=find(max(abs(pc_result(1,:)))==abs(pc_result(1,:)));
tic
for vi=1:num_sou%速度i=L/5:L/2    i=1:L
    disp([num2str(vi/num_sou*100),' %']);
    V(vi);
    for i=1:L%初始距离单元
    indexM=round((0:M-1)*Tr*(-V(vi))/delt_R)+i;
    find0=find(indexM<1|indexM>1500);
    boola(vi,i)=isempty(find0);    
%     pc_result(index_All);
    if boola(vi,i)
        index_All=indexM+(0:M-1)*L;
        fd_t=2*V(vi)/lamda;%%搜索速度的匹配多普勒
        fdd(vi)=fd_t;
        %s_t=sum((pc_result(index_All).*exp(1j*2*pi*fd_t.*(0:M-1)*Tr)));%abs
%         real_Gv=real(pc_result(index_All)).*(cos(2*pi*fd_t.*(0:M-1)*Tr))-imag(pc_result(index_All)).*(sin(2*pi*fd_t.*(0:M-1)*Tr));
%         imag_Gv=real(pc_result(index_All)).*(sin(2*pi*fd_t.*(0:M-1)*Tr))+imag(pc_result(index_All)).*(cos(2*pi*fd_t.*(0:M-1)*Tr));
%         s_t=sum(real_Gv)+1j*sum(imag_Gv);
%         Gv(vi,i)=s_t;
        Gv(vi,i)=sum(pc_result(index_All).*exp(1j*2*pi*fd_t.*(0:M-1)*Tr));%s_t;
        s_t=0;
    end
    end 
end
toc
% save data1008.mat
% load data915.mat
% load data1008.mat
%%%画图
% figure(),
% surf(real(echo))
% view(0,90)



%%%MTD
% MTD_fft=zeros(length(t),pusle_num);
% MTD_bu_fft=zeros(length(t),pusle_num);
% MTD_fft=(fft(pc_result,[],1));
% figure()
% mesh((abs(MTD_fft)));
% title('未补偿相参积累')
% xlabel('距离')
% ylabel('多普勒单元')
figure()
[X,Y]=meshgrid(Range,V);
mesh(Y,X,abs(Gv))
xlabel('速度m/s')
ylabel('距离m')
% title('RFT图')
[Sudu,Juli]=find(abs(Gv)==max(max(abs(Gv))));
figure()
plot(V,(abs(Gv(:,Juli))));%10*log
xlabel('速度m/s')
ylabel('幅度')
aaaaaa=abs(Gv);
figure
plot(Range,abs(pc_result(:,1)))
axis([150e3-100*delt_R,150e3+100*delt_R,0,1500])
xlabel('距离m')
ylabel('幅度')
grid on
hold on
plot(Range,abs(pc_result(:,100))/max(abs(pc_result(:,100)))*1500,'r')
xlabel('距离m')
ylabel('幅度')
axis([150e3-100*delt_R,150e3+100*delt_R,0,1500])
legend('第1个回波脉压','第100个回波脉压')
%%加窗
% Gv_win=zeros(size(Gv));

Win=hann(num_sou );
for j=1:L
    Gv_win(:,j)=abs(Gv(:,j)).*Win;
end


% max_Gv=max(max(abs(Gv)));
% win_Gv=abs(Gv)./max_Gv;
% inex_1=find(win_Gv<1);
% win_Gv(inex_1)=win_Gv(inex_1)*1;
% Gv_win=win_Gv.*Gv_win;

Win2=zeros(size(Gv));
Win2(150:250,:)=0.2;
Win2(350:400,:)=0.2;
Win2(401,:)=1;
Win2(402:450,:)=0.5;
Win2(550:600,:)=20;
index_0=find(Win2==0);
Win2(index_0)=0.1;
Gv_win=Win2.*Gv_win;


figure()
[X,Y]=meshgrid(Range,V);
mesh(Y,X,abs(Gv_win))
xlabel('速度m/s')
ylabel('距离m')
zlabel('幅度')
% axis([-1000,1000,150e3-100*delt_R,150e3+100*delt_R,0,1500e2])

[Sudu,Juli]=find(abs(Gv_win)==max(max(abs(Gv_win))));
figure()
plot(V,(abs(Gv_win(:,Juli))));%10*log
xlabel('速度m/s')
ylabel('幅度')