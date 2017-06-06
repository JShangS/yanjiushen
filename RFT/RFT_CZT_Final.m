close all
clear all
clc
%用RFT方法,类似于MTD，用最佳多普勒去匹配，但是会带来量化误差和旁瓣
%%自己写的Chirp-Zmatlab中的一次只能就行一次CZT，不能够将采样间隔进行变化，所以
%%在此基础上进行的改进,能不能发个论文这个？
% f0=1e9;
fc=100e6;%载频1000e6;
B=4e6;%带宽1e6;
Tao=128e-6;%脉宽127e-6;
Fs=4*B;%采样频率B;
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
mu=B/Tao;%条频率
C=3e8;
R_max=C*Tao/2;
lamda=C/fc;
delt_R=C/(2*Fs);%%采样距离单元
R_start=100*delt_R;
PRF=400;
Tr=1/PRF;
Vr_start=200;%初始速度1700;
aa=0;
Vb=lamda/(2*Tr);%%第一盲速
pusle_num=1024;%脉冲数2048;
PV=PRF*lamda/4;
L=length(t);
M=pusle_num;
    %%
% %     变量声明
    echo=zeros(M,L);%回波
    echo_fft=zeros(M,L);%频域回波
    pc_result_fft=zeros(M,L);%脉压频域信号
    pc_result=zeros(M,L);%脉压时域信号
%%
%%%每个脉冲间隔后的速度、加速度
for i=1:pusle_num
    Vr(i)=Vr_start+aa*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%回拨延迟
end
%%
%%脉压系数
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj((ht_t));
    ht_fft=(fft(ht));%fftshift 
%%
% %%%噪声
% noise=wgn(pusle_num,L,16)+1i*wgn(pusle_num,L,16);%高斯噪声，wgn（m,n,p）,p是以dBw为单位的高斯噪声

%%
% % 脉压
SNR=10;
for i=1:pusle_num
   echo(i,:)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));
%    echo(i,:)=awgn(echo(i,:),SNR);%%加噪声
end
tic
for i=1:pusle_num
   echo_fft(i,:)=(fft(echo(i,:)));%fftshift
   pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
   pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%快时间FFT
end
toc
bu0=zeros(M,L);
pc_bu0=cat(1,pc_result_fft,bu0);
pc_bu0_z=pc_bu0.';
%%MTD
tic
MTD_fft=(fft(pc_result,[],1));
toc
%%
% % 画回波图
% figure()
% surf(real(echo))
% view(0,90)
% % 脉压后原始回波
figure(), 
% mesh(abs(pc_result))
imagesc (abs(pc_result))
xlabel('快时间')
ylabel('慢时间')
view(0,90)
title('原始回波脉压')
grid on
% juli_kuadu=find(max(abs(pc_result(1,:)))==abs(pc_result(1,:)))-find(max(abs(pc_result(pusle_num,:)))==abs(pc_result(pusle_num,:)))
% figure(), 
% mesh(abs(pc_result_fft))
%%
% % RFT_CZT法
K=2;
fai=lamda*B/(L*C);
delta_V=lamda/(2*M*Tr);
%%系数生成
V=(K(1)-0.5)*Vb:Vb/M:(K(end)+0.5)*Vb-Vb/M;%%速度搜索
num_sou=M;
Sr=zeros(num_sou,L);
m=1:M;
m=m.';
a=1-fai.*(1:L);
Cn=exp(-1j*2*pi*a./num_sou);
Cn=Cn.';
kk = ( (-M+1):M ).';
kk2 = (kk .^ 2) ./ 2;
for i=1:length(Cn)
    ww(:,i)= Cn(i).^ (kk2); 
end
v=1 ./ ww;
fv = fft(v);
nfft = 2*M;%2^nextpow2(M+M-1);
% v=zeros(nfft,L);
% K=2;%模糊因子
sp_ALL=[];
L_i=(0:L-1);
Bu_mL=m*L_i;
[~,n]=size(pc_result_fft);
tic
for Nk=1:length(K)
%     disp([num2str(Nk/length(K)*100),' %']);
    %%乘以补偿因子补偿多普勒模糊
    BU=exp(-1j*2*pi*K(Nk).*Bu_mL.*lamda*B./(L*C));%%补偿因子
    x= pc_result_fft.*BU;
    nn = (0:(M-1))';
    y = x .* ww(M+nn,:);
    fy = fft(y, nfft);%默认是列FFT
    fy = fy .* fv;
    g  = ifft(fy);
    sp=g(M:2*M-1,:).* ww(M:2*M-1,:);%(M:2*M-1,:)
    sp=(ifft(sp,[],2));
%     mesh(fftshift(abs(x_fft_czt)))
    sp_ALL=cat(1,sp_ALL,sp);
end
toc
figure()
[Y,X]=meshgrid((1:L),V');
mesh(Y,X,(abs(sp_ALL)))%sqrt
sp_jiao=ifft(sp_ALL,[],1);
figure()
imagesc((abs(sp_jiao)))
% V_all=Vb*(K(1)-1):delta_V:Vb*(K(end))-delta_V;
% % [Y_ALL,X_ALL]=meshgrid((1:L),V_all);
% figure()
% mesh((abs(sp_ALL)))
% % mesh((abs(sp_ALL)))
% % view(-90,0)
% % % % imagesc((1:L),V,abs(sp))
% xlabel('距离','FontSize',30)
% ylabel('速度','FontSize',30)
% zlabel('幅度','FontSize',30)
% title('快速RFT图')
% [hang_all,lie_all]=find(abs(sp_ALL)==max(max(abs(sp_ALL))));
% figure()
% plot(V_all,abs(sp_ALL(:,lie_all)));
% xlabel('速度')
% ylabel('幅度')
% title('RFT速度刨面图')
% figure()
% mesh((abs(MTD_fft)));
% title('未补偿相参积累')
% xlabel('距离')
% ylabel('多普勒单元')
% zlabel('幅度')
ww_z=ww.';