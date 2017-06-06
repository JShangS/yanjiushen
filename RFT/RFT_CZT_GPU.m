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
Fs=1*B;%采样频率B;
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
mu=B/Tao;%条频率
C=3e8;
R_max=C*Tao/2;
lamda=C/fc;
delt_R=C/(2*Fs);%%采样距离单元
R_start=100*delt_R;
PRF=500;
Tr=1/PRF;
Vr_start=1200;%初始速度1700;
aa=0;
Vb=lamda/(2*Tr);%%第一盲速
pusle_num=256;%脉冲数2048;
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
% K=-4:4;
for i=1:pusle_num
   echo(i,:)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));
%    echo(i,:)=awgn(echo(i,:),SNR);%%加噪声
end
% tic
for i=1:pusle_num
   echo_fft(i,:)=(fft(echo(i,:)));%fftshift
   pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
   pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%快时间FFT
end
% toc
bu0=zeros(M,L);
pc_bu0=cat(1,pc_result_fft,bu0);
pc_bu0_z=pc_bu0.';
%%MTD
tic
MTD_fft=(fft(pc_result,[],1));
MTD_time=toc
%%
% % 画回波图
% figure()
% surf(real(echo))
% view(0,90)
% % % 脉压后原始回波
% figure(), 
% % mesh(abs(pc_result))
% imagesc (abs(pc_result))
% xlabel('快时间')
% ylabel('慢时间')
% view(0,90)
% title('原始回波脉压')
% grid on
% juli_kuadu=find(max(abs(pc_result(1,:)))==abs(pc_result(1,:)))-find(max(abs(pc_result(pusle_num,:)))==abs(pc_result(pusle_num,:)))
% figure(), 
% mesh(abs(pc_result_fft))
%%
% % RFT_CZT法
fai=lamda*B/(L*C);
delta_V=lamda/(2*M*Tr);
%%系数生成
% V=Vb-delta_V:-delta_V:0;%%速度搜索
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=3;
x=zeros(M,L*length(K));
nn = (0:(M-1))';
ww_t=ww(M+nn,:);
ww_all=zeros(M,L*length(K));
fv_all=zeros(2*M,L*length(K));
for Nk=1:length(K)
    BU=exp(-1j*2*pi*K(Nk).*Bu_mL.*lamda*B./(L*C));%%补偿因子
    x_t= pc_result_fft.*BU;
    x(:,(Nk-1)*L+1:Nk*L)=x_t;
    ww_all(:,(Nk-1)*L+1:Nk*L)=ww_t;
    fv_all(:,(Nk-1)*L+1:Nk*L)=fv;
end
%%%%%%%%%%%%%GPU%%%%%%%%%%%%%%%%%%%%%

fv_GPU=gpuArray(fv_all);
ww_GPU=gpuArray(ww_all);
x_GPU=gpuArray(x);
tic
    y_GPU = x_GPU .* ww_GPU;
%     time_diancheng1=toc
    
%     tic
    fy_GPU = fft(y_GPU, nfft);%默认是列FFT
%     time_fft=toc
%     
%     tic
    fy_GPU = fy_GPU .* fv_GPU;
%     time_diancheng2=toc
    
%     tic
    g_GPU  = ifft(fy_GPU);
%     time_ifft1=toc
%     tic
    sp_GPU = g_GPU(M:2*M-1,:).* ww_GPU;
%     time_dicheng3=toc
    
%     tic
%     sp_GPU=sp_GPU.';
%     toc
    tic
%     sp_GPU = (ifft(sp_GPU));
    sp_GPU = (ifft(sp_GPU,[],2));
    toc
%     time_ifft2=toc
time_all=toc
    tic
    sp=gather(sp_GPU);
    time_trans=toc

    [hang,lie]=size(sp);
    trans_band=hang*lie*16/1024/1024/time_trans%传输带宽Mb/s   
    mesh((abs(sp)))