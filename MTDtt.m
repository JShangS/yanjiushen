close all
clear all
clc
fc=1e9;
% fc=0;%载频
B=2e6;%带宽
Tao=100e-6;%脉宽
Fs=2*B;%采样频率
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
mu=B/Tao;%条频率
R_start=0.5e4;
C=3e8;
lamda=C/fc;
delt_R=C/(2*Fs);%%采样距离单元

PRF=1000;
Tr=1/PRF;
Vr_start=300;%初始速度
a=0;
pusle_num=128;%脉冲数
PV=PRF*lamda/4;
L=length(t);
R_cell = L-round( R_start/delt_R);
Hz=linspace(0,Fs,L);
M=pusle_num;
    %%%变量声明
    echo=zeros(M,L);%回波
    echo_bu=zeros(M,L);%keystone变换用回波
    echo_fft=zeros(M,L);%频域回波
    echo_bu_fft=zeros(M,L);%keystone变换用频域回波
    pc_result_fft=zeros(M,L);%脉压频域信号
    pc_result=zeros(M,L);%脉压时域信号
    pc_result_bu=zeros(M,L);%keystone变换后频时域脉压
    pc_result_bu_t=zeros(M,L);%keystone变换后频时域脉压临时变量
    MTD=zeros(M,L);%未作keystone的MTD
    MTD_bu_fft=zeros(M,L);%作keystone的MTD

for i=1:pusle_num
    Vr(i)=Vr_start+a*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%回拨延迟
end
% % 脉压系数
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj(fliplr(ht_t));
    ht_fft=fftshift(fft(ht));%
figure()
plot(Hz,abs(ht_fft))
%%脉压
tic
for i=1:pusle_num
   echo(i,:)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));
%    echo(i,:) = awgn(echo(i,:), -10);
   echo_bu(i,:)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));
   echo_fft(i,:)=fftshift(fft(echo(i,:)));
   echo_bu_fft(i,:)=fftshift(fft(echo_bu(i,:)));
   pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
   pc_result_fft(i,:)=fftshift(fft(pc_result(i,:)));%fftshift
end

%% MTI
for i=2:pusle_num-1
    MTI_tt(i-1,:) = 2*pc_result(i,:)-pc_result(i-1,:)-pc_result(i+1,:);
end
MTI =fftshift( fft(MTI_tt,[],1));

figure
imagesc((abs(MTI)))
figure()
plot(Hz,abs(echo_fft(1,:)));

figure()
surf(real(echo));
% view(0,90)

figure(), 
imagesc(abs(pc_result));
xlabel('快时间')
ylabel('慢时间')
title('原始回波脉压')
grid on
% axis([400,550,0,140,0,800])
% view(0,90)
% MTD=(fft(pc_result,[],1));
V_MTD = linspace(-PV,PV,M);
for i =1:L
    for i_v = 1:length(V_MTD)
        fd_t = 2*V_MTD(i_v)/lamda;
        mtd_ht = exp(1j*2*pi*fd_t*Tr*(0:M-1)).';
        road = pc_result(:,i);
        MTD(i_v,i) = sum(mtd_ht.*road);
    end
end

toc
[Range_MTD, V_MTD]=meshgrid(1:L, linspace(-PV,PV,M));
figure()
mesh(Range_MTD, V_MTD,(abs(MTD)));
title('未补偿相参积累')
xlabel('距离')
ylabel('多普勒单元')
figure
plot((abs(MTD(:,66))))
% axis([500,800,0,64])
%%
V = -600:5:600; %搜索的速度
RFT = zeros(length(V), L);
for i = 1:L  %%初始在的距离单元
    i/L
    for i_v = 1:length(V)
        fd_t = 2*V(i_v)/lamda; %
        rft_ht = exp(1j*2*pi*fd_t*Tr*(0:M-1));
        %% 找到目标运动的斜线
        r = (i - round( V(i_v)* Tr*(0:M-1)/delt_R));
        if isempty(find((r<1) )) && isempty(find((r>L) ))  
            for m =1:M
                road(m) = pc_result(m,r(m));
            end
            RFT(i_v, i) = sum(road.*rft_ht);
        end
    end
end
[Range,VV ] = meshgrid(1:L,V);
figure
mesh(Range, VV, abs(RFT))
xlabel('初始距离')
ylabel('速度')
% imagesc(1:L,V, abs(RFT))
