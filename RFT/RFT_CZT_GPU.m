close all
clear all
clc
%��RFT����,������MTD������Ѷ�����ȥƥ�䣬���ǻ�������������԰�
%%�Լ�д��Chirp-Zmatlab�е�һ��ֻ�ܾ���һ��CZT�����ܹ�������������б仯������
%%�ڴ˻����Ͻ��еĸĽ�,�ܲ��ܷ������������
% f0=1e9;
fc=100e6;%��Ƶ1000e6;
B=4e6;%����1e6;
Tao=128e-6;%����127e-6;
Fs=1*B;%����Ƶ��B;
t=-Tao/2:1/Fs:Tao/2-1/Fs;%����ʱ��
mu=B/Tao;%��Ƶ��
C=3e8;
R_max=C*Tao/2;
lamda=C/fc;
delt_R=C/(2*Fs);%%�������뵥Ԫ
R_start=100*delt_R;
PRF=500;
Tr=1/PRF;
Vr_start=1200;%��ʼ�ٶ�1700;
aa=0;
Vb=lamda/(2*Tr);%%��һä��
pusle_num=256;%������2048;
PV=PRF*lamda/4;
L=length(t);
M=pusle_num;
    %%
% %     ��������
    echo=zeros(M,L);%�ز�
    echo_fft=zeros(M,L);%Ƶ��ز�
    pc_result_fft=zeros(M,L);%��ѹƵ���ź�
    pc_result=zeros(M,L);%��ѹʱ���ź�
%%
%%%ÿ������������ٶȡ����ٶ�
for i=1:pusle_num
    Vr(i)=Vr_start+aa*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%�ز��ӳ�
end
%%
%%��ѹϵ��
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj((ht_t));
    ht_fft=(fft(ht));%fftshift 
%%
% %%%����
% noise=wgn(pusle_num,L,16)+1i*wgn(pusle_num,L,16);%��˹������wgn��m,n,p��,p����dBwΪ��λ�ĸ�˹����

%%
% % ��ѹ
SNR=10;
% K=-4:4;
for i=1:pusle_num
   echo(i,:)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));
%    echo(i,:)=awgn(echo(i,:),SNR);%%������
end
% tic
for i=1:pusle_num
   echo_fft(i,:)=(fft(echo(i,:)));%fftshift
   pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
   pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%��ʱ��FFT
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
% % ���ز�ͼ
% figure()
% surf(real(echo))
% view(0,90)
% % % ��ѹ��ԭʼ�ز�
% figure(), 
% % mesh(abs(pc_result))
% imagesc (abs(pc_result))
% xlabel('��ʱ��')
% ylabel('��ʱ��')
% view(0,90)
% title('ԭʼ�ز���ѹ')
% grid on
% juli_kuadu=find(max(abs(pc_result(1,:)))==abs(pc_result(1,:)))-find(max(abs(pc_result(pusle_num,:)))==abs(pc_result(pusle_num,:)))
% figure(), 
% mesh(abs(pc_result_fft))
%%
% % RFT_CZT��
fai=lamda*B/(L*C);
delta_V=lamda/(2*M*Tr);
%%ϵ������
% V=Vb-delta_V:-delta_V:0;%%�ٶ�����
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
% K=2;%ģ������
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
    BU=exp(-1j*2*pi*K(Nk).*Bu_mL.*lamda*B./(L*C));%%��������
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
    fy_GPU = fft(y_GPU, nfft);%Ĭ������FFT
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
    trans_band=hang*lie*16/1024/1024/time_trans%�������Mb/s   
    mesh((abs(sp)))