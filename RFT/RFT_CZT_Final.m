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
Fs=4*B;%����Ƶ��B;
t=-Tao/2:1/Fs:Tao/2-1/Fs;%����ʱ��
mu=B/Tao;%��Ƶ��
C=3e8;
R_max=C*Tao/2;
lamda=C/fc;
delt_R=C/(2*Fs);%%�������뵥Ԫ
R_start=100*delt_R;
PRF=400;
Tr=1/PRF;
Vr_start=200;%��ʼ�ٶ�1700;
aa=0;
Vb=lamda/(2*Tr);%%��һä��
pusle_num=1024;%������2048;
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
for i=1:pusle_num
   echo(i,:)=exp(-1j*2*pi*(mu/2*(t+delt_t(i)).^2)+-1j*2*pi*(fc)*(delt_t(i)));
%    echo(i,:)=awgn(echo(i,:),SNR);%%������
end
tic
for i=1:pusle_num
   echo_fft(i,:)=(fft(echo(i,:)));%fftshift
   pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
   pc_result_fft(i,:)=(fft(pc_result(i,:)));%fftshift%��ʱ��FFT
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
% % ���ز�ͼ
% figure()
% surf(real(echo))
% view(0,90)
% % ��ѹ��ԭʼ�ز�
figure(), 
% mesh(abs(pc_result))
imagesc (abs(pc_result))
xlabel('��ʱ��')
ylabel('��ʱ��')
view(0,90)
title('ԭʼ�ز���ѹ')
grid on
% juli_kuadu=find(max(abs(pc_result(1,:)))==abs(pc_result(1,:)))-find(max(abs(pc_result(pusle_num,:)))==abs(pc_result(pusle_num,:)))
% figure(), 
% mesh(abs(pc_result_fft))
%%
% % RFT_CZT��
K=2;
fai=lamda*B/(L*C);
delta_V=lamda/(2*M*Tr);
%%ϵ������
V=(K(1)-0.5)*Vb:Vb/M:(K(end)+0.5)*Vb-Vb/M;%%�ٶ�����
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
tic
for Nk=1:length(K)
%     disp([num2str(Nk/length(K)*100),' %']);
    %%���Բ������Ӳ���������ģ��
    BU=exp(-1j*2*pi*K(Nk).*Bu_mL.*lamda*B./(L*C));%%��������
    x= pc_result_fft.*BU;
    nn = (0:(M-1))';
    y = x .* ww(M+nn,:);
    fy = fft(y, nfft);%Ĭ������FFT
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
% xlabel('����','FontSize',30)
% ylabel('�ٶ�','FontSize',30)
% zlabel('����','FontSize',30)
% title('����RFTͼ')
% [hang_all,lie_all]=find(abs(sp_ALL)==max(max(abs(sp_ALL))));
% figure()
% plot(V_all,abs(sp_ALL(:,lie_all)));
% xlabel('�ٶ�')
% ylabel('����')
% title('RFT�ٶ�����ͼ')
% figure()
% mesh((abs(MTD_fft)));
% title('δ������λ���')
% xlabel('����')
% ylabel('�����յ�Ԫ')
% zlabel('����')
ww_z=ww.';