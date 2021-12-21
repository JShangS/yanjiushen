close all
clear all
clc
fc=1e9;
% fc=0;%��Ƶ
B=2e6;%����
Tao=100e-6;%����
Fs=2*B;%����Ƶ��
t=-Tao/2:1/Fs:Tao/2-1/Fs;%����ʱ��
mu=B/Tao;%��Ƶ��
R_start=0.5e4;
C=3e8;
lamda=C/fc;
delt_R=C/(2*Fs);%%�������뵥Ԫ

PRF=1000;
Tr=1/PRF;
Vr_start=300;%��ʼ�ٶ�
a=0;
pusle_num=128;%������
PV=PRF*lamda/4;
L=length(t);
R_cell = L-round( R_start/delt_R);
Hz=linspace(0,Fs,L);
M=pusle_num;
    %%%��������
    echo=zeros(M,L);%�ز�
    echo_bu=zeros(M,L);%keystone�任�ûز�
    echo_fft=zeros(M,L);%Ƶ��ز�
    echo_bu_fft=zeros(M,L);%keystone�任��Ƶ��ز�
    pc_result_fft=zeros(M,L);%��ѹƵ���ź�
    pc_result=zeros(M,L);%��ѹʱ���ź�
    pc_result_bu=zeros(M,L);%keystone�任��Ƶʱ����ѹ
    pc_result_bu_t=zeros(M,L);%keystone�任��Ƶʱ����ѹ��ʱ����
    MTD=zeros(M,L);%δ��keystone��MTD
    MTD_bu_fft=zeros(M,L);%��keystone��MTD

for i=1:pusle_num
    Vr(i)=Vr_start+a*Tr*(i-1);
    fd(i)=2*Vr(i)/lamda;
    delt_t(i)=2*(R_start+Vr(i)*Tr*(i-1))/C;%�ز��ӳ�
end
% % ��ѹϵ��
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj(fliplr(ht_t));
    ht_fft=fftshift(fft(ht));%
figure()
plot(Hz,abs(ht_fft))
%%��ѹ
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
xlabel('��ʱ��')
ylabel('��ʱ��')
title('ԭʼ�ز���ѹ')
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
title('δ������λ���')
xlabel('����')
ylabel('�����յ�Ԫ')
figure
plot((abs(MTD(:,66))))
% axis([500,800,0,64])
%%
V = -600:5:600; %�������ٶ�
RFT = zeros(length(V), L);
for i = 1:L  %%��ʼ�ڵľ��뵥Ԫ
    i/L
    for i_v = 1:length(V)
        fd_t = 2*V(i_v)/lamda; %
        rft_ht = exp(1j*2*pi*fd_t*Tr*(0:M-1));
        %% �ҵ�Ŀ���˶���б��
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
xlabel('��ʼ����')
ylabel('�ٶ�')
% imagesc(1:L,V, abs(RFT))
