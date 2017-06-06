close all
clear all
clc
% % %��RFT����,������MTD������Ѷ�����ȥƥ�䣬���ǻ�������������԰�
% f0=1e9;
%%��������
fc=100e6;%��Ƶ
B=4e6;%����
Tao=128e-6;%����
Fs=1*B;%����Ƶ��
Ts=1/Fs;
t=-Tao/2:1/Fs:Tao/2-1/Fs;%����ʱ��
mu=B/Tao;%��Ƶ��
C=3e8;
% R0=80e3;%̽���λ��.
delt_R=C/(2*Fs)%%�������뵥Ԫ
R_start1=100*delt_R;%203
R_start2=253*delt_R;
R_start3=359*delt_R;
lamda=C/fc;
delt_R=C/(2*Fs);%%�������뵥Ԫ
PRF=500;
Tr=1/PRF;
Vr_start1=300;%��ʼ�ٶ�
Vr_start2=1200;%��ʼ�ٶ�
Vr_start3=-2000;%��ʼ�ٶ�
a=0;
a1=50;%50
a2=20;%-50
a3=90;%90
pusle_num=256;%������
PV=PRF*lamda/4;
L=length(t);
M=pusle_num;
PSLR1=M*lamda/2/delt_R
    %��������
    echo=zeros(L,M);%�ز�
    echo_fft=zeros(L,M);%Ƶ��ز�
    pc_result_fft=zeros(L,M);%��ѹƵ���ź�
    pc_result=zeros(L,M);%��ѹʱ���ź�
for i=1:pusle_num
    Vr1(i)=Vr_start1+a1*Tr*(i-1);
    fd1(i)=2*Vr1(i)/lamda;
    delt_t1(i)=2*(R_start1+Vr1(i)*Tr*(i-1))/C;%�ز��ӳ�
        Vr2(i)=Vr_start2+a2*Tr*(i-1);
    fd2(i)=2*Vr2(i)/lamda;
    delt_t2(i)=2*(R_start2+Vr2(i)*Tr*(i-1))/C;%�ز��ӳ�
        Vr3(i)=Vr_start3+a3*Tr*(i-1);
    fd3(i)=2*Vr3(i)/lamda;
    delt_t3(i)=2*(R_start3+Vr3(i)*Tr*(i-1))/C;%�ز��ӳ�
end
tic
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2)).';
    ht=conj((ht_t));%fliplr
    ht_fft=fft(ht);
% % ��ѹ
A1=1;
A2=1.1;
A3=1.2;
SNR=-20;
for i=1:pusle_num
   echo1(:,i)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
   echo2(:,i)=A2*exp(-1j*2*pi*(mu/2*(t+delt_t2(i)).^2)+-1j*2*pi*(fc)*(delt_t2(i)));
   echo3(:,i)=A3*exp(-1j*2*pi*(mu/2*(t+delt_t3(i)).^2)+-1j*2*pi*(fc)*(delt_t3(i)));
   echo(:,i)=echo1(:,i)+echo2(:,i)+echo3(:,i);
   echo(:,i)=awgn(echo(:,i),SNR);%%������
   echo_fft(:,i)=(fft(echo(:,i)));
   pc_result(:,i)=ifft((echo_fft(:,i).*ht_fft));
   pc_result_fft(:,i)=(fft(pc_result(:,i)));
% %    ��һ��
%    pc_result(i,:)=pc_result(i,:)./max(abs(pc_result(i,:)));
end
toc
%%����ѹ�������ȫ�̾���
R0=74e3;%��ʼ74km
% R1=80e3;%��ֹ80km
% PC_All=zeros(round((R1-R0)/delt_R),256);
% Pc=importdata('d:\Pc.txt');
% Gv_cpu=importdata('d:\Gv_cpu.txt');
% figure
% mesh(abs(Pc))
% aaaa=abs(pc_result)';
figure(), 
% mesh(abs(pc_result))
rang=(R0-L/2*delt_R:delt_R:R0+L/2*delt_R-delt_R)/1000;
imagesc((0:M-1)*Tr,rang,abs(pc_result))
xlabel('��λ���ʱ��(s)');
ylabel('��������(Km)');
% imagesc(abs(pc_result))
% xlabel('������')
% ylabel('���뵥Ԫ')
% title('ԭʼ�ز���ѹ')
% grid on
%%MTD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_V=lamda/(2*M*Tr);
vb=lamda/(2*Tr);%��һä��
tic
MTD=(fft(pc_result,[],2));
toc
% figure()
% mesh(abs(MTD));
% xlabel('')
V_MTD=-vb/2:delta_V:vb/2-delta_V;
rang_MTD=rang;
[X_MTD,Y_MTD]=meshgrid(V_MTD,rang_MTD);
figure
mesh(X_MTD,Y_MTD,(abs(MTD)));
xlabel('�ٶ�(m/s)')
ylabel('����(Km)')
zlabel('����')
zlabel('��һ������')
% ylabel('Searching range (Km)','FontSize',20);
% xlabel('Radial veocity (m/s)','FontSize',20);
% zlabel('Normalized RFT outputs','FontSize',20);
% axis([-40,40,64,84,0,1])
% RFT��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=round(1024/pusle_num);%%���ä������
V=-4*vb:delta_V:4*vb-delta_V;%%�ٶ�����
num_sou=length(V);
Gv=zeros(num_sou,L);
s_t=0;%%�������ʱ��
% R00=find(max(abs(pc_result(1,:)))==abs(pc_result(1,:)));
a_t=[50,20,90];
tic
for index_a=1:length(a_t)%���ٶ�����
    for vi=1:num_sou%�ٶ�i=L/5:L/2    i=1:L
        disp([num2str((vi+num_sou*(index_a-1))/num_sou/length(a_t)*100),' %']);
        for i=1:L%��ʼ���뵥Ԫ
        indexM=round((0:M-1)*Tr*(-V(vi))/delt_R)+i;
        index_All=indexM+(0:M-1)*L;
            if index_All>1&index_All<M*L
                pc_result(index_All);
%             fd_t=2*V(vi)/lamda;%%�����ٶȵ�ƥ�������
            fd_t=2*(V(vi)+a_t(index_a)*(0:M-1)*Tr)/lamda;
%             fdd(vi)=fd_t;
            find0=find(indexM<1);
            boola(vi,i)=isempty(find0);
                if boola(vi,i)
                %s_t=sum((pc_result(index_All).*exp(1j*2*pi*fd_t.*(0:M-1)*Tr)));%abs
        %         real_Gv=real(pc_result(index_All)).*(cos(2*pi*fd_t.*(0:M-1)*Tr))-imag(pc_result(index_All)).*(sin(2*pi*fd_t.*(0:M-1)*Tr));
        %         imag_Gv=real(pc_result(index_All)).*(sin(2*pi*fd_t.*(0:M-1)*Tr))+imag(pc_result(index_All)).*(cos(2*pi*fd_t.*(0:M-1)*Tr));
        %         s_t=sum(real_Gv)+1j*sum(imag_Gv);
        %         Gv(vi,i)=s_t;
                    Gv(vi,i,index_a)=sum(pc_result(index_All).*exp(1j*2*pi*fd_t.*(0:M-1)*Tr));%s_t;
                    s_t=0;
                end   
            end 
        end
    end
end
Gv_max=zeros(num_sou,L);%ѡ����
for i=1:length(a_t)
      bool_gv=abs(Gv(:,:,i))>Gv_max;
      Gv_max=Gv_max.*(1-bool_gv)+abs(Gv(:,:,i)).*bool_gv;
%     figure()
%     mesh(Y,X,fftshift(abs(Gv(:,:,i))))
end
toc
% save data1008.mat
% load data915.mat
% load data1008.mat
%%%��ͼ
% figure(),
% surf(real(echo))
% view(0,90)



%%%MTD
% MTD_fft=zeros(length(t),pusle_num);
% MTD_bu_fft=zeros(length(t),pusle_num);
% MTD_fft=(fft(pc_result,[],1));
% figure()
% mesh((abs(MTD)));
% title('δ������λ���')
% xlabel('����')
% ylabel('�����յ�Ԫ')

figure()
[X,Y]=meshgrid(rang,V);
mesh(Y,X,(abs(Gv_max))*1.3)%fftshift
% mesh(Y,X,(abs(Gv(:,50:end))));
xlabel('�ٶ�(m/s)')
ylabel('����(Km)')
zlabel('����')
% zlabel('��һ������')


% title('RFTͼ')
% [Sudu,Juli]=find(abs(Gv)==max(max(abs(Gv))));
% index_max=find(abs(Gv)<max(max(abs(Gv))));
% Gv_max=zeros(size(Gv));
% Gv_max(index_max)=abs(Gv(index_max))*0.3;
% Gv_max(Sudu,Juli)=abs(Gv(Sudu,Juli));
% figure()
% plot(V,(abs(Gv(:,Juli))));%10*log
% aaaaaa=abs(Gv);
% figure
% plot(abs(pc_result(:,1)))
