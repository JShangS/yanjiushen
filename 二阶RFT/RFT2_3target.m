close all
clear all
clc
% % %用RFT方法,类似于MTD，用最佳多普勒去匹配，但是会带来量化误差和旁瓣
% f0=1e9;
%%参数计算
fc=100e6;%载频
B=4e6;%带宽
Tao=128e-6;%脉宽
Fs=1*B;%采样频率
Ts=1/Fs;
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
mu=B/Tao;%条频率
C=3e8;
% R0=80e3;%探测的位置.
delt_R=C/(2*Fs)%%采样距离单元
R_start1=100*delt_R;%203
R_start2=253*delt_R;
R_start3=359*delt_R;
lamda=C/fc;
delt_R=C/(2*Fs);%%采样距离单元
PRF=500;
Tr=1/PRF;
Vr_start1=300;%初始速度
Vr_start2=1200;%初始速度
Vr_start3=-2000;%初始速度
a=0;
a1=50;%50
a2=20;%-50
a3=90;%90
pusle_num=256;%脉冲数
PV=PRF*lamda/4;
L=length(t);
M=pusle_num;
PSLR1=M*lamda/2/delt_R
    %变量声明
    echo=zeros(L,M);%回波
    echo_fft=zeros(L,M);%频域回波
    pc_result_fft=zeros(L,M);%脉压频域信号
    pc_result=zeros(L,M);%脉压时域信号
for i=1:pusle_num
    Vr1(i)=Vr_start1+a1*Tr*(i-1);
    fd1(i)=2*Vr1(i)/lamda;
    delt_t1(i)=2*(R_start1+Vr1(i)*Tr*(i-1))/C;%回拨延迟
        Vr2(i)=Vr_start2+a2*Tr*(i-1);
    fd2(i)=2*Vr2(i)/lamda;
    delt_t2(i)=2*(R_start2+Vr2(i)*Tr*(i-1))/C;%回拨延迟
        Vr3(i)=Vr_start3+a3*Tr*(i-1);
    fd3(i)=2*Vr3(i)/lamda;
    delt_t3(i)=2*(R_start3+Vr3(i)*Tr*(i-1))/C;%回拨延迟
end
tic
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2)).';
    ht=conj((ht_t));%fliplr
    ht_fft=fft(ht);
% % 脉压
A1=1;
A2=1.1;
A3=1.2;
SNR=-20;
for i=1:pusle_num
   echo1(:,i)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
   echo2(:,i)=A2*exp(-1j*2*pi*(mu/2*(t+delt_t2(i)).^2)+-1j*2*pi*(fc)*(delt_t2(i)));
   echo3(:,i)=A3*exp(-1j*2*pi*(mu/2*(t+delt_t3(i)).^2)+-1j*2*pi*(fc)*(delt_t3(i)));
   echo(:,i)=echo1(:,i)+echo2(:,i)+echo3(:,i);
   echo(:,i)=awgn(echo(:,i),SNR);%%加噪声
   echo_fft(:,i)=(fft(echo(:,i)));
   pc_result(:,i)=ifft((echo_fft(:,i).*ht_fft));
   pc_result_fft(:,i)=(fft(pc_result(:,i)));
% %    归一化
%    pc_result(i,:)=pc_result(i,:)./max(abs(pc_result(i,:)));
end
toc
%%把脉压结果加入全程距离
R0=74e3;%起始74km
% R1=80e3;%终止80km
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
xlabel('相参积累时间(s)');
ylabel('搜索距离(Km)');
% imagesc(abs(pc_result))
% xlabel('脉冲数')
% ylabel('距离单元')
% title('原始回波脉压')
% grid on
%%MTD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_V=lamda/(2*M*Tr);
vb=lamda/(2*Tr);%第一盲速
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
xlabel('速度(m/s)')
ylabel('距离(Km)')
zlabel('幅度')
zlabel('归一化幅度')
% ylabel('Searching range (Km)','FontSize',20);
% xlabel('Radial veocity (m/s)','FontSize',20);
% zlabel('Normalized RFT outputs','FontSize',20);
% axis([-40,40,64,84,0,1])
% RFT法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=round(1024/pusle_num);%%最大盲速因子
V=-4*vb:delta_V:4*vb-delta_V;%%速度搜索
num_sou=length(V);
Gv=zeros(num_sou,L);
s_t=0;%%求积分临时量
% R00=find(max(abs(pc_result(1,:)))==abs(pc_result(1,:)));
a_t=[50,20,90];
tic
for index_a=1:length(a_t)%加速度搜索
    for vi=1:num_sou%速度i=L/5:L/2    i=1:L
        disp([num2str((vi+num_sou*(index_a-1))/num_sou/length(a_t)*100),' %']);
        for i=1:L%初始距离单元
        indexM=round((0:M-1)*Tr*(-V(vi))/delt_R)+i;
        index_All=indexM+(0:M-1)*L;
            if index_All>1&index_All<M*L
                pc_result(index_All);
%             fd_t=2*V(vi)/lamda;%%搜索速度的匹配多普勒
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
Gv_max=zeros(num_sou,L);%选大处理
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
%%%画图
% figure(),
% surf(real(echo))
% view(0,90)



%%%MTD
% MTD_fft=zeros(length(t),pusle_num);
% MTD_bu_fft=zeros(length(t),pusle_num);
% MTD_fft=(fft(pc_result,[],1));
% figure()
% mesh((abs(MTD)));
% title('未补偿相参积累')
% xlabel('距离')
% ylabel('多普勒单元')

figure()
[X,Y]=meshgrid(rang,V);
mesh(Y,X,(abs(Gv_max))*1.3)%fftshift
% mesh(Y,X,(abs(Gv(:,50:end))));
xlabel('速度(m/s)')
ylabel('距离(Km)')
zlabel('幅度')
% zlabel('归一化幅度')


% title('RFT图')
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
