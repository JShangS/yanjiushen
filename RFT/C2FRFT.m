close all
clear all
clc
%快速RFT跳区间搜多目标问题，目标速度不在一个速度区间内
% f0=1e9;
fc=1000e6;%载频%2000
B=1e6;%带宽%4
Tao=127e-6;%脉宽
Fs=B;%采样频率
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
mu=B/Tao;%条频率
C=3e8;
R_max=C*Tao/2;
delt_R=C/(2*Fs);%%采样距离单元
R_start1=87*delt_R;%目标1初始位置%170
R_start2=27*delt_R;%目标2初始位置%60
R_start3=47*delt_R;%目标3初始位置lamda=C/fc;%80
lamda=C/fc;
PRF=1000;
Tr=1/PRF;
Vr_start1=570;%1初始速度
Vr_start2=1060;%2初始速度
Vr_start3=2050; %3初始速度
Vb=lamda/(2*Tr);%%第一盲速
pusle_num=1024;%脉冲数
PV=PRF*lamda/4;
L=length(t);
M=pusle_num;
    %%
% %     变量声明
    echo=zeros(M,L);%回波
    echo1=zeros(M,L);%1回波
    echo2=zeros(M,L);%2回波
    echo_fft=zeros(M,L);%频域回波
    pc_result_fft=zeros(M,L);%脉压频域信号
    pc_result=zeros(M,L);%脉压时域信号
%%
%%%每个脉冲间隔后的速度、加速度
for i=1:pusle_num
    Vr1(i)=Vr_start1;
    fd1(i)=2*Vr1(i)/lamda;
    delt_t1(i)=2*(R_start1+Vr1(i)*Tr*(i-1))/C;%回拨延迟
    
    Vr2(i)=Vr_start2;
    fd2(i)=2*Vr2(i)/lamda;
    delt_t2(i)=2*(R_start2+Vr2(i)*Tr*(i-1))/C;%回拨延迟
    
    Vr3(i)=Vr_start3;
    fd3(i)=2*Vr3(i)/lamda;
    delt_t3(i)=2*(R_start3+Vr3(i)*Tr*(i-1))/C;%回拨延迟
end
%%
%%脉压系数
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj(fliplr(ht_t));
    ht_fft=fftshift(fft(ht));%fftshift
%%
%%脉压
SNR=-10;
A1=1; %%1.8%1.45
A2=1;%%1.25%1.28
A3=1;%%1.25%1.05
%K=40:60;
for i=1:pusle_num
    echo1(i,:)=A1*exp(-1j*2*pi*(mu/2*(t+delt_t1(i)).^2)+-1j*2*pi*(fc)*(delt_t1(i)));
    echo2(i,:)=A2*exp(-1j*2*pi*(mu/2*(t+delt_t2(i)).^2)+-1j*2*pi*(fc)*(delt_t2(i)));
    echo3(i,:)=A3*exp(-1j*2*pi*(mu/2*(t+delt_t3(i)).^2)+-1j*2*pi*(fc)*(delt_t3(i)));
    echo(i,:)=echo1(i,:)+echo2(i,:)+echo3(i,:);
    echo(i,:)=awgn(echo(i,:),SNR);%%加噪声
    echo_fft(i,:)=fftshift(fft(echo(i,:)));%fftshift
    pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
    pc_result_fft(i,:)=fftshift(fft(pc_result(i,:)));%fftshift%快时间FFT
end
%%
% % 画回波图
figure()
mesh(real(echo))
view(0,90)
% % 脉压后原始回波
figure(), 
% mesh(abs(pc_result))
imagesc (abs(pc_result))
xlabel('距离单元')
ylabel('脉冲数')
view(0,90)
title('原始回波脉压')
grid on
%%
% % RFT_CZT法
fai=lamda*B/(L*C);
delta_V=lamda/(2*M*Tr);
num_sou=M;
Sr=zeros(num_sou,L);
m=1:M;
m=m.';
a=1-fai.*(1:L);
Cn=exp(-1j*2*pi*a./num_sou);
sp_ALL=[];
sp_ALL_yuanshi=[];
P=30;
K=-P:P;
% K=[4,7,13];
Np=3;%搜索间隔
Pfa=10^-2; % 粗检测虚警概率
tic
count=1;
mubiao_num=1;%目标数
mubiao=[1;0;0;0;0;0];%目标记录,第一行为目标数，第二行为目标速度单元即行号，
%第三行为目标最大幅度所属模糊因子，第四行为所搜到的目标最大幅度
%第五行为跳跃搜索是第几个模糊因子%第六航是目标所在距离单元即列号
for Nk=1:Np:length(K)
    for i=1:L
    %%乘以补偿因子补偿多普勒模糊
    BU=exp(-1j*2*pi*K(Nk)*i.*m*lamda*B./(L*C));%%补偿因子
    pc_result_fft_buchang= pc_result_fft(:,i).*BU;
    Sr(:,i)=czt((pc_result_fft_buchang),num_sou,Cn(i));
    end
    sp(:,:,count)=(ifft((Sr),[],2));%ifftshift
    sp_yuanshi(:,:,count)=sp(:,:,count);%保留一份原始的
%     %%第一次检测 global_cfar
    cankao_danyuan=abs(sp(:,:,count));
    Mean_cankao=(mean(mean(abs(cankao_danyuan))));%%参考单元均值
    alpha=(M*L)*(Pfa^(-1/(M*L))-1);%求alpha
    T=alpha*Mean_cankao;%%求门限;
        bijiao_v=abs(sp(:,:,count))>T;
        sp(:,:,count)=sp(:,:,count).*bijiao_v;
        if sum(sum(sp))~=0
        %寻找潜在目标
        max_sp_count=max(max(abs(sp(:,:,count))));
        [hang_mubiao,lie_mubiao]=find(abs(sp(:,:,count))==max_sp_count);
        if count==1
           mubiao(2,1)=hang_mubiao;
           mubiao(3,1)=K(Nk);
           mubiao(4,1)=max_sp_count;
           mubiao(5,1)=count;
           mubiao(6,1)=lie_mubiao;
        else
            if abs(mubiao(2,mubiao_num)-hang_mubiao)<30%判断为同一目标的副瓣&&abs(mubiao(6,mubiao_num)-lie_mubiao)<20 
               if mubiao(4,mubiao_num)< max(max(abs(sp(:,:,count)))) % 目标这个副瓣比之前那个高
                  mubiao(2,mubiao_num)=hang_mubiao;
                  mubiao(3,mubiao_num)=K(Nk); 
                  mubiao(4,mubiao_num)=max(max(abs(sp(:,:,count))));
                  mubiao(5,mubiao_num)=count;
                  mubiao(6,mubiao_num)=lie_mubiao;
               end 
            else %不是同一个副瓣,则有两种情况，1是另一目标副瓣，2是噪声,暂且当目标
                mubiao_num=mubiao_num+1;
                mubiao=[mubiao,[mubiao_num;0;0;0;0;0]];
                mubiao(2,mubiao_num)=hang_mubiao;
                mubiao(3,mubiao_num)=K(Nk); 
                mubiao(4,mubiao_num)=max(max(abs(sp(:,:,count))));
                mubiao(5,mubiao_num)=count;
                mubiao(6,mubiao_num)=lie_mubiao;
            end
        end
        sp_ALL=cat(1,sp_ALL,sp(:,:,count));
        sp_ALL_yuanshi=cat(1,sp_ALL_yuanshi,sp_yuanshi(:,:,count));
        kj(count)=K(Nk);%%%%可检测到的模糊因子
        count=count+1;
        end
end%%搜索模糊因子RFT
mubiao_real=[];%真实目标
[m_hang,m_lie]=size(mubiao);%潜在目标
[sp_hang,sp_lei,sp_ye]=size(sp);
for i=1:m_lie
    for j=1:sp_ye
        mubiao1(i,j)=abs(sp(mubiao(2,i),mubiao(6,i),j));%一列一个目标的主副瓣变化规律
    end
end
for i=1:m_lie
    figure
    stem(mubiao1(i,:))
end


% %%剔除假目标，剔除规则1，看幅度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %step1,根据幅度提出假目标
bijiao2=mubiao(4,:)>0.5*M*B*Tao;%根据幅度提出噪声假目标
[~,tichu]=find(bijiao2==0);
mubiao_real=mubiao;
mubiao_real(:,tichu)=[];


% %%剔除假目标，剔除规则2，看平均幅度幅度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% for i=1:mubiao_num%粗检后的目标数，还没剔除
%     index_i=find(sp_ALL(:,mubiao(6,i))~=0&sp_ALL(:,mubiao(6,i))~=mubiao(4,i))%找到该目标距离单元处不为零的值&sp_ALL(:,mubiao(6,i))~=mubiao(4,i)
%     if length(index_i)==1
%        mubiao_real=[mubiao_real,mubiao(:,i)];
%     else
%         if mubiao(4,i)>1.5*mean(abs(sp_ALL(index_i,mubiao(6,i))))%如果该目标最大峰值大于其所在距离单元其他峰值的两倍
%             mubiao_real=[mubiao_real,mubiao(:,i)];
%         end
%     end
%     
% end
% %%剔除假目标，剔除规则3，看现有目标的平均幅度幅度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:mubiao_num%粗检后的目标数，还没剔除
%     [a_peak,index_peak]=findpeaks(abs(sp_ALL(:,mubiao(6,i))));%找到峰值
%     index_it=find(abs(sp_ALL(index_peak,mubiao(6,i)))==max(abs(sp_ALL(index_peak,mubiao(6,i)))));%找到不是可能主瓣所在区间的峰值
%     max_sp_ALL_i=abs(sp_ALL(index_peak(index_it),mubiao(6,i)));
%     if length(index_peak)~=1%如果只有一个主瓣则认为是目标主瓣
%         index_peak(index_it)=[];
%         if max_sp_ALL_i>2*mean(abs(sp_ALL(index_peak,mubiao(6,i))));
%             mubiao_real=[mubiao_real,mubiao(:,i)];
%         end
%     else
%         mubiao_real=[mubiao_real,mubiao(:,i)];
%     end
% end
%%%%剔除规则三，根据副瓣变化规律%%%%%%%%%%%%%%%%%%%%
% flag=1;%真实目标的标识，1为真实，0为假
% for i=1:mubiao_num%粗检后的目标数，还没剔除
%     [a_peak,index_peak]=findpeaks(abs(sp_ALL(:,mubiao(6,i))));%找到峰值
%     [~,index_in_peak]=max(a_peak);%最大峰值在峰值中的位置。
%     if index_in_peak==length(a_peak)%峰值在最后一个位置
%         if index_in_peak==1%总共就一个峰值
%            flag=1;%为真实目标 
%         else%%峰值在最后一个位置且不是1个峰值
%             for i_peak=1:length(a_peak)-1%从第一个开始比较
%                 if a_peak(i_peak)>1.1*a_peak(i_peak+1)
%                     flag=0;%前面的峰值比后面峰值大为假目标
%                 end
%             end
%         end 
%     elseif index_in_peak==1%%在第一个位置
%         if index_in_peak==length(a_peak)%%就一个
%             flag=1;%为真实目标 
%         else
%             for i_peak=1:length(a_peak)-1%从第一个开始比较
%                 if a_peak(i_peak)<1.1*a_peak(i_peak+1)
%                     flag=0;%前面的峰值比后面峰值小为假目标
%                 end
%             end
%         end
%     else%%剩下情况就是，至少有3个瓣
%         if abs(a_peak(index_in_peak)-mean(a_peak))/mean(a_peak)<0.2%噪声，电平基本上一样 
%             flag=0;
%         end
%     end
%     if flag
%        mubiao_real=[mubiao_real,mubiao(:,i)];%为真实目标 
%     end
%     flag=1;
% end
[M_mubiao,N_mubiao]=size(mubiao_real);

%先把选出来的进行选大
sp_final=zeros(M,L);%没做cfar的时候
sp_cfar=zeros(M,L);%做了cfar以后
for j=1:N_mubiao%真实目标数
    bijiao_t=sp_final>sp(:,:,mubiao_real(5,j));
    sp_final=sp_final.*bijiao_t+sp(:,:,mubiao_real(5,j)).*(ones(M,L)-bijiao_t);
end
%计算目标真实的速度
%进行选大处理
for j=1:N_mubiao%真实目标数
    K_j=[mubiao_real(3,j)-1,mubiao_real(3,j)+1];%计算目标剩余模糊因子
    for NK_j=1:length(K_j)
        for i=1:L
        %%乘以补偿因子补偿多普勒模糊
        BU=exp(-1j*2*pi*K_j(NK_j)*i.*m*lamda*B./(L*C));%%补偿因子
        pc_result_fft_buchang= pc_result_fft(:,i).*BU;
        Sr(:,i)=czt((pc_result_fft_buchang),num_sou,Cn(i));
        end
        sp1=(ifft((Sr),[],2));%ifftshift
        bijiao_t=sp_final>sp1;
        sp_final=sp_final.*bijiao_t+sp1.*(ones(M,L)-bijiao_t);
    end
end 
toc


%最后检测结果
max_sp=max(abs(sp_final),[],2);
max_sp=sort(max_sp,'descend');
sp_cfar=sp_final.*(abs(sp_final)>=max_sp(N_mubiao));
% Mean_cankao_cfar=(mean(mean(abs(sp_final))));%%参考单元均值
% alpha=(M*L)*(Pfa^(-1/(M*L))-1);%求alpha
% T=alpha*Mean_cankao_cfar;%%求门限;
% bijiao_cfar=abs(sp_final)>T;
% sp_cfar=sp_final.*bijiao_cfar;
figure
mesh(abs(sp_cfar))
% figure()
% mesh(abs(sp_ALL))
% view(90,0)
figure
mesh(abs(sp_final))
% figure
% mesh(abs(sp_cfar))
% P_label_t=(linspace(min(kj),max(kj),length(sp_ALL)));
% [P_label,R_label]=meshgrid(1:L,P_label_t);
% figure()
% mesh(P_label,R_label,(abs(sp_ALL)))
% % view(90,0)
% % grid off
% figure()
% plot(P_label_t,(abs(sp_ALL(:,mubiao(6,2)))),'k-')
% figure()
% plot(P_label_t,(abs(sp_ALL(:,97))),'k-')
% figure()
% plot(abs(sp_ALL(:,47)))
%粗检测
% figure()
% mesh(abs(sp_yuanshi(:,:,10)))
% figure()
% mesh(abs(sp(:,:,10)))
% sp_final=awgn(sp_cfar,-25,'measured');

