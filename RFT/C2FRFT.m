close all
clear all
clc
%����RFT�������Ѷ�Ŀ�����⣬Ŀ���ٶȲ���һ���ٶ�������
% f0=1e9;
fc=1000e6;%��Ƶ%2000
B=1e6;%����%4
Tao=127e-6;%����
Fs=B;%����Ƶ��
t=-Tao/2:1/Fs:Tao/2-1/Fs;%����ʱ��
mu=B/Tao;%��Ƶ��
C=3e8;
R_max=C*Tao/2;
delt_R=C/(2*Fs);%%�������뵥Ԫ
R_start1=87*delt_R;%Ŀ��1��ʼλ��%170
R_start2=27*delt_R;%Ŀ��2��ʼλ��%60
R_start3=47*delt_R;%Ŀ��3��ʼλ��lamda=C/fc;%80
lamda=C/fc;
PRF=1000;
Tr=1/PRF;
Vr_start1=570;%1��ʼ�ٶ�
Vr_start2=1060;%2��ʼ�ٶ�
Vr_start3=2050; %3��ʼ�ٶ�
Vb=lamda/(2*Tr);%%��һä��
pusle_num=1024;%������
PV=PRF*lamda/4;
L=length(t);
M=pusle_num;
    %%
% %     ��������
    echo=zeros(M,L);%�ز�
    echo1=zeros(M,L);%1�ز�
    echo2=zeros(M,L);%2�ز�
    echo_fft=zeros(M,L);%Ƶ��ز�
    pc_result_fft=zeros(M,L);%��ѹƵ���ź�
    pc_result=zeros(M,L);%��ѹʱ���ź�
%%
%%%ÿ������������ٶȡ����ٶ�
for i=1:pusle_num
    Vr1(i)=Vr_start1;
    fd1(i)=2*Vr1(i)/lamda;
    delt_t1(i)=2*(R_start1+Vr1(i)*Tr*(i-1))/C;%�ز��ӳ�
    
    Vr2(i)=Vr_start2;
    fd2(i)=2*Vr2(i)/lamda;
    delt_t2(i)=2*(R_start2+Vr2(i)*Tr*(i-1))/C;%�ز��ӳ�
    
    Vr3(i)=Vr_start3;
    fd3(i)=2*Vr3(i)/lamda;
    delt_t3(i)=2*(R_start3+Vr3(i)*Tr*(i-1))/C;%�ز��ӳ�
end
%%
%%��ѹϵ��
    ht_t=exp(-1j*2*pi*(mu/2*(t).^2));
    ht=conj(fliplr(ht_t));
    ht_fft=fftshift(fft(ht));%fftshift
%%
%%��ѹ
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
    echo(i,:)=awgn(echo(i,:),SNR);%%������
    echo_fft(i,:)=fftshift(fft(echo(i,:)));%fftshift
    pc_result(i,:)=(ifft(echo_fft(i,:).*ht_fft));%ifftshift
    pc_result_fft(i,:)=fftshift(fft(pc_result(i,:)));%fftshift%��ʱ��FFT
end
%%
% % ���ز�ͼ
figure()
mesh(real(echo))
view(0,90)
% % ��ѹ��ԭʼ�ز�
figure(), 
% mesh(abs(pc_result))
imagesc (abs(pc_result))
xlabel('���뵥Ԫ')
ylabel('������')
view(0,90)
title('ԭʼ�ز���ѹ')
grid on
%%
% % RFT_CZT��
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
Np=3;%�������
Pfa=10^-2; % �ּ���龯����
tic
count=1;
mubiao_num=1;%Ŀ����
mubiao=[1;0;0;0;0;0];%Ŀ���¼,��һ��ΪĿ�������ڶ���ΪĿ���ٶȵ�Ԫ���кţ�
%������ΪĿ������������ģ�����ӣ�������Ϊ���ѵ���Ŀ��������
%������Ϊ��Ծ�����ǵڼ���ģ������%��������Ŀ�����ھ��뵥Ԫ���к�
for Nk=1:Np:length(K)
    for i=1:L
    %%���Բ������Ӳ���������ģ��
    BU=exp(-1j*2*pi*K(Nk)*i.*m*lamda*B./(L*C));%%��������
    pc_result_fft_buchang= pc_result_fft(:,i).*BU;
    Sr(:,i)=czt((pc_result_fft_buchang),num_sou,Cn(i));
    end
    sp(:,:,count)=(ifft((Sr),[],2));%ifftshift
    sp_yuanshi(:,:,count)=sp(:,:,count);%����һ��ԭʼ��
%     %%��һ�μ�� global_cfar
    cankao_danyuan=abs(sp(:,:,count));
    Mean_cankao=(mean(mean(abs(cankao_danyuan))));%%�ο���Ԫ��ֵ
    alpha=(M*L)*(Pfa^(-1/(M*L))-1);%��alpha
    T=alpha*Mean_cankao;%%������;
        bijiao_v=abs(sp(:,:,count))>T;
        sp(:,:,count)=sp(:,:,count).*bijiao_v;
        if sum(sum(sp))~=0
        %Ѱ��Ǳ��Ŀ��
        max_sp_count=max(max(abs(sp(:,:,count))));
        [hang_mubiao,lie_mubiao]=find(abs(sp(:,:,count))==max_sp_count);
        if count==1
           mubiao(2,1)=hang_mubiao;
           mubiao(3,1)=K(Nk);
           mubiao(4,1)=max_sp_count;
           mubiao(5,1)=count;
           mubiao(6,1)=lie_mubiao;
        else
            if abs(mubiao(2,mubiao_num)-hang_mubiao)<30%�ж�ΪͬһĿ��ĸ���&&abs(mubiao(6,mubiao_num)-lie_mubiao)<20 
               if mubiao(4,mubiao_num)< max(max(abs(sp(:,:,count)))) % Ŀ����������֮ǰ�Ǹ���
                  mubiao(2,mubiao_num)=hang_mubiao;
                  mubiao(3,mubiao_num)=K(Nk); 
                  mubiao(4,mubiao_num)=max(max(abs(sp(:,:,count))));
                  mubiao(5,mubiao_num)=count;
                  mubiao(6,mubiao_num)=lie_mubiao;
               end 
            else %����ͬһ������,�������������1����һĿ�긱�꣬2������,���ҵ�Ŀ��
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
        kj(count)=K(Nk);%%%%�ɼ�⵽��ģ������
        count=count+1;
        end
end%%����ģ������RFT
mubiao_real=[];%��ʵĿ��
[m_hang,m_lie]=size(mubiao);%Ǳ��Ŀ��
[sp_hang,sp_lei,sp_ye]=size(sp);
for i=1:m_lie
    for j=1:sp_ye
        mubiao1(i,j)=abs(sp(mubiao(2,i),mubiao(6,i),j));%һ��һ��Ŀ���������仯����
    end
end
for i=1:m_lie
    figure
    stem(mubiao1(i,:))
end


% %%�޳���Ŀ�꣬�޳�����1��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %step1,���ݷ��������Ŀ��
bijiao2=mubiao(4,:)>0.5*M*B*Tao;%���ݷ������������Ŀ��
[~,tichu]=find(bijiao2==0);
mubiao_real=mubiao;
mubiao_real(:,tichu)=[];


% %%�޳���Ŀ�꣬�޳�����2����ƽ�����ȷ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% for i=1:mubiao_num%�ּ���Ŀ��������û�޳�
%     index_i=find(sp_ALL(:,mubiao(6,i))~=0&sp_ALL(:,mubiao(6,i))~=mubiao(4,i))%�ҵ���Ŀ����뵥Ԫ����Ϊ���ֵ&sp_ALL(:,mubiao(6,i))~=mubiao(4,i)
%     if length(index_i)==1
%        mubiao_real=[mubiao_real,mubiao(:,i)];
%     else
%         if mubiao(4,i)>1.5*mean(abs(sp_ALL(index_i,mubiao(6,i))))%�����Ŀ������ֵ���������ھ��뵥Ԫ������ֵ������
%             mubiao_real=[mubiao_real,mubiao(:,i)];
%         end
%     end
%     
% end
% %%�޳���Ŀ�꣬�޳�����3��������Ŀ���ƽ�����ȷ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:mubiao_num%�ּ���Ŀ��������û�޳�
%     [a_peak,index_peak]=findpeaks(abs(sp_ALL(:,mubiao(6,i))));%�ҵ���ֵ
%     index_it=find(abs(sp_ALL(index_peak,mubiao(6,i)))==max(abs(sp_ALL(index_peak,mubiao(6,i)))));%�ҵ����ǿ���������������ķ�ֵ
%     max_sp_ALL_i=abs(sp_ALL(index_peak(index_it),mubiao(6,i)));
%     if length(index_peak)~=1%���ֻ��һ����������Ϊ��Ŀ������
%         index_peak(index_it)=[];
%         if max_sp_ALL_i>2*mean(abs(sp_ALL(index_peak,mubiao(6,i))));
%             mubiao_real=[mubiao_real,mubiao(:,i)];
%         end
%     else
%         mubiao_real=[mubiao_real,mubiao(:,i)];
%     end
% end
%%%%�޳������������ݸ���仯����%%%%%%%%%%%%%%%%%%%%
% flag=1;%��ʵĿ��ı�ʶ��1Ϊ��ʵ��0Ϊ��
% for i=1:mubiao_num%�ּ���Ŀ��������û�޳�
%     [a_peak,index_peak]=findpeaks(abs(sp_ALL(:,mubiao(6,i))));%�ҵ���ֵ
%     [~,index_in_peak]=max(a_peak);%����ֵ�ڷ�ֵ�е�λ�á�
%     if index_in_peak==length(a_peak)%��ֵ�����һ��λ��
%         if index_in_peak==1%�ܹ���һ����ֵ
%            flag=1;%Ϊ��ʵĿ�� 
%         else%%��ֵ�����һ��λ���Ҳ���1����ֵ
%             for i_peak=1:length(a_peak)-1%�ӵ�һ����ʼ�Ƚ�
%                 if a_peak(i_peak)>1.1*a_peak(i_peak+1)
%                     flag=0;%ǰ��ķ�ֵ�Ⱥ����ֵ��Ϊ��Ŀ��
%                 end
%             end
%         end 
%     elseif index_in_peak==1%%�ڵ�һ��λ��
%         if index_in_peak==length(a_peak)%%��һ��
%             flag=1;%Ϊ��ʵĿ�� 
%         else
%             for i_peak=1:length(a_peak)-1%�ӵ�һ����ʼ�Ƚ�
%                 if a_peak(i_peak)<1.1*a_peak(i_peak+1)
%                     flag=0;%ǰ��ķ�ֵ�Ⱥ����ֵСΪ��Ŀ��
%                 end
%             end
%         end
%     else%%ʣ��������ǣ�������3����
%         if abs(a_peak(index_in_peak)-mean(a_peak))/mean(a_peak)<0.2%��������ƽ������һ�� 
%             flag=0;
%         end
%     end
%     if flag
%        mubiao_real=[mubiao_real,mubiao(:,i)];%Ϊ��ʵĿ�� 
%     end
%     flag=1;
% end
[M_mubiao,N_mubiao]=size(mubiao_real);

%�Ȱ�ѡ�����Ľ���ѡ��
sp_final=zeros(M,L);%û��cfar��ʱ��
sp_cfar=zeros(M,L);%����cfar�Ժ�
for j=1:N_mubiao%��ʵĿ����
    bijiao_t=sp_final>sp(:,:,mubiao_real(5,j));
    sp_final=sp_final.*bijiao_t+sp(:,:,mubiao_real(5,j)).*(ones(M,L)-bijiao_t);
end
%����Ŀ����ʵ���ٶ�
%����ѡ����
for j=1:N_mubiao%��ʵĿ����
    K_j=[mubiao_real(3,j)-1,mubiao_real(3,j)+1];%����Ŀ��ʣ��ģ������
    for NK_j=1:length(K_j)
        for i=1:L
        %%���Բ������Ӳ���������ģ��
        BU=exp(-1j*2*pi*K_j(NK_j)*i.*m*lamda*B./(L*C));%%��������
        pc_result_fft_buchang= pc_result_fft(:,i).*BU;
        Sr(:,i)=czt((pc_result_fft_buchang),num_sou,Cn(i));
        end
        sp1=(ifft((Sr),[],2));%ifftshift
        bijiao_t=sp_final>sp1;
        sp_final=sp_final.*bijiao_t+sp1.*(ones(M,L)-bijiao_t);
    end
end 
toc


%�������
max_sp=max(abs(sp_final),[],2);
max_sp=sort(max_sp,'descend');
sp_cfar=sp_final.*(abs(sp_final)>=max_sp(N_mubiao));
% Mean_cankao_cfar=(mean(mean(abs(sp_final))));%%�ο���Ԫ��ֵ
% alpha=(M*L)*(Pfa^(-1/(M*L))-1);%��alpha
% T=alpha*Mean_cankao_cfar;%%������;
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
%�ּ��
% figure()
% mesh(abs(sp_yuanshi(:,:,10)))
% figure()
% mesh(abs(sp(:,:,10)))
% sp_final=awgn(sp_cfar,-25,'measured');

