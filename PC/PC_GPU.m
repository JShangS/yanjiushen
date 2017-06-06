clc
clear all
close all
%%matlab��CUDA��ϱ�̵ĳ���
system('nvcc -c PC.cu -ccbin "D:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin"')
mex PC.cpp PC.obj  -lcudart -lcufft -L"C:\Program Files\NVIDIA GPU ...
Computing Toolkit\CUDA\v7.0\lib\x64"
%%��ʼ���ò���
fc=2000e6;%��Ƶ
B=4e6;%����
Tao=50e-6;%����
Fs=1*B;%����Ƶ��
t=-Tao/2:1/Fs:Tao/2-1/Fs;%����ʱ��
mu=B/Tao;%��Ƶ��
C=3e8;
R_max=C*Tao/2;
R_start=0;%%0ʱ��Ӧ���һ�����뵥Ԫ��L-mod(R_start/delt_R,L)��ΪĿ�����ڳ�ʼ���뵥Ԫ
lamda=C/fc;
delt_R=C/(2*Fs);%%�������뵥Ԫ
PRF=200;
Tr=1/PRF;
Vb=lamda/(2*Tr);%%��һä��
Vr_start=25;%��ʼ�ٶ�
PV=PRF*lamda/4;
L=length(t);
h_t=exp(-1j*2*pi*(mu/2*(t).^2));%%���ݺ���
M=1;
echo=h_t;%�ز�
tic
echo_fft=fft(echo);
h=conj(fliplr(h_t));
h_fft=(fft(h));
for i=1:M
    pc(i,:)=ifftshift(ifft(echo_fft.*h_fft));%ifftshift
end
time1=toc
echo=echo.';
t_echo=echo;
% for i=1:127
%     echo=[echo,t_echo];
% end
plot((abs(pc)));

realh=real(h);imagh=imag(h);realecho=real(echo);imagecho=imag(echo);
tic
[realpc_GPU,imagpc_GPU]=PC(realh,imagh,realecho,imagecho);%%GPU��ѹ
time2=toc
pc_GPU=realpc_GPU+1i*imagpc_GPU;
figure
plot(fftshift(abs(pc_GPU)))
