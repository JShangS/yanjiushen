clc
clear all
close all
%%matlab和CUDA混合编程的尝试
system('nvcc -c PC.cu -ccbin "D:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin"')
mex PC.cpp PC.obj  -lcudart -lcufft -L"C:\Program Files\NVIDIA GPU ...
Computing Toolkit\CUDA\v7.0\lib\x64"
%%开始设置参数
fc=2000e6;%载频
B=4e6;%带宽
Tao=50e-6;%脉宽
Fs=1*B;%采样频率
t=-Tao/2:1/Fs:Tao/2-1/Fs;%脉冲时间
mu=B/Tao;%条频率
C=3e8;
R_max=C*Tao/2;
R_start=0;%%0时对应最后一个距离单元，L-mod(R_start/delt_R,L)即为目标所在初始距离单元
lamda=C/fc;
delt_R=C/(2*Fs);%%采样距离单元
PRF=200;
Tr=1/PRF;
Vb=lamda/(2*Tr);%%第一盲速
Vr_start=25;%初始速度
PV=PRF*lamda/4;
L=length(t);
h_t=exp(-1j*2*pi*(mu/2*(t).^2));%%传递函数
M=1;
echo=h_t;%回波
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
[realpc_GPU,imagpc_GPU]=PC(realh,imagh,realecho,imagecho);%%GPU脉压
time2=toc
pc_GPU=realpc_GPU+1i*imagpc_GPU;
figure
plot(fftshift(abs(pc_GPU)))
