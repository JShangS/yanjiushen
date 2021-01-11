clear
clc
T =10e-6;  %LFM�ź�����10us
fc = 0;  %����Ƶ��
TBP = 150;  %ʱ������
K = TBP/T^2;  %���Ե�Ƶ��
B = K*T;  %����
oversamplerate=5;  %��������
fs = oversamplerate*B;  
Ts=1/fs;
t=(-T/2:Ts:T/2);  %������λ��
L = length(t);
t0 = 0;  %����ʱ��Ϊ0


sr = rectpuls((t-t0)/T).*exp(1i*pi*K*(t-t0).^2);  %����Ŀ��ز�
h = rectpuls(t/T).*exp(-1i*pi*K*(t).^2);  %ƥ���˲���
sout = conv(sr,h,'same');  %ʱ����  same��ȡ�м䲿��
figure()
plot(abs(sout))

Sout = 20*log10(abs(sout)/max(abs(sout)));  %��һ����ת����dB תdBʱ����ȡ20 ����ȡ10
% IRW = 0.886/(K*T) %����IRW
% Figure 1: �������Ե�Ƶ�źŵ�ƥ���˲�
% figure(1)
% subplot(2,2,1);
% plot(t/1e-6,real(sr));
% xlabel('t/us');ylabel('���� ');title('(a)ԭʼ�ź�ʵ�� ');
% 
% subplot(2,2,2);
% % plot(t/1e-6,Sout);
% plot(abs(sout));
% % axis([-0.4 0.4 -40 1.25]);
% xlabel('t/us');ylabel('���� /dB');title('(b) ��ѹ(����չ) ');
% 
% subplot(2,2,3);
% plot(t/1e-6,abs(sout));
% % axis([-1 1 -5 600]);
% xlabel('t/us');ylabel('����');title('(c) ѹ������ź�');
% 
% subplot(2,2,4);
% sout=real(sout);
% plot(t/1e-6,angle(sout))
% % axis([-1 1 -5 5]);
% xlabel('t/us');ylabel('����');title('(d) ѹ�����źŵ���λ������չ��')

% figure(2)
% 
% sssout=T*sinc(K*T*(t-t