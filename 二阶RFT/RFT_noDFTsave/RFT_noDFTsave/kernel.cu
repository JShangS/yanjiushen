//�������ȴ洢��DFTϵ����ͨ��ע�����ı�����DFT�洢��ʽ
//
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <cublas.h>
#include <cufft.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Windows.h>
#define PI 3.14159265358979323846
#define M 256//���������
///kernel ����
//���
__global__ void MulVector(cufftComplex *a, cufftComplex *b, cufftComplex *c,int size)
{
	int id_by = blockIdx.y;//���ڵĿ������,�ڼ�����
	int id_tx = threadIdx.x;//���ڿ���̺߳ź��ᣬ�ڼ������뵥Ԫ
	int id_ty = threadIdx.y;//���ڿ���̺߳����ᣬ���ڵڼ�������
	int threadPerBlock_x = blockDim.x;//���ڳߴ�ĺ��᳤�ȣ���L
	int threadPerBlock_y = blockDim.y;//���ڳߴ�����᳤�ȣ���N

	int index_real = (id_tx+threadPerBlock_x*id_ty)+id_by*threadPerBlock_x*threadPerBlock_y; //����������ֵ=�ڿ���+����Ŀ���
	int index_ht = id_tx;//���ݺ���������
	int index_echo = index_real;// �ز�����������
	if (index_real < size)
	{
		c[index_real].x = a[index_ht].x * b[index_echo].x - a[index_ht].y * b[index_echo].y;//blockDim.x;//index_ht;//
		c[index_real].y = a[index_ht].x * b[index_echo].y + a[index_ht].y * b[index_echo].x;//blockDim.y;//index_real;//
	}
}
//���
__global__ void ChuVector(cufftComplex *a, int size,int L)
{
    int id_by = blockIdx.y;//���ڵĿ������,�ڼ�����
	int id_tx = threadIdx.x;//���ڿ���̺߳ź��ᣬ�ڼ������뵥Ԫ
	int id_ty = threadIdx.y;//���ڿ���̺߳����ᣬ���ڵڼ�������
	int threadPerBlock_x = blockDim.x;//���ڳߴ�ĺ��᳤�ȣ���L
	int threadPerBlock_y = blockDim.y;//���ڳߴ�����᳤�ȣ���N
	int index_real = (id_tx+threadPerBlock_x*id_ty)+id_by*threadPerBlock_x*threadPerBlock_y; //����������ֵ=�ڿ���+����Ŀ���
	if (index_real < size)
	{
		a[index_real].x = a[index_real].x / L;
		a[index_real].y = a[index_real].y / L;
	}
}
//��׼RFT///////////////////////////////////////////////////////
//                        ����ز���              RFT�����     �����ٶȲ������ٶ������ľ��뵥Ԫ���ֵ,�����ظ����,���뵥Ԫ��С������
__global__ void RFT(cufftComplex *pc, cufftComplex *Gv, float Vi, int L, int SP, float Tr , float delt_R, float lamda,cufftComplex *d_DFT, float *d_offset, float V_offset)
{
	cufftComplex Sum={0,0};
	int id_bx = blockIdx.x;//���ڵĿ�ź���,�ڼ�����
	int id_by = blockIdx.y;//���ڵĿ������,�ڼ�����
	int id_tx = threadIdx.x;//���ڿ���̺߳ź���,��ʼ���뵥Ԫ
	int id_ty = threadIdx.y;//���ڿ���̺߳�����,���ڵڼ����ٶ�
	int threadPerBlock_x = blockDim.x;//���ڳߴ�ĺ��᳤,��ʼ����
	int threadPerBlock_y = blockDim.y;//���ڳߴ�����᳤��,һ��������ѵ��ٶȸ���
	int BlockPerGid_x= gridDim.x;//grid�еĿ��������
	int BlockPerGid_y= gridDim.y;//grid�еĿ��������
	//__shared__ cufftComplex Pc_share[M];
	//����������ֵ=      �ڿ���                      +           ����Ŀ���        *           ��Ĵ�С
	int index_real =(id_tx + id_ty*threadPerBlock_x) + (id_bx + id_by*BlockPerGid_x) * (threadPerBlock_x*threadPerBlock_y);
	int V_index = id_ty + threadPerBlock_y * (id_bx + id_by*BlockPerGid_x) ;//�ٶ�����
	float V=Vi*V_index-V_offset;//���߳��������ٶ�
	int Strat_R=id_tx;//��ʼ���뵥Ԫ
	int maxoffset=floor(V*M*Tr/delt_R+0.5);//�������߶�
	float Sum_x=0;
	float Sum_y=0;
	int   offset=0;
	if (Strat_R-maxoffset>=0 && index_real<SP*M*L && Strat_R<L)//û���������������ͼ���
	{
		float fd = 2*V/lamda;//û��DFTϵ��
		for(int i=0;i<M;++i)
		{
			offset=floor(V*i*Tr/delt_R+0.5);
			Strat_R=Strat_R-d_offset[i+id_bx*M+id_by*M*M];
			
			Sum_x=cos(2*PI*fd*i*Tr);//DFTϵ�����߳������
			Sum_y=sin(2*PI*fd*i*Tr);//DFTϵ�����߳������
			Sum.x+=pc[Strat_R+i*L].x * Sum_x - pc[Strat_R+i*L].y * Sum_y;//i;//
			Sum.y+=pc[Strat_R+i*L].x * Sum_y + pc[Strat_R+i*L].y * Sum_x;//i;//
			//�������ȴ洢��DFTϵ����ͨ��ע�����ı�����DFT�洢��ʽ
			//Sum.x+=pc[Strat_R+i*L].x * d_DFT[i+id_bx*M+id_by*M*M].x - pc[Strat_R+i*L].y * d_DFT[i+id_bx*M+id_by*M*M].y;//i;//
			//Sum.y+=pc[Strat_R+i*L].x * d_DFT[i+id_bx*M+id_by*M*M].y + pc[Strat_R+i*L].y * d_DFT[i+id_bx*M+id_by*M*M].x;//i;//
			//Strat_R=id_tx;
		}
		Gv[index_real]=Sum;
		Sum.x= 0;
		Sum.y= 0;
	}
	else 
	{
		Gv[index_real].x=0;
		Gv[index_real].y=0;
	}
	Sum.x= 0;
	Sum.y= 0;
}
//ת��
__global__ void ChangeVector(cufftComplex *a,cufftComplex *b,int L)//����ת��
{
	__shared__ cufftComplex share_block[600];//�潻���õ�
	//cufftComplex temp;//���õ���ʱ���������ڼĴ�������
	int id_tx = threadIdx.x;//���ڿ���̺߳ź���
	int id_bx = blockIdx.x;//���ڵĿ������
	int threadPerBlock_x = blockDim.x;//���ڳߴ�ĺ��᳤��
	if(id_tx < L && id_bx < 2 * M)
	{
		int index_out = id_tx + id_bx*threadPerBlock_x; //���������
		share_block[id_tx] = a[index_out];
	}

	__syncthreads();//���Ϊ������߳�ͬ��

	id_tx = threadIdx.x;//���ڿ���̺߳ź���
	id_bx = blockIdx.x;//���ڵĿ������
	int threadPerBlock_y = blockDim.y;//���ڳߴ�����᳤��
	if(id_tx < L && id_bx < 2 * M)
	{
		int index_in = id_tx*(2*M)+id_bx;//����������
		b[index_in] = share_block[id_tx];
	}

	//temp = a[index_out];
	//a[index_out] = a[index_in];
	//a[index_in] = temp;

}
int main()
{
	//��ʱ��������
	LARGE_INTEGER fp_cpu;//cpu��Ƶ
	QueryPerformanceFrequency(&fp_cpu);//��ȡ��Ƶ
		///GPU��Ϣ˵��
	int MTPB = 1024;//��block���֧���߳���
	int MaxThreadBlockSize[3]={1024,1024,64};//���Ŀ����֧꣬��1024*1024*64����
		/*int dev=1;
	cudaDeviceProp prop;
	cudaGetDevice(&dev);
	printf("GPU�ͺ�:%s\n",prop.name);
	printf("��block���֧���߳���:%d\n",prop.maxThreadsPerBlock);*/
	///���������͸�ֵ
	float fc=100e6;//��Ƶ
	float B=4e6;//����
	float Tao=128e-6;//����
	float Fs=1*B;//����Ƶ��
	float Ts=1/Fs;
	//printf("Ts=%0.8f\n",Ts);
	float mu=B/Tao;//��Ƶ��
	float C=3e8;//����
	float delt_R=C/(2*Fs);//����ֱ���
	float R_start=359*delt_R;//��ʼ����
	float lamda=C/fc;//����
	float PRF=500;//�����ظ�����
	float Tr=1/PRF;//�����ظ�ʱ��
	float Vr=2000;//��ʼ�ٶ�
	printf("Vr=%f\n",Vr);
	float a=0;//���ٶ�
	int    L=Tao*Fs;//���뵥Ԫ���������������
	printf("L=%d\n",L);
	int    N=MTPB/L;//�����߳�����ά����
	printf("N=%d\n",N);
	float *t=(float*)malloc(sizeof(float)*L);//��ʱ��
	float *delt_t=(float*)malloc(sizeof(float)*M);//�ӳ�ʱ��


	for(int i=0;i<L;++i)
	{
		t[i]=-Tao/2+Ts*i;	
		//t[i]=Ts*i;	
	   // printf("t[%d]=%0.8f\n",i+1,t[i]);
	}
	
	///��������/////////////////////////////////////////////
	
	cufftComplex  *h_echo,*d_echo; //������,�豸�λز�ʱ��
	cufftComplex  *h_echo_fft,*d_echo_fft; //������,�豸�λز�Ƶ��

	cufftComplex  *h_ht,*d_ht;    //������,�豸����ѹʱ��ϵ��
	cufftComplex  *h_ht_fft,*d_ht_fft;    //������,�豸����ѹʱ��ϵ��

	cufftComplex  *h_pc,*d_pc;   //������,�豸��ѹʱ����
	cufftComplex  *h_pc_fft,*d_pc_fft;//������,�豸��ѹƵ����

	//�����ڴ棬�Դ�////////////////////////
	cudaError_t cudaStatus;//״̬��¼
    h_echo=(cufftComplex *)malloc(sizeof(cufftComplex )*M*L);//�����˻ز�ʱ��
	h_echo_fft=(cufftComplex *)malloc(sizeof(cufftComplex )*M*L);//�����˻ز�Ƶ��

	h_ht=(cufftComplex *)malloc(sizeof(cufftComplex )*L);//��������ѹʱ��ϵ��
	h_ht_fft=(cufftComplex *)malloc(sizeof(cufftComplex )*L);//��������ѹƵ��ϵ��

	h_pc=(cufftComplex *)malloc(sizeof(cufftComplex )*M*L);//��������ѹʱ����
	h_pc_fft=(cufftComplex *)malloc(sizeof(cufftComplex )*M*L);//��������ѹƵ����

    cudaStatus=cudaMalloc((void**)&d_echo,sizeof(cufftComplex )*M*L);//�豸�λز�ʱ�򿪱��Դ�	    
    if (cudaStatus != cudaSuccess) {
        printf( "d_echo cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_echo_fft,sizeof(cufftComplex )*M*L);//�豸�λز�Ƶ�򿪱��Դ�
	if (cudaStatus != cudaSuccess) {
        printf( "d_echo_fft cudaMalloc failed!\n");
		return 1;
    }

	cudaStatus=cudaMalloc((void**)&d_pc_fft,sizeof(cufftComplex )*M*L);//�豸�λز���ѹƵ��
	if (cudaStatus != cudaSuccess) {
        printf( "d_pc_fft cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_pc,sizeof(cufftComplex )*M*L);//�豸�λز���ѹʱ��
	if (cudaStatus != cudaSuccess) {
        printf( "d_pc cudaMalloc failed!\n");
		return 1;
    }

	cudaStatus=cudaMalloc((void**)&d_ht,sizeof(cufftComplex )*L);//�豸�λز���ѹʱ��ϵ���Դ�
	
	if (cudaStatus != cudaSuccess) {
        printf( "d_ht cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_ht_fft,sizeof(cufftComplex )*L);//�豸�λز���ѹƵ��ϵ���Դ�
	if (cudaStatus != cudaSuccess) {
        printf( "d_ht_fft cudaMalloc failed!\n");
		return 1;
    }
	
	///��ʼ��ֵ///////
	///��������ѹʱ��ϵ��
	for(int i=0;i<=L-1;++i)///
	{
		h_ht[i].x=cos(2*PI*(mu/2)*t[i]*t[i]);
		h_ht[i].y=sin(2*PI*(mu/2)*t[i]*t[i]);
	}
	//���ݺ�����ת����
	/*cufftComplex t_h;//�����õ���ʱ����
	for (int i = 0; i < L / 2; ++i)//���ݺ���
	{
		if (L - i - 1 >= 0)
		{
			t_h.x = h_ht[i].x;
			h_ht[i].x = h_ht[L - i - 1].x;
			h_ht[L - i - 1].x = t_h.x;
			t_h.y = -1 * h_ht[i].y;//����Ӹ���
			h_ht[i].y = -1 * h_ht[L - i - 1].y;
			h_ht[L - i - 1].y = t_h.y;
		}
	}*/
   //for(int i=0;i<=L-1;++i)printf("h_ht[%d]=%0.8f+%0.8fi\n",i+1,h_ht[i].x,h_ht[i].y);
	///�ز���ֵ
	for(int j=0;j<M;++j)//�ӳټ���
	{
			Vr=Vr+a*Tr*j;
		    delt_t[j]=2*(R_start+Vr*Tr*j)/C;
	}
	for(int i=0;i<M;++i)//��ʱ��
	{
		for(int j=0;j<L; ++j)//��ʱ��
		{
			h_echo[j+i*L].x=cos(2*PI*(mu/2*(t[j]+delt_t[i])*(t[j]+delt_t[i])+fc*delt_t[i]));
			h_echo[j+i*L].y=-sin(2*PI*(mu/2*(t[j]+delt_t[i])*(t[j]+delt_t[i])+fc*delt_t[i]));
			//printf("h_echo[%d][%d]=%0.8f+%0.8fi  ",i+1,j+1,h_echo[j+i*L].x,h_echo[j+i*L].y);
			//if((j+1)%3==0) printf("\n");
		}
		//printf("\n");
	}
	//for(int i=L;i<=2*L-1;++i)printf("h_echo[%d]=%0.8f+%0.8fi\n",i+1,h_echo[i].x,h_echo[i].y);
	//���������ݿ������豸��
	cudaStatus=cudaMemcpy(d_echo,h_echo,sizeof(cufftComplex )*M*L,cudaMemcpyHostToDevice);//�ز����������˵��豸��
	if (cudaStatus != cudaSuccess) {
        printf("h_echo->d_echo cudaMemcpy failed!\n");
		return 1;
    }
	cudaStatus=cudaMemcpy(d_ht,h_ht,sizeof(cufftComplex )*L,cudaMemcpyHostToDevice);// ��ѹʱ��ϵ��,���������˵��豸��
    if (cudaStatus != cudaSuccess) {
        printf("h_ht->d_ht cudaMemcpy failed!\n");
		return 1;
    }
	//�Իز�����ѹϵ����fft�任��Ƶ����Ƶ�������ѹ
	//����fft�ƻ�
	LARGE_INTEGER b_pc1,e_pc1,b_pc2,e_pc2;//��ѹ��ʱ
	float time_pc1,time_pc2;
	cufftHandle plan_ML,plan_L;//��L��M����fft����L��1����fft�ƻ�
	cufftResult Result_fft_ML,Result_fft_L;//���ִ�н��������ֵ���ͻ�����cudaError_t�����˹���
	cufftPlan1d(&plan_ML,L,CUFFT_C2C,M);//���ûز�fft�ƻ�
	QueryPerformanceCounter(&b_pc1);
	Result_fft_ML=cufftExecC2C(plan_ML,d_echo,d_echo_fft,CUFFT_FORWARD);//�Իز�����FFT
	//��FFT��0��fsǧ����matlab��fftshift������ݱȶ�
	/*if (Result_fft_ML != cudaSuccess)
	{
		printf("�ز�fft����!\n");	
		return 1;
	}*/
	/*cudaStatus=cudaMemcpy(h_echo_fft,d_echo,sizeof(cufftComplex )*M*L,cudaMemcpyDeviceToHost);//�ز�Ƶ�������豸�ε�������
	if (cudaStatus != cudaSuccess) {
        printf("d_echo_fft->h_echo_fft cudaMemcpy failed!\n");
		return 1;
    }*/
	///
	/*for(int i=0;i<M;++i)//��ʱ��
	{
		for(int j=0;j<L; ++j)//��ʱ��
		{
			printf("h_echo_fft[%d][%d]=%0.8f+%0.8fi  ",i+1,j+1,h_echo_fft[j+i*L].x,h_echo_fft[j+i*L].y);
			if((j+1)%2==0) printf("\n");
		}
		printf("\n");
	}*/
	//for(int i=0;i<=L-1;++i)printf("h_echo_fft[%d]=%0.8f+%0.8fi\n",i+1,h_echo_fft[i].x,h_echo_fft[i].y);
	cufftPlan1d(&plan_L,L,CUFFT_C2C,1);
	Result_fft_L=cufftExecC2C(plan_L,d_ht,d_ht_fft,CUFFT_FORWARD);//����ѹϵ������FFT CUFFT_INVERSE CUFFT_FORWARD
	if (Result_fft_L != cudaSuccess)
	{
		printf("��ѹϵ��fft����!\n");	
		return 1;
	}
	/*cudaStatus=cudaMemcpy(h_ht_fft,d_ht_fft,sizeof(cufftComplex )*L,cudaMemcpyDeviceToHost);//��ѹϵ��Ƶ�������豸�ε�������
	if (cudaStatus != cudaSuccess) {
        printf("d_ht_fft->h_ht_fft cudaMemcpy failed!\n");
		return 1;
    }
	for(int j=0;j<L; ++j)
	{
			printf("h_ht_fft[%d]=%0.8f+%0.8fi \n ",j+1,h_ht_fft[j].x,h_ht_fft[j].y);
	}*/
	//������ѹ/////////////////////////////////////////////
	//GPU���
	dim3 blcok(1,M);//ÿ������һ�����뵥Ԫ
	dim3 threadPerBlock(L,1);//ÿ�����е�һ���߳���һ������
	MulVector <<<blcok,threadPerBlock>>>(d_ht_fft, d_echo_fft, d_pc_fft, L*M);//����ˣ���֪��Ϊʲô���к���
	//������Ƿ���ȷ/////////////////////
	/*cudaStatus=cudaMemcpy(h_pc_fft,d_pc_fft,sizeof(cufftComplex )*M*L,cudaMemcpyDeviceToHost);//��ѹƵ�������豸�ε�������
    if (cudaStatus != cudaSuccess) {
        printf("d_pc_fft->h_pc_fft cudaMemcpy failed!\n");
		return 1;
    }*/
	/*for(int i=0;i<M;++i)
	{
		for(int j=0;j<L; ++j)
		{
			printf("h_pc_fft[%d][%d]=%0.8f+%0.8fi  ",i+1,j+1,h_pc_fft[j+i*L].x,h_pc_fft[j+i*L].y);
			printf("\n");
		}
		printf("\n");
	}*/

	Result_fft_ML=cufftExecC2C(plan_ML,d_pc_fft,d_pc,CUFFT_INVERSE);//��ѹIFFT
	if (Result_fft_ML != cudaSuccess)
	{
		printf("�ز�fft����!\n");	
		return 1;
	}
	ChuVector << <blcok, threadPerBlock >> >(d_pc, M*L,L);//�����
	QueryPerformanceCounter(&e_pc1);
	time_pc1=(float)(e_pc1.QuadPart-b_pc1.QuadPart)/(float)fp_cpu.QuadPart;
	QueryPerformanceCounter(&b_pc2);
	cudaStatus=cudaMemcpy(h_pc,d_pc,sizeof(cufftComplex )*M*L,cudaMemcpyDeviceToHost);//��ѹƵ�������豸�ε�������
	QueryPerformanceCounter(&e_pc2);
	time_pc2=(float)(e_pc2.QuadPart-b_pc2.QuadPart)/(float)fp_cpu.QuadPart;
	printf("��ѹ��ʱ:%f s\n",time_pc1);
	printf("���ݴ�����ʱ:%f s\n",time_pc2);
	printf("����ʱ:%f s\n\n",time_pc1+time_pc2);
    /*if (cudaStatus != cudaSuccess) {
        printf("d_pc_fft->h_pc_fft cudaMemcpy failed!\n");
		return 1;
    }
	for(int i=33;i<34;++i)
	{
		for(int j=0;j<L; ++j)
		{
			printf("h_pc[%d][%d]=%0.8f+%0.8fi  ",i+1,j+1,h_pc[j+i*L].x,h_pc[j+i*L].y);
			printf("\n");
		}
		printf("\n");
	}*/
	/////////��ѹ�������////////////////////////////////////
	float *h_abs_pc=(float*)malloc(sizeof(float)*M*L);
	for(int i=0; i<M*L; ++i)h_abs_pc[i]=sqrt(h_pc[i].x*h_pc[i].x+h_pc[i].y*h_pc[i].y);
	FILE *fp_pc;
	fp_pc=fopen("d:/Pc.txt","w");
	for(int i=0; i<M; ++i)
	{
		for(int j=0; j<L; ++j)
		{
			fprintf(fp_pc,"%0.8f\t",h_abs_pc[j+i*L]);
		}
		fprintf(fp_pc,"\n");
	}
	//////////////////////////////////////////��ѹ���////////////////////////////////////////////////
	///////////////////////////////////////////RFT///////////////////////////////////////////////////////
	//������������
	float Vb = lamda/2/Tr;//ä��
	//printf("Vb=%f\n",Vb);
	float Vi = Vb/M; //�ٶ���������,�ɵ�
	//printf("Vi=%f\n",Vi);
	int    SP =400;//1024*6/M;   //������ä���������//
	float V_offset =SP/2*Vb;
	printf("SP=%d\n",SP);
	printf("�ٶ�����[%f,%f]\n",-SP/2*Vb,SP/2*Vb);

	/////����SP�������ٶȡ���Ӧ�Ķ����մ洢�ã��ڼ����Ӧ��DFTϵ���ٴ洢��//////////////
	float *h_fd, *h_offset, *d_offset;
	h_fd=(float *)malloc(sizeof(float )*M*SP);//������fd�洢
	h_offset=(float *)malloc(sizeof(float )*M*SP*M);//������û���ٶ���,ÿ�������ظ����ڵ�ƫ������
	cudaStatus=cudaMalloc((void**)&d_offset,sizeof(cuComplex )*M*SP*M);//�豸��û���ٶ���,ÿ�������ظ����ڵ�ƫ������
	if (cudaStatus != cudaSuccess) {
        printf("d_offset cudaMalloc failed!\n");
		return 1;
    }
	////����fd,ƫ����//////////
	for(int i=0;i<M*SP;++i)
	{
		h_fd[i] = 2*Vi*i/lamda;
		for(int j=0; j<M; ++j)
		{
			h_offset[i*M+j]=floor(Vi*i*j*Tr/delt_R+0.5);
		}
		
	}
	cudaStatus=cudaMemcpy(d_offset,h_offset,sizeof(float)*M*SP*M,cudaMemcpyHostToDevice);//��������������ƫ�������豸��
		if (cudaStatus != cudaSuccess) {
        printf("h_offset->d_offset  cudaMemcpy failed!\n");
		return 1;
    }
	///////////////����DFTϵ��///////////////
	cuComplex *h_DFT,*d_DFT;
	h_DFT=(cuComplex *)malloc(sizeof(cuComplex)*M*SP*M);//������DFT�洢
	cudaStatus=cudaMalloc((void**)&d_DFT,sizeof(cuComplex )*M*SP*M);//�豸��DFT�洢
	for(int i=0; i<M*SP; ++i)
	{
		for(int j=0; j<M; ++j)
		{
			h_DFT[i*M+j].x=cos(2*PI*h_fd[i]*j*Tr);
			h_DFT[i*M+j].y=sin(2*PI*h_fd[i]*j*Tr);
		}	
	}
	cudaStatus=cudaMemcpy(d_DFT,h_DFT,sizeof(cufftComplex)*M*SP*M,cudaMemcpyHostToDevice);//��������������DFTϵ�����豸��
	if (cudaStatus != cudaSuccess) {
        printf("h_DFT->d_DFT  cudaMemcpy failed!\n");
		return 1;
    }
	free(h_DFT);
	free(h_fd);

	////////////////////////////////////////CPU_RFT//////////////////////////////////////////////////////////
	/*int    Strat_R=0;
	float Sum_x,Sum_y;
	long float offset;
	float2 Sum_cpu;
	cufftComplex *Gv_cpu=(cufftComplex*)malloc(sizeof(cufftComplex)*SP*M*L);//V*L,�ٶ�-���뵥Ԫ��ά
	LARGE_INTEGER b_cpu,e_cpu;//��ʱ
	float time_cpu;
	float fd_cpu;
	float *V=(float*)malloc(sizeof(float)*SP*M);
	int    index_real;
	QueryPerformanceCounter(&b_cpu);
	for(int i=0; i<SP*M; ++i)//�ٶ�
	{
		V[i]=Vi*i-V_offset;
		fd_cpu=2*V[i]/lamda;//������
		//printf("V=%f\t",V[i]);
		//printf("fd_cpu=%f\t",fd_cpu);
		for(int j=0; j<L; ++j)//���뵥Ԫ
		{
			Sum_cpu.x=0;
			Sum_cpu.y=0;
			Strat_R=j;
			//printf("Strat_R=%d\n",Strat_R);
			if(Strat_R -floor(V[i]*Tr*M/delt_R+0.5)>=10 && Strat_R + V[i]*Tr*M/delt_R<L-10)
			{
				for(int ti=0; ti<M; ++ti)
				{
					offset=floor(V[i]*ti*Tr/delt_R+0.5);
					Strat_R=Strat_R-floor(offset);
					Sum_x=cos(2*PI*fd_cpu*ti*Tr);
					Sum_y=sin(2*PI*fd_cpu*ti*Tr);
					index_real=Strat_R+ti*L;
					Strat_R=j;
					//printf("offset=%d\t index_real=%d\t\n",offset,index_real+1);
					Sum_cpu.x+=h_pc[index_real].x * Sum_x - h_pc[index_real].y * Sum_y;//i;//
					Sum_cpu.y+=h_pc[index_real].x * Sum_y + h_pc[index_real].y * Sum_x;//i;//
					//printf("h_pc=%0.8f+%0.8fi\n",h_pc[index_real].x,h_pc[index_real].y);
					//printf("Sum=%0.8f+%0.8fi\n",Sum_x,Sum_y);
					//Gv_cpu[j+i*L].x+=h_pc[Strat_R+ti*L].x * Sum_x - h_pc[Strat_R+ti*L].y *Sum_y;//i;//
					//Gv_cpu[j+i*L].y+=h_pc[Strat_R+ti*L].x * Sum_y + h_pc[Strat_R+ti*L].y * Sum_x;//i;//������ô�ۼӣ������������ͷ���
				}	
			}
			Gv_cpu[j+i*L]=Sum_cpu;
		}	
	}

	QueryPerformanceCounter(&e_cpu);
	time_cpu=(float(e_cpu.QuadPart-b_cpu.QuadPart))/(float)(fp_cpu.QuadPart);
	printf("CPU-RFT����ʱ��:%f s\n\n",time_cpu);
	fd_cpu=2*Vi*(SP*M-1)/lamda;*/
	//for(int i=0; i<M;++i)printf("exp=%f+%fi\n",cos(2*PI*fd_cpu*i*Tr),sin(2*PI*fd_cpu*i*Tr));////fd
	//for(int i=0; i<SP*M;++i)printf("V=%0.8f\n",V[i]);

	//for(int i=0; i<SP*M;++i)printf("Maxoffset=%0.8f\n",0 -int( V[i]*Tr*M/delt_R));
	/*for(int i=0; i<1; ++i)//SP*M
	{
		for(int j=0; j<L; ++j)printf("Gv_cpu[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,Gv_cpu[j+i*L].x,Gv_cpu[j+i*L].y);
	}*/
	////////������RFT�ı����/////////////////////////////////////////
	/*float *abs_Gv_cpu=(float*)malloc(sizeof(float)*SP*M*L);
	for(int i=0; i<SP*M*L; ++i)//SP*M
	{
		abs_Gv_cpu[i]=sqrt((Gv_cpu[i].x)*(Gv_cpu[i].x)+(Gv_cpu[i].y)*(Gv_cpu[i].y));
		if(abs_Gv_cpu[i]>1e6)abs_Gv_cpu[i]=0;
	}
	FILE *fp_Gv_cpu;
	fp_Gv_cpu=fopen("d:/GV_cpu.txt","w");
	for(int i=0; i<SP*M; ++i)//SP*M
	{
		for(int j=0; j<L; ++j)//SP*M
		{
			fprintf(fp_Gv_cpu,"%f\t",abs_Gv_cpu[j+i*L]);
		}
		fprintf(fp_Gv_cpu,"\n");
	}*/
	////////////////////////////////�����ֻ������ѹ���ʱ���Ƶ����
    cudaFree(d_echo);
    cudaFree(d_echo_fft);
	cudaFree(d_ht);
	cudaFree(d_ht_fft);
	/////////////////////////////////////////////////GPU_RFT��ʼ//////////////////////////////////////////////
	//RFT��������Ϳ��ٿռ�
	cufftComplex *h_Gv, *d_Gv;//�����ˣ��豸�α�׼RFT�������
	h_Gv=(cufftComplex*)malloc(sizeof(cufftComplex)*SP*M*L);
	cudaStatus=cudaMalloc((void**)&d_Gv,sizeof(cufftComplex)*SP*M*L);
	if (cudaStatus != cudaSuccess)
	{
		printf("d_Gv cudaMalloc fail!\n");
		return 1;
	}
	float DataQ=(float)sizeof(cufftComplex)*SP*M*L/1024/1024;//RFT���������
	printf("������%fMB\n\n",DataQ);
	dim3 block_s(M,SP);//M,SP
	dim3 threadPerBlock_s(L,1);
	
	LARGE_INTEGER b1,b2,e1,e2;//��ʼʱ�䣬��������ʱ�䣬��������ʱ��
	float time1,time2;
	cudaEvent_t RFT_start,RFT_end;
	float time_RFT;
	cudaEventCreate(&RFT_start);
	cudaEventCreate(&RFT_end);
	cudaEventRecord(RFT_start,0);
	//QueryPerformanceCounter(&b1);
	RFT<<<block_s,threadPerBlock_s>>>(d_pc,d_Gv,Vi,L,SP,Tr,delt_R,lamda,d_DFT,d_offset,V_offset);
	cudaEventRecord(RFT_end,0);
	cudaEventSynchronize(RFT_end);
	cudaEventElapsedTime(&time_RFT,RFT_start,RFT_end);
	printf("GPU-RFT����ʱ��event��ʱ:%0.8f s\n",time_RFT/1000);
	cudaEventDestroy(RFT_start);
	cudaEventDestroy(RFT_end);
	//QueryPerformanceCounter(&e1);
	//time1=(float)(e1.QuadPart-b1.QuadPart)/(float)fp_cpu.QuadPart;
	//printf("GPU-RFT����ʱ��:%0.8f s\n",time1);

	//Ҫ���ͷŵ���������� 
	cudaFree(d_pc);//
	cudaFree(d_pc_fft);
	////
	//QueryPerformanceCounter(&b2);
	cudaEvent_t RFT_t_start,RFT_t_end;//����
	float time_RFT_trans;//����ʱ��
	cudaEventCreate(&RFT_t_start);
	cudaEventCreate(&RFT_t_end);
	cudaEventRecord(RFT_t_start,0);
	cudaStatus=cudaMemcpy(h_Gv,d_Gv,sizeof(cufftComplex)*SP*M*L,cudaMemcpyDeviceToHost);// RFT���,���������˵��豸��
	cudaEventRecord(RFT_t_end,0);
	cudaEventSynchronize(RFT_t_end);
	cudaEventElapsedTime(&time_RFT_trans,RFT_t_start,RFT_t_end);
	cudaEventDestroy(RFT_t_start);
	cudaEventDestroy(RFT_t_end);
	//QueryPerformanceCounter(&e2);
    if (cudaStatus != cudaSuccess) {
        printf("d_Gv->h_Gv cudaMemcpy failed!\n Error Code:%d",cudaStatus);
		return 1;
    }
	printf("RFT���ݴ���ʱ��event��ʱ:%0.8f s\n",time_RFT_trans/1000);
	printf("RFT����ʱevent��ʱ:%0.8f s\n",(time_RFT_trans+time_RFT)/1000);
	//time2=(float)(e2.QuadPart-b2.QuadPart)/(float)fp_cpu.QuadPart;
	//printf("���ݴ���ʱ��:%0.8f s\n",time2);
	//printf("����ʱ:%0.8f s\n",(time2+time1));
	float TransSpeed =DataQ/time_RFT_trans;
	//int TransSpeed =DataQ/time2;
	printf("�����ٶ�Ϊ%.2fMB s\n\n",TransSpeed*1000);

//	float Speedup=time_cpu/(time_RFT_trans+time_RFT)*1000;
	//printf("���ٱ�:%f\n",Speedup);
	/*for(int i=10*M-2; i<10*M-1; ++i)//SP*M
	{
		for(int j=0; j<L; ++j)printf("h_Gv[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,h_Gv[j+i*L].x,h_Gv[j+i*L].y);
	}*/
	/////����Դ�///////////////////////////////////////////////////////////
	//cudaFree(d_Gv);
	/////////////////////////////////////////////////////////////////////////////�ı����/////
	float *h_abs_Gv=(float*)malloc(sizeof(float)*SP*M*L);
	for(int i=0; i<SP*M*L; ++i)//SP*M
	{
		h_abs_Gv[i]=sqrt((h_Gv[i].x)*(h_Gv[i].x)+(h_Gv[i].y)*(h_Gv[i].y));
		if(h_abs_Gv[i]<0) h_abs_Gv[i]=0;//h_abs_Gv[i]>M*Fs*Tao||
	}
	FILE *fp_Gv;
	fp_Gv=fopen("d:/GV.txt","w");
	for(int i=0; i<SP*M; ++i)//SP*M
	{
		for(int j=0; j<L; ++j)//SP*M
		{
			fprintf(fp_Gv,"%.2f\t",h_abs_Gv[j+i*L]);
		}
		fprintf(fp_Gv,"\n");
	}
	//free(Gv_cpu);
	//free(h_abs_Gv);
	//free(h_Gv);
    return 0;
}

