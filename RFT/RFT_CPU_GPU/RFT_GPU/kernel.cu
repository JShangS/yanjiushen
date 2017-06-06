
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
#define M 512//���������
///kernel ����
//���
__global__ void MulVector(cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c,int size)
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
//x.*ww
__global__ void MulVector_xw(cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c,int size)
{
	int id_by = blockIdx.y;//���ڵĿ������,�ڼ�����,�ڼ�����ʱ����뵥Ԫ
	int id_bx = blockIdx.x;//��ĺ����ʾ�ڼ���ģ������
	int id_tx = threadIdx.x;//���ڿ���̺߳ź��ᣬ�ڼ������뵥Ԫ
	int threadPerBlock_x = blockDim.x;//���ڳߴ�ĺ��᳤�ȣ���2*M
	int threadPerBlock_y = blockDim.y;//���ڳߴ�����᳤�ȣ���1
	int index_x=id_tx+id_by*threadPerBlock_x;//xҪ��˵�����ֵ
	int index_y=index_x;//y���Ҫ���������ֵ
	int index_w=id_tx+size-1+id_by*threadPerBlock_x;//wҪ��˵�����ֵ
	if (id_tx < size)
	{
		c[index_y].x = a[index_x].x * b[index_w].x - a[index_x].y * b[index_w].y;//blockDim.x;//index_ht;//
		c[index_y].y = a[index_x].x * b[index_w].y + a[index_x].y * b[index_w].x;//blockDim.y;//index_real;//
	}
}
//���
__global__ void ChuVector(cufftDoubleComplex *a, int size,int L)
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
__global__ void RFT(cufftDoubleComplex *pc, cufftDoubleComplex *Gv, double Vi, int L, int SP, double Tr , double delt_R, double lamda,cufftDoubleComplex *d_DFT, double *d_offset)
{
	cufftDoubleComplex Sum={0,0};
	int id_bx = blockIdx.x;//���ڵĿ�ź���,�ڼ�����
	int id_by = blockIdx.y;//���ڵĿ������,�ڼ�����
	int id_tx = threadIdx.x;//���ڿ���̺߳ź���,��ʼ���뵥Ԫ
	int id_ty = threadIdx.y;//���ڿ���̺߳�����,���ڵڼ����ٶ�
	int threadPerBlock_x = blockDim.x;//���ڳߴ�ĺ��᳤,��ʼ����
	int threadPerBlock_y = blockDim.y;//���ڳߴ�����᳤��,һ��������ѵ��ٶȸ���
	int BlockPerGid_x= gridDim.x;//grid�еĿ��������
	int BlockPerGid_y= gridDim.y;//grid�еĿ��������
	//__shared__ cufftDoubleComplex Pc_share[M];
	//����������ֵ=      �ڿ���                      +           ����Ŀ���        *           ��Ĵ�С
	int index_real =(id_tx + id_ty*threadPerBlock_x) + (id_bx + id_by*BlockPerGid_x) * (threadPerBlock_x*threadPerBlock_y);
	int V_index = id_ty + threadPerBlock_y * (id_bx + id_by*BlockPerGid_x) ;//�ٶ�����
	double V=Vi*V_index;//���߳��������ٶ�
	int Strat_R=id_tx;//��ʼ���뵥Ԫ
	int maxoffset=floor(V*M*Tr/delt_R+0.5);//�������߶�
	double Sum_x=0;
	double Sum_y=0;
	int   offset=0;
	if (Strat_R-maxoffset>=0 && index_real<SP*M*L && Strat_R<L)//û���������������ͼ���
	{
		//double fd = 2*V/lamda;
		for(int i=0;i<M;++i)
		{
			//offset=floor(V*i*Tr/delt_R+0.5);
			Strat_R=Strat_R-d_offset[i+id_bx*M+id_by*M*M];
			/*
			Sum_x=cos(2*PI*fd*i*Tr);
			Sum_y=sin(2*PI*fd*i*Tr);
			Sum.x+=pc[Strat_R+i*L].x * Sum_x - pc[Strat_R+i*L].y * Sum_y;//i;//
			Sum.y+=pc[Strat_R+i*L].x * Sum_y + pc[Strat_R+i*L].y * Sum_x;//i;//
			*/
			Sum.x+=pc[Strat_R+i*L].x * d_DFT[i+id_bx*M+id_by*M*M].x - pc[Strat_R+i*L].y * d_DFT[i+id_bx*M+id_by*M*M].y;//i;//
			Sum.y+=pc[Strat_R+i*L].x * d_DFT[i+id_bx*M+id_by*M*M].y + pc[Strat_R+i*L].y * d_DFT[i+id_bx*M+id_by*M*M].x;//i;//
			Strat_R=id_tx;
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
__global__ void ChangeVector(cufftDoubleComplex *a,cufftDoubleComplex *b,int L)//����ת��
{
	__shared__ cufftDoubleComplex share_block[600];//�潻���õ�
	//cufftDoubleComplex temp;//���õ���ʱ���������ڼĴ�������
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
	double fc=100e6;//��Ƶ
	double B=4e6;//����
	double Tao=128e-6;//����
	double Fs=2*B;//����Ƶ��
	double Ts=1/Fs;
	//printf("Ts=%0.8f\n",Ts);
	double mu=B/Tao;//��Ƶ��
	double C=3e8;//����
	double delt_R=C/(2*Fs);//����ֱ���
	double R_start=500*delt_R;//��ʼ����
	double lamda=C/fc;//����
	double PRF=500;//�����ظ�����
	double Tr=1/PRF;//�����ظ�ʱ��
	double Vr=1200;//��ʼ�ٶ�
	printf("Vr=%f\n",Vr);
	double a=0;//���ٶ�
	int    L=Tao*Fs;//���뵥Ԫ���������������
	printf("L=%d\n",L);
	int    N=MTPB/L;//�����߳�����ά����
	printf("N=%d\n",N);
	double *t=(double*)malloc(sizeof(double)*L);//��ʱ��
	double *delt_t=(double*)malloc(sizeof(double)*M);//�ӳ�ʱ��


	for(int i=0;i<L;++i)
	{
		t[i]=-Tao/2+Ts*i;	
		//t[i]=Ts*i;	
	   // printf("t[%d]=%0.8f\n",i+1,t[i]);
	}
	
	///��������/////////////////////////////////////////////
	
	cufftDoubleComplex  *h_echo,*d_echo; //������,�豸�λز�ʱ��
	cufftDoubleComplex  *h_echo_fft,*d_echo_fft; //������,�豸�λز�Ƶ��

	cufftDoubleComplex  *h_ht,*d_ht;    //������,�豸����ѹʱ��ϵ��
	cufftDoubleComplex  *h_ht_fft,*d_ht_fft;    //������,�豸����ѹʱ��ϵ��

	cufftDoubleComplex  *h_pc,*d_pc;   //������,�豸��ѹʱ����
	cufftDoubleComplex  *h_pc_fft,*d_pc_fft;//������,�豸��ѹƵ����

	//�����ڴ棬�Դ�////////////////////////
	cudaError_t cudaStatus;//״̬��¼
    h_echo=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*M*L);//�����˻ز�ʱ��
	h_echo_fft=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*M*L);//�����˻ز�Ƶ��

	h_ht=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*L);//��������ѹʱ��ϵ��
	h_ht_fft=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*L);//��������ѹƵ��ϵ��

	h_pc=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*M*L);//��������ѹʱ����
	h_pc_fft=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*M*L);//��������ѹƵ����

    cudaStatus=cudaMalloc((void**)&d_echo,sizeof(cufftDoubleComplex )*M*L);//�豸�λز�ʱ�򿪱��Դ�	    
    if (cudaStatus != cudaSuccess) {
        printf( "d_echo cudaMalloc failed!\n %d \n",cudaStatus);
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_echo_fft,sizeof(cufftDoubleComplex )*M*L);//�豸�λز�Ƶ�򿪱��Դ�
	if (cudaStatus != cudaSuccess) {
        printf( "d_echo_fft cudaMalloc failed!\n");
		return 1;
    }

	cudaStatus=cudaMalloc((void**)&d_pc_fft,sizeof(cufftDoubleComplex )*M*L);//�豸�λز���ѹƵ��
	if (cudaStatus != cudaSuccess) {
        printf( "d_pc_fft cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_pc,sizeof(cufftDoubleComplex )*M*L);//�豸�λز���ѹʱ��
	if (cudaStatus != cudaSuccess) {
        printf( "d_pc cudaMalloc failed!\n");
		return 1;
    }

	cudaStatus=cudaMalloc((void**)&d_ht,sizeof(cufftDoubleComplex )*L);//�豸�λز���ѹʱ��ϵ���Դ�
	
	if (cudaStatus != cudaSuccess) {
        printf( "d_ht cudaMalloc failed! \n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_ht_fft,sizeof(cufftDoubleComplex )*L);//�豸�λز���ѹƵ��ϵ���Դ�
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
	/*cufftDoubleComplex t_h;//�����õ���ʱ����
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
	cudaStatus=cudaMemcpy(d_echo,h_echo,sizeof(cufftDoubleComplex )*M*L,cudaMemcpyHostToDevice);//�ز����������˵��豸��
	if (cudaStatus != cudaSuccess) {
        printf("h_echo->d_echo cudaMemcpy failed!\n");
		return 1;
    }
	cudaStatus=cudaMemcpy(d_ht,h_ht,sizeof(cufftDoubleComplex )*L,cudaMemcpyHostToDevice);// ��ѹʱ��ϵ��,���������˵��豸��
    if (cudaStatus != cudaSuccess) {
        printf("h_ht->d_ht cudaMemcpy failed!\n");
		return 1;
    }
	//�Իز�����ѹϵ����fft�任��Ƶ����Ƶ�������ѹ
	//����fft�ƻ�
	LARGE_INTEGER b_pc1,e_pc1,b_pc2,e_pc2;//��ѹ��ʱ
	double time_pc1,time_pc2;
	
	cufftHandle plan_ML,plan_L;//��L��M����fft����L��1����fft�ƻ�
	cufftResult Result_fft_ML,Result_fft_L;//���ִ�н��������ֵ���ͻ�����cudaError_t�����˹���
	cufftPlan1d(&plan_ML,L,CUFFT_Z2Z,M);//���ûز�fft�ƻ�
	QueryPerformanceCounter(&b_pc1);
	Result_fft_ML=cufftExecZ2Z(plan_ML,d_echo,d_echo_fft,CUFFT_FORWARD);//�Իز�����FFT
	//��FFT��0��fsǧ����matlab��fftshift������ݱȶ�
	/*if (Result_fft_ML != cudaSuccess)
	{
		printf("�ز�fft����!\n");	
		return 1;
	}*/
	/*cudaStatus=cudaMemcpy(h_echo_fft,d_echo,sizeof(cufftDoubleComplex )*M*L,cudaMemcpyDeviceToHost);//�ز�Ƶ�������豸�ε�������
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
	cufftPlan1d(&plan_L,L,CUFFT_Z2Z,1);
	Result_fft_L=cufftExecZ2Z(plan_L,d_ht,d_ht_fft,CUFFT_FORWARD);//����ѹϵ������FFT CUFFT_INVERSE CUFFT_FORWARD
	if (Result_fft_L != cudaSuccess)
	{
		printf("��ѹϵ��fft����!\n");	
		return 1;
	}
	/*cudaStatus=cudaMemcpy(h_ht_fft,d_ht_fft,sizeof(cufftDoubleComplex )*L,cudaMemcpyDeviceToHost);//��ѹϵ��Ƶ�������豸�ε�������
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
	/*cudaStatus=cudaMemcpy(h_pc_fft,d_pc_fft,sizeof(cufftDoubleComplex )*M*L,cudaMemcpyDeviceToHost);//��ѹƵ�������豸�ε�������
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

	Result_fft_ML=cufftExecZ2Z(plan_ML,d_pc_fft,d_pc,CUFFT_INVERSE);//��ѹIFFT
	if (Result_fft_ML != cudaSuccess)
	{
		printf("�ز�fft����!\n");	
		return 1;
	}
	ChuVector << <blcok, threadPerBlock >> >(d_pc, M*L,L);//�����
	QueryPerformanceCounter(&e_pc1);
	time_pc1=(double)(e_pc1.QuadPart-b_pc1.QuadPart)/(double)fp_cpu.QuadPart;
	QueryPerformanceCounter(&b_pc2);
	cudaStatus=cudaMemcpy(h_pc,d_pc,sizeof(cufftDoubleComplex )*M*L,cudaMemcpyDeviceToHost);//��ѹƵ�������豸�ε�������
	QueryPerformanceCounter(&e_pc2);
	time_pc2=(double)(e_pc2.QuadPart-b_pc2.QuadPart)/(double)fp_cpu.QuadPart;
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
	double *h_abs_pc=(double*)malloc(sizeof(double)*M*L);
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
	double Vb = lamda/2/Tr;//ä��
	//printf("Vb=%f\n",Vb);
	double Vi = Vb/M; //�ٶ���������
	//printf("Vi=%f\n",Vi);
	int    SP =8;//1024*6/M;   //������ä���������//��������һ��ͺ�����GTX650�Կ������Կ������Ƴ���
	printf("SP=%d\n",SP);


	/////����SP�������ٶȡ���Ӧ�Ķ����մ洢�ã��ڼ����Ӧ��DFTϵ���ٴ洢��//////////////
	double *h_fd, *h_offset, *d_offset;
	h_fd=(double *)malloc(sizeof(double )*M*SP);//������fd�洢
	h_offset=(double *)malloc(sizeof(double )*M*SP*M);//������û���ٶ���,ÿ�������ظ����ڵ�ƫ������
	cudaStatus=cudaMalloc((void**)&d_offset,sizeof(cuDoubleComplex )*M*SP*M);//�豸��û���ٶ���,ÿ�������ظ����ڵ�ƫ������
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
	cudaStatus=cudaMemcpy(d_offset,h_offset,sizeof(double)*M*SP*M,cudaMemcpyHostToDevice);//��������������ƫ�������豸��
		if (cudaStatus != cudaSuccess) {
        printf("h_offset->d_offset  cudaMemcpy failed!\n");
		return 1;
    }
	///////////////����DFTϵ��///////////////
	cuDoubleComplex *h_DFT,*d_DFT;
	h_DFT=(cuDoubleComplex *)malloc(sizeof(cuDoubleComplex)*M*SP*M);//������DFT�洢
	cudaStatus=cudaMalloc((void**)&d_DFT,sizeof(cuDoubleComplex )*M*SP*M);//�豸��DFT�洢
	for(int i=0; i<M*SP; ++i)
	{
		for(int j=0; j<M; ++j)
		{
			h_DFT[i*M+j].x=cos(2*PI*h_fd[i]*j*Tr);
			h_DFT[i*M+j].y=sin(2*PI*h_fd[i]*j*Tr);
		}	
	}
	cudaStatus=cudaMemcpy(d_DFT,h_DFT,sizeof(cufftDoubleComplex)*M*SP*M,cudaMemcpyHostToDevice);//��������������DFTϵ�����豸��
	if (cudaStatus != cudaSuccess) {
        printf("h_DFT->d_DFT  cudaMemcpy failed!\n");
		return 1;
    }
	free(h_DFT);
	free(h_fd);

	////////////////////////////////////////CPU_RFT//////////////////////////////////////////////////////////
	/*int    Strat_R=0;
	double Sum_x,Sum_y;
	long double offset;
	double2 Sum_cpu;
	cufftDoubleComplex *Gv_cpu=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*SP*M*L);//V*L,�ٶ�-���뵥Ԫ��ά
	LARGE_INTEGER b_cpu,e_cpu;//��ʱ
	double time_cpu;
	double fd_cpu;
	double *V=(double*)malloc(sizeof(double)*SP*M);
	int    index_real;
	QueryPerformanceCounter(&b_cpu);
	for(int i=0; i<SP*M; ++i)//�ٶ�
	{
		V[i]=Vi*i;
		fd_cpu=2*V[i]/lamda;//������
		//printf("V=%f\t",V[i]);
		//printf("fd_cpu=%f\t",fd_cpu);
		for(int j=0; j<L; ++j)//���뵥Ԫ
		{
			Sum_cpu.x=0;
			Sum_cpu.y=0;
			Strat_R=j;
			//printf("Strat_R=%d\n",Strat_R);
			if(Strat_R -floor(V[i]*Tr*M/delt_R+0.5)>=0 && Strat_R + V[i]*Tr*M/delt_R<L)
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
	time_cpu=(double(e_cpu.QuadPart-b_cpu.QuadPart))/(double)(fp_cpu.QuadPart);
	printf("CPU-RFT����ʱ��:%f s\n\n",time_cpu);
	fd_cpu=2*Vi*(SP*M-1)/lamda;*/
	//for(int i=0; i<M;++i)printf("exp=%f+%fi\n",cos(2*PI*fd_cpu*i*Tr),sin(2*PI*fd_cpu*i*Tr));////fd
	//for(int i=0; i<SP*M;++i)printf("V=%0.8f\n",V[i]);

	//for(int i=0; i<SP*M;++i)printf("Maxoffset=%0.8f\n",0 -int( V[i]*Tr*M/delt_R));
	/*for(int i=0; i<1; ++i)//SP*M
	{
		for(int j=0; j<L; ++j)printf("Gv_cpu[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,Gv_cpu[j+i*L].x,Gv_cpu[j+i*L].y);
	}
	////////������RFT�ı����/////////////////////////////////////////
	double *abs_Gv_cpu=(double*)malloc(sizeof(double)*SP*M*L);
	for(int i=0; i<SP*M*L; ++i)//SP*M
	{
		abs_Gv_cpu[i]=sqrt((Gv_cpu[i].x)*(Gv_cpu[i].x)+(Gv_cpu[i].y)*(Gv_cpu[i].y));
		//if(abs_Gv_cpu[i]>1e5)abs_Gv_cpu[i]=0;
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
	/////////////////////////////////////////////CZT_RFT_GPU��ʼ/////////////////////////////////////////////////////
	//ϵ�������͸�ֵ////
	/*//////�μ�matlab��czt��ʵ��
	double fai=lamda*B/L/C;//��ʱ��Ƶ��ֱ���
	double *fai_a=(double*)malloc(sizeof(double)*L);//1-fai�Ľ��
	//cufftDoubleComplex *w=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*L);//��ͬ��ʱ��Ƶ�ʵ㴦�Ļ���cztϵ��
	for(int i=0;i<L;++i)//ͬһ����ʱ�䴦��һ��fai_a
	{
		fai_a[i]=1.0f-fai*(i+1.0f);
		//w=exp(-1j*2*pi*fai_a/L);
		//w[i].x=cos(2*PI*fai_a[i]/M);
		//w[i].y=-sin(2*PI*fai_a[i]/M);
		//printf("fai_a[%d]=%0.8f\n",i+1,fai_a[i]);
	}
	double *kk=(double*)malloc(sizeof(double)*2*M-1);//�ݴ�
	double *kk2=(double*)malloc(sizeof(double)*2*M-1);//�ݴ�ƽ��
	for(int i=0;i<2*M;++i)
	{
		kk[i]=-1.0f*M+i+1;
		kk2[i]=kk[i]*kk[i]/2.0f;
		//printf("kk2[%d]=%0.8f\n",i+1,kk2[i]);
	}
	cufftDoubleComplex *h_ww=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*L*2*M);//��ͬ��ʱ��Ƶ�ʵ㴦��cztϵ��
	cufftDoubleComplex *h_v=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*L*2*M);//��ͬ��ʱ��Ƶ�ʵ㴦��czt 1./ϵ��
	///ww��L�У�2*M�У�����ÿһ�ж�Ӧһ����ʱ��Ƶ�㣬��ͬһ����ʱ��Ƶ�㴦��M������ز���CZT�任
	for(int i=0;i<L;++i)
	{
		//w=exp(-1j*2*pi*fai_a/L);
		//ww=w.^kk2
		//v=1./ww;
		for(int j=0; j<2*M; ++j)///���������matlba��ת�õĹ�ϵ 
		{
			//int indx_a=int(i/2/M);//ͬһ����ʱ��Ƶ��
			//int indx_kk2=i%(2*M);//��2*Mѭ��
			//printf("indx_a=%d,indx_kk2=%d\n",indx_a,indx_kk2);
			h_ww[i*2*M+j].x=cos(2*PI*fai_a[i]*kk2[j]/M);
			h_ww[i*2*M+j].y=-sin(2*PI*fai_a[i]*kk2[j]/M);
			h_v[i*2*M+j].x=cos(2*PI*fai_a[i]*kk2[j]/M);
			h_v[i*2*M+j].y=sin(2*PI*fai_a[i]*kk2[j]/M);
		}

		
	}
	
	cufftDoubleComplex *d_ww;
	
	cudaStatus=cudaMalloc((void**)&d_ww,sizeof(cufftDoubleComplex)*2*M*L);//��GPU�˿���һ�齫�źŽ�����ά��ԭ����M*Lת�ü���ά��L*M*2������ѹƵ��ֵ��d_x
	if (cudaStatus != cudaSuccess) {
        printf("d_ww cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMemcpy(d_ww,h_ww,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyHostToDevice);//�������ݵ�x����czt
	if (cudaStatus != cudaSuccess) {
        printf("h_ww->d_ww cudaMemcpy failed!\n");
		return 1;
    }
	/////��һ��λ�ã��Ѻ���M��wֵ��ǰ��M��wֵλ�û���//////////////
	cufftDoubleComplex *d_ww_change;
	cudaStatus=cudaMalloc((void**)&d_ww_change,sizeof(cufftDoubleComplex)*2*M*L);//��GPU�˿���һ�齫�źŽ�����ά��ԭ����M*Lת�ü���ά��L*M*2������ѹƵ��ֵ��d_x
	if (cudaStatus != cudaSuccess) {
        printf("d_ww_change cudaMalloc failed!\n");
		return 1;
    }
	for(int i=0; i<L;++i)
	{
		//cudaMemcpy(&d_ww_change[i*2*M],&d_ww[i*2*M+M],sizeof(cufftDoubleComplex)*M,cudaMemcpyDeviceToDevice);
		//cudaMemcpy(&d_ww_change[i*2*M+M],&d_ww[i*2*M],sizeof(cufftDoubleComplex)*M,cudaMemcpyDeviceToDevice);
	}
	cudaMemcpy(d_ww,d_ww_change,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyDeviceToDevice);
	cudaStatus=cudaMemcpy(h_ww,d_ww,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyDeviceToHost);*/
	//////////////////////////////////////////////////////////////////
	//���ww�Ƿ���ȷ,û��2016/10/11///////////////////
	/*for(int i=0;i<1;++i)//��ʱ�䣬��
	{
		for(int j=255;j<275;++j)//��ʱ�䣬��
		{
			printf("v[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,h_v[i*2*M+j].x,h_v[i*2*M+j].y);
		}
		
	}*/

	/*cufftDoubleComplex *d_v;
	cudaStatus=cudaMalloc((void**)&d_v,sizeof(cufftDoubleComplex)*2*M*L);//��GPU�˿���һ�齫�źŽ�����ά��ԭ����M*Lת�ü���ά��L*M*2������ѹƵ��ֵ��d_x
	if (cudaStatus != cudaSuccess) {
        printf("d_v cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMemcpy(d_v,h_v,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyHostToDevice);//�������ݵ�x����czt
	if (cudaStatus != cudaSuccess) {
        printf("h_v->d_v cudaMemcpy failed!\n");
		return 1;
    }
	//��չx��FFTǰ��׼��
	cufftDoubleComplex *zeros_x=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*1*M*L);//һ��M��L�е�0���󣬸�d_x���θ�ֵ
	for(int i = 0; i < 1 *M * L; ++i)
	{
		zeros_x[i].x = 0;
		zeros_x[i].y = 0;
		//if(i%100==0) printf("%0.8f+%0.8fi\n",zeros_x[i].x,zeros_x[i].y);

	}
	cufftDoubleComplex *d_x;//pc_fft����չ����
	cudaStatus=cudaMalloc((void**)&d_x,sizeof(cufftDoubleComplex)*2*M*L);//��GPU�˿���һ�齫�źŽ�����ά��ԭ����M*Lת�ü���ά��L*M*2������ѹƵ��ֵ��d_x
	if (cudaStatus != cudaSuccess) {
        printf("d_x cudaMalloc failed!\n");
		return 1;
    }
	//��Ƶ����ѹ�������
	//��0����
	cudaStatus=cudaMemcpy(&d_x[0],d_pc_fft,sizeof(cufftDoubleComplex)*M*L,cudaMemcpyDeviceToDevice);//�������ݵ�x����czt
	if (cudaStatus != cudaSuccess) {
        printf("d_pc_fft->d_x cudaMemcpy failed!\n");
		return 1;
    }//��0����x����
	cudaStatus=cudaMemcpy(&d_x[M*L],zeros_x,sizeof(cufftDoubleComplex)*1*M*L,cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
        printf("zeros_x->d_x cudaMemcpy failed!\n");
		return 1;
    }*/
    
	
	///////////ת��//////////////////////////////////////��ʼCZT//////////////////////////
	/*cufftDoubleComplex *d_x_change;//ת�ú��d_x
	cudaMalloc((void**)&d_x_change,sizeof(cufftDoubleComplex)*L*2*M);//
	dim3 blcok_chang(2*M, 1);
	dim3 threadPerBlock_chang(L, 1);
	////////////��ʱ��������////////////////////
	LARGE_INTEGER b_czt_change,e_czt_change,b_czt,e_czt,b_diancheng,e_diancheng,b_fft,e_fft;//��ʱ,chang:ת��
	LARGE_INTEGER b_dianchu,e_dianchu, b_ifft, e_ifft;
	double time_czt_change,time_czt,time_diancheng,time_fft,time_ifft,time_dianchu;
	cudaEvent_t stop_e_change,start_e_change;
	float time_e_change;
	cudaEventCreate(&start_e_change);
	cudaEventCreate(&stop_e_change);
	//QueryPerformanceCounter(&b_czt_change);
	cudaEventRecord(start_e_change,0);
	ChangeVector << <blcok_chang, threadPerBlock_chang >> >(d_x, d_x_change,L);//����ת��
	cudaEventRecord(stop_e_change,0);
	cudaEventSynchronize(stop_e_change);
	cudaEventElapsedTime(&time_e_change,start_e_change,stop_e_change);

	cudaEventDestroy(stop_e_change);
	cudaEventDestroy(start_e_change);
	//QueryPerformanceCounter(&e_czt_change);
	//time_czt_change =(double(e_czt_change.QuadPart-b_czt_change.QuadPart))/(double)(fp_cpu.QuadPart);*/

	///���ת���Ƿ�ɹ�////////////////////////
	/*cufftDoubleComplex *h_x_change=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex)*2*M*L);//
	cudaMemcpy(h_x_change,d_x_change,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyDeviceToHost);
	for(int i=0;i<1;++i)//��ʱ�䣬��
	{
		for(int j=245;j<265;++j)//��ʱ�䣬��
		{
			printf("h_x_change[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,h_x_change[i*L+j].x,h_x_change[i*L+j].y);
		}
		
	}*/
	////////////////////////////////��ʼ����FFTʵ��CZT/////////////////////////////
	//x(:,0:M-1).*ww(:,M:2*M-1)
	////���һ��//////////////////��GPU����Ҫ�����õ��������ٿռ�/////////////
	/*cudaFree(d_x);
	cufftDoubleComplex *d_y;//x.*ww���
	cudaStatus=cudaMalloc((void**)&d_y,sizeof(cufftDoubleComplex)*2*M*L);//��GPU�˿���һ�齫�źŽ�����ά��ԭ����M*Lת�ü���ά��L*M*2������ѹƵ��ֵ��d_x
	//printf("һ��ϵ����Ҫ�Ĵ洢�ռ�Ϊ%dMB\n",sizeof(cufftDoubleComplex)*2*M*L/1024/1024);
	if (cudaStatus != cudaSuccess) {
        printf("d_y cudaMalloc failed!\n");
		return 1;
    }
	//cudaStatus=cudaMemcpy(d_y,d_x_change,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyDeviceToDevice);//
	if (cudaStatus != cudaSuccess) {
        printf("d_x->d_y cudaMemcpy failed!\n");
		return 1;
    }
	dim3 blcok_xw(L,1);
	dim3 threadPerBlock_xw(2*M, 1);
	//���x.*w�Ƿ���ȷ
	cufftDoubleComplex *h_y=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex)*2*M*L);//
	cudaStatus=cudaMemcpy(h_y,d_y,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyDeviceToHost);//
	if (cudaStatus != cudaSuccess) {
        printf("d_y->h_y cudaMemcpy failed!\n");
		return 1;
    }*/
	/*for(int i=0;i<1;++i)//��ʱ�䣬��
	{
		for(int j=245;j<265;++j)//��ʱ�䣬��
		{
			printf("h_y[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,h_y[i*L+j].x,h_y[i*L+j].y);
		}
	}*/

	/*cufftDoubleComplex *d_fy,*d_fv,*d_ify,*d_sp,*d_sp_change;//x.*ww��fft���
	cudaStatus=cudaMalloc((void**)&d_fy,sizeof(cufftDoubleComplex)*2*M*L);
	if (cudaStatus != cudaSuccess) {
        printf("d_fy cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_fv,sizeof(cufftDoubleComplex)*2*M*L);
	if (cudaStatus != cudaSuccess) {
        printf("d_fv cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_ify,sizeof(cufftDoubleComplex)*2*M*L);
	if (cudaStatus != cudaSuccess) {
        printf("d_ify cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_sp,sizeof(cufftDoubleComplex)*M*L);
	if (cudaStatus != cudaSuccess) {
        printf("d_sp cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_sp_change,sizeof(cufftDoubleComplex)*M*L);
	if (cudaStatus != cudaSuccess) {
        printf("d_sp_change cudaMalloc failed!\n");
		return 1;
    }
	////////fft��Plan�Ͳ���׼��
	cufftHandle Plan_fy;//��y��fft
	cufftPlan1d(&Plan_fy,2*M,CUFFT_Z2Z,L);
	

	cufftHandle Plan_ify;//��y��ifft
	cufftPlan1d(&Plan_ify,2*M,CUFFT_Z2Z,L);
	cufftHandle Plan_ifsp;//��sp��ifft
	cufftPlan1d(&Plan_ifsp,L,CUFFT_Z2Z,M);
	dim3 blcok_czt(1,L);//ÿ������һ�����뵥Ԫ
	dim3 threadPerBlock_czt(2*M,1);//ÿ�����е�һ���߳���һ������
	dim3 blcok_sp(1,M);//ÿ�����ǿ�ʱ��
	dim3 threadPerBlock_sp(L,1);//ÿ�����е�һ���߳���һ������

	

	cudaEvent_t start_event_diancheng,stop_event_diancheng,start_e_czt,stop_e_czt;
	cudaEvent_t start_e_fft,stop_e_fft;
	float time_e_diancheng,time_e_czt;
	float time_e_fft;
	cudaEventCreate(&start_e_czt);
	cudaEventCreate(&stop_e_czt);

	cudaEventCreate(&start_e_fft);
	cudaEventCreate(&stop_e_fft);
	//QueryPerformanceCounter(&b_czt);//czt��ʱ��ʼ
	//QueryPerformanceCounter(&b_diancheng);//��˼�ʱ��ʼ
	cudaEventRecord(start_e_czt,0);
	//for (int i=0;i<1;++i)
	//{
		MulVector_xw << <blcok_xw, threadPerBlock_xw >> >(d_x,d_ww,d_y,M);//���x.*w
		//QueryPerformanceCounter(&e_diancheng);//��˼�ʱ����
		cufftExecZ2Z(Plan_fy,d_v,d_fv,CUFFT_FORWARD);
		//QueryPerformanceCounter(&b_fft);//fft��ʱ��ʼ
		cudaEventRecord(start_e_fft,0);
		cufftExecZ2Z(Plan_fy,d_y,d_fy,CUFFT_FORWARD);//fft
		cudaEventRecord(stop_e_fft,0);
		
		//cudaEventSynchronize(stop_e_fft);
		//QueryPerformanceCounter(&e_fft);//fft��ʱ����

		MulVector <<<blcok_czt,threadPerBlock_czt>>>(d_fy, d_fv, d_fy, 2*L*M);//����ˣ���֪��Ϊʲô���к���

		//QueryPerformanceCounter(&b_ifft);//ifft��ʱ��ʼ
		cufftExecZ2Z(Plan_fy,d_fy,d_ify,CUFFT_INVERSE);//ifft���껹Ҫ��Ҫ������Ϊ2*M�ĳ���
		//QueryPerformanceCounter(&e_ifft);//ifft��ʱ����

		//QueryPerformanceCounter(&b_dianchu);//�����ʱ��ʼ
		//ChuVector << <blcok_czt, threadPerBlock_czt >> >(d_ify, 2*M*L,2*M);//�����
		//QueryPerformanceCounter(&e_dianchu);//�����ʱ����

		MulVector_xw << <blcok_xw, threadPerBlock_xw >> >(d_ify,d_ww,d_sp,M);//ify(:,M:2*M-1).*w(:,M:2*M-1)

		ChangeVector<< <blcok_chang, threadPerBlock_chang >> >(d_sp,d_sp_change,M);//ת��

		cufftExecZ2Z(Plan_ifsp,d_sp_change,d_sp_change,CUFFT_INVERSE);
		cufftDoubleComplex *h_sp;//czt���ս��
		h_sp=(cufftDoubleComplex * )malloc(sizeof(cufftDoubleComplex)*M*L);
		cudaStatus=cudaMemcpy(h_sp,d_sp_change,sizeof(cufftDoubleComplex)*M*L,cudaMemcpyDeviceToHost);//
		if (cudaStatus != cudaSuccess) {
			printf("d_sp_change->h_sp cudaMemcpy failed!\n");
			return 1;
		 }
		//////�ı����/////////////////////////////////
	double *h_abs_sp=(double*)malloc(sizeof(double)*M*L);
	for(int i=0; i<M*L; ++i)//SP*M
	{
		h_abs_sp[i]=sqrt((h_sp[i].x)*(h_sp[i].x)+(h_sp[i].y)*(h_sp[i].y));
		if(h_abs_sp[i]>M*Fs*Tao||h_abs_sp[i]<0) h_abs_sp[i]=0;
	}
	FILE *fp_sp;
	fp_sp=fopen("d:/sp.txt","w");
	for(int i=0; i<M; ++i)//SP*M
	{
		for(int j=0; j<L; ++j)//SP*M
		{
			fprintf(fp_sp,"%f\t",h_abs_sp[j+i*L]);
		}
		fprintf(fp_sp,"\n");
	}
	/////////////////////////////////////////////
		//ChuVector << <blcok_sp, threadPerBlock_sp >> >(d_sp_change, M*L,L);//�����
	//}
	
	cudaEventRecord(stop_e_czt,0);
	cudaEventSynchronize(stop_e_czt);
	cudaEventElapsedTime(&time_e_czt,start_e_czt,stop_e_czt);
	cudaEventElapsedTime(&time_e_fft,start_e_fft,stop_e_fft);
	cudaEventDestroy(start_e_czt);
	cudaEventDestroy(stop_e_czt);
	QueryPerformanceCounter(&e_czt);//czt��ʱ����

	//time_czt = (double(e_czt.QuadPart-b_czt.QuadPart))/(double)(fp_cpu.QuadPart);
	//time_diancheng = (double(e_diancheng.QuadPart-b_diancheng.QuadPart))/(double)(fp_cpu.QuadPart);
	//time_dianchu = (double(e_dianchu.QuadPart-b_dianchu.QuadPart))/(double)(fp_cpu.QuadPart);
	//time_ifft = (double(e_ifft.QuadPart-b_ifft.QuadPart))/(double)(fp_cpu.QuadPart);
	//time_fft = (double(e_fft.QuadPart-b_fft.QuadPart))/(double)(fp_cpu.QuadPart);
	//printf("һ��ת��ʱ��=%0.8fs,CZTʱ��=%0.8fs\n",time_czt_change,time_czt);
	//printf("һ�ε��ʱ��=%0.8fs\n",time_diancheng);
	//printf("һ�ε��ʱ��=%0.8fs\n",time_dianchu);
	//printf("һ��fftʱ��=%0.8fs\n",time_fft);
	//printf("һ��ifftʱ��=%0.8fs\n",time_ifft);
	//printf("��GPU�ڹ���ʱ=%0.8fs\n",time_czt_change+time_czt);
	printf("һ��fftʱ��event��ʱ=%0.8fs\n",time_e_fft/1000);
	printf("һ��ת��ʱ��event��ʱ=%0.8fs\n",time_e_change/1000);
	printf("��GPU�ڹ���ʱevent��ʱ=%0.8fs\n",time_e_change/1000+time_e_czt/1000);
	////////////////////////////////////////////CZT_RFT_GPU����//////////////////////////////////////////////////
	/////////////////////////���CZT��ʹ�õ�ϵ��////////////////////////////////////////////////////////////////
	cudaFree(d_ww);
    cudaFree(d_v);
	cudaFree(d_fv);
	cudaFree(d_fy);
	cudaFree(d_ify);
	cudaFree(d_sp);
	cudaFree(d_sp_change);
	cudaFree(d_x_change);
	cudaFree(d_y);
	free(h_echo);
	free(h_echo_fft);
	free(h_ht);
	free(h_ht_fft);
	free(h_pc);
	free(h_pc_fft);
	free(h_v);
	free(h_ww);*/
	/////////////////////////////////////////////////GPU_RFT��ʼ//////////////////////////////////////////////
	//RFT��������Ϳ��ٿռ�
	cufftDoubleComplex *h_Gv, *d_Gv;//�����ˣ��豸�α�׼RFT�������
	h_Gv=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*SP*M*L);
	cudaStatus=cudaMalloc((void**)&d_Gv,sizeof(cufftDoubleComplex)*SP*M*L);
	if (cudaStatus != cudaSuccess)
	{
		printf("d_Gv cudaMalloc fail!\n");
		return 1;
	}
	double DataQ=(double)sizeof(cufftDoubleComplex)*SP*M*L/1024/1024;//RFT���������
	printf("������%fMB\n\n",DataQ);
	dim3 block_s(M,SP);
	dim3 threadPerBlock_s(L,1);
	
	LARGE_INTEGER b1,b2,e1,e2;//��ʼʱ�䣬��������ʱ�䣬��������ʱ��
	double time1,time2;
	cudaEvent_t RFT_start,RFT_end;
	float time_RFT;
	cudaEventCreate(&RFT_start);
	cudaEventCreate(&RFT_end);
	cudaEventRecord(RFT_start,0);
	//QueryPerformanceCounter(&b1);
	RFT<<<block_s,threadPerBlock_s>>>(d_pc,d_Gv,Vi,L,SP,Tr,delt_R,lamda,d_DFT,d_offset);
	cudaEventRecord(RFT_end,0);
	cudaEventSynchronize(RFT_end);
	cudaEventElapsedTime(&time_RFT,RFT_start,RFT_end);
	printf("GPU-RFT����ʱ��event��ʱ:%0.8f s\n",time_RFT/1000);
	cudaEventDestroy(RFT_start);
	cudaEventDestroy(RFT_end);
	//QueryPerformanceCounter(&e1);
	//time1=(double)(e1.QuadPart-b1.QuadPart)/(double)fp_cpu.QuadPart;
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
	cudaStatus=cudaMemcpy(h_Gv,d_Gv,sizeof(cufftDoubleComplex)*SP*M*L,cudaMemcpyDeviceToHost);// RFT���,���������˵��豸��
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
	//time2=(double)(e2.QuadPart-b2.QuadPart)/(double)fp_cpu.QuadPart;
	//printf("���ݴ���ʱ��:%0.8f s\n",time2);
	//printf("����ʱ:%0.8f s\n",(time2+time1));
	double TransSpeed =DataQ/time_RFT_trans;
	//int TransSpeed =DataQ/time2;
	printf("�����ٶ�Ϊ%.2fMB s\n\n",TransSpeed*1000);

	//double Speedup=time_cpu/(time_RFT_trans+time_RFT)*1000;
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

