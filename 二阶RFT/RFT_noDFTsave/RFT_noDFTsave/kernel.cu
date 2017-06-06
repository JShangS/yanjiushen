//调用事先存储的DFT系数，通过注释来改变两种DFT存储方式
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
#define M 256//脉冲积累数
///kernel 函数
//点乘
__global__ void MulVector(cufftComplex *a, cufftComplex *b, cufftComplex *c,int size)
{
	int id_by = blockIdx.y;//所在的块号纵轴,第几个块
	int id_tx = threadIdx.x;//所在块的线程号横轴，第几个距离单元
	int id_ty = threadIdx.y;//所在块的线程号纵轴，块内第几个脉冲
	int threadPerBlock_x = blockDim.x;//块内尺寸的横轴长度，即L
	int threadPerBlock_y = blockDim.y;//块内尺寸的纵轴长度，即N

	int index_real = (id_tx+threadPerBlock_x*id_ty)+id_by*threadPerBlock_x*threadPerBlock_y; //真正的索引值=在块内+跨过的块数
	int index_ht = id_tx;//传递函数索引号
	int index_echo = index_real;// 回波函数索引号
	if (index_real < size)
	{
		c[index_real].x = a[index_ht].x * b[index_echo].x - a[index_ht].y * b[index_echo].y;//blockDim.x;//index_ht;//
		c[index_real].y = a[index_ht].x * b[index_echo].y + a[index_ht].y * b[index_echo].x;//blockDim.y;//index_real;//
	}
}
//点除
__global__ void ChuVector(cufftComplex *a, int size,int L)
{
    int id_by = blockIdx.y;//所在的块号纵轴,第几个块
	int id_tx = threadIdx.x;//所在块的线程号横轴，第几个距离单元
	int id_ty = threadIdx.y;//所在块的线程号纵轴，块内第几个脉冲
	int threadPerBlock_x = blockDim.x;//块内尺寸的横轴长度，即L
	int threadPerBlock_y = blockDim.y;//块内尺寸的纵轴长度，即N
	int index_real = (id_tx+threadPerBlock_x*id_ty)+id_by*threadPerBlock_x*threadPerBlock_y; //真正的索引值=在块内+跨过的块数
	if (index_real < size)
	{
		a[index_real].x = a[index_real].x / L;
		a[index_real].y = a[index_real].y / L;
	}
}
//标准RFT///////////////////////////////////////////////////////
//                        输入回波，              RFT结果，     搜索速度步长，速度搜索的距离单元最大值,脉冲重复间隔,距离单元大小，波长
__global__ void RFT(cufftComplex *pc, cufftComplex *Gv, float Vi, int L, int SP, float Tr , float delt_R, float lamda,cufftComplex *d_DFT, float *d_offset, float V_offset)
{
	cufftComplex Sum={0,0};
	int id_bx = blockIdx.x;//所在的块号横轴,第几个块
	int id_by = blockIdx.y;//所在的块号纵轴,第几个块
	int id_tx = threadIdx.x;//所在块的线程号横轴,初始距离单元
	int id_ty = threadIdx.y;//所在块的线程号纵轴,块内第几个速度
	int threadPerBlock_x = blockDim.x;//块内尺寸的横轴长,初始距离
	int threadPerBlock_y = blockDim.y;//块内尺寸的纵轴长度,一个块可以搜的速度个数
	int BlockPerGid_x= gridDim.x;//grid中的块坐标横轴
	int BlockPerGid_y= gridDim.y;//grid中的块坐标横轴
	//__shared__ cufftComplex Pc_share[M];
	//真正的索引值=      在块内                      +           跨过的块数        *           块的大小
	int index_real =(id_tx + id_ty*threadPerBlock_x) + (id_bx + id_by*BlockPerGid_x) * (threadPerBlock_x*threadPerBlock_y);
	int V_index = id_ty + threadPerBlock_y * (id_bx + id_by*BlockPerGid_x) ;//速度索引
	float V=Vi*V_index-V_offset;//该线程搜索的速度
	int Strat_R=id_tx;//初始距离单元
	int maxoffset=floor(V*M*Tr/delt_R+0.5);//最大距离走动
	float Sum_x=0;
	float Sum_y=0;
	int   offset=0;
	if (Strat_R-maxoffset>=0 && index_real<SP*M*L && Strat_R<L)//没超出最大搜索距离就计算
	{
		float fd = 2*V/lamda;//没存DFT系数
		for(int i=0;i<M;++i)
		{
			offset=floor(V*i*Tr/delt_R+0.5);
			Strat_R=Strat_R-d_offset[i+id_bx*M+id_by*M*M];
			
			Sum_x=cos(2*PI*fd*i*Tr);//DFT系数在线程里计算
			Sum_y=sin(2*PI*fd*i*Tr);//DFT系数在线程里计算
			Sum.x+=pc[Strat_R+i*L].x * Sum_x - pc[Strat_R+i*L].y * Sum_y;//i;//
			Sum.y+=pc[Strat_R+i*L].x * Sum_y + pc[Strat_R+i*L].y * Sum_x;//i;//
			//调用事先存储的DFT系数，通过注释来改变两种DFT存储方式
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
//转置
__global__ void ChangeVector(cufftComplex *a,cufftComplex *b,int L)//矩阵转置
{
	__shared__ cufftComplex share_block[600];//存交换用的
	//cufftComplex temp;//换用的临时变量，是在寄存器中吗？
	int id_tx = threadIdx.x;//所在块的线程号横轴
	int id_bx = blockIdx.x;//所在的块号纵轴
	int threadPerBlock_x = blockDim.x;//块内尺寸的横轴长度
	if(id_tx < L && id_bx < 2 * M)
	{
		int index_out = id_tx + id_bx*threadPerBlock_x; //输出的索引
		share_block[id_tx] = a[index_out];
	}

	__syncthreads();//标红为嘛，块内线程同步

	id_tx = threadIdx.x;//所在块的线程号横轴
	id_bx = blockIdx.x;//所在的块号纵轴
	int threadPerBlock_y = blockDim.y;//块内尺寸的纵轴长度
	if(id_tx < L && id_bx < 2 * M)
	{
		int index_in = id_tx*(2*M)+id_bx;//这里有问题
		b[index_in] = share_block[id_tx];
	}

	//temp = a[index_out];
	//a[index_out] = a[index_in];
	//a[index_in] = temp;

}
int main()
{
	//计时参数设置
	LARGE_INTEGER fp_cpu;//cpu主频
	QueryPerformanceFrequency(&fp_cpu);//获取主频
		///GPU信息说明
	int MTPB = 1024;//单block最大支持线程数
	int MaxThreadBlockSize[3]={1024,1024,64};//最大的块坐标，支持1024*1024*64个块
		/*int dev=1;
	cudaDeviceProp prop;
	cudaGetDevice(&dev);
	printf("GPU型号:%s\n",prop.name);
	printf("单block最大支持线程数:%d\n",prop.maxThreadsPerBlock);*/
	///变量声明和赋值
	float fc=100e6;//载频
	float B=4e6;//带宽
	float Tao=128e-6;//带宽
	float Fs=1*B;//采样频率
	float Ts=1/Fs;
	//printf("Ts=%0.8f\n",Ts);
	float mu=B/Tao;//调频率
	float C=3e8;//光速
	float delt_R=C/(2*Fs);//距离分辨率
	float R_start=359*delt_R;//初始距离
	float lamda=C/fc;//波长
	float PRF=500;//脉冲重复周期
	float Tr=1/PRF;//脉冲重复时间
	float Vr=2000;//初始速度
	printf("Vr=%f\n",Vr);
	float a=0;//加速度
	int    L=Tao*Fs;//距离单元数，距离采样点数
	printf("L=%d\n",L);
	int    N=MTPB/L;//块中线程纵轴维数；
	printf("N=%d\n",N);
	float *t=(float*)malloc(sizeof(float)*L);//快时间
	float *delt_t=(float*)malloc(sizeof(float)*M);//延迟时间


	for(int i=0;i<L;++i)
	{
		t[i]=-Tao/2+Ts*i;	
		//t[i]=Ts*i;	
	   // printf("t[%d]=%0.8f\n",i+1,t[i]);
	}
	
	///变量声明/////////////////////////////////////////////
	
	cufftComplex  *h_echo,*d_echo; //主机端,设备段回波时域
	cufftComplex  *h_echo_fft,*d_echo_fft; //主机端,设备段回波频域

	cufftComplex  *h_ht,*d_ht;    //主机端,设备段脉压时域系数
	cufftComplex  *h_ht_fft,*d_ht_fft;    //主机端,设备段脉压时域系数

	cufftComplex  *h_pc,*d_pc;   //主机端,设备脉压时域结果
	cufftComplex  *h_pc_fft,*d_pc_fft;//主机端,设备脉压频域结果

	//开辟内存，显存////////////////////////
	cudaError_t cudaStatus;//状态记录
    h_echo=(cufftComplex *)malloc(sizeof(cufftComplex )*M*L);//主机端回波时域
	h_echo_fft=(cufftComplex *)malloc(sizeof(cufftComplex )*M*L);//主机端回波频域

	h_ht=(cufftComplex *)malloc(sizeof(cufftComplex )*L);//主机端脉压时域系数
	h_ht_fft=(cufftComplex *)malloc(sizeof(cufftComplex )*L);//主机端脉压频域系数

	h_pc=(cufftComplex *)malloc(sizeof(cufftComplex )*M*L);//主机端脉压时域结果
	h_pc_fft=(cufftComplex *)malloc(sizeof(cufftComplex )*M*L);//主机端脉压频域结果

    cudaStatus=cudaMalloc((void**)&d_echo,sizeof(cufftComplex )*M*L);//设备段回波时域开辟显存	    
    if (cudaStatus != cudaSuccess) {
        printf( "d_echo cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_echo_fft,sizeof(cufftComplex )*M*L);//设备段回波频域开辟显存
	if (cudaStatus != cudaSuccess) {
        printf( "d_echo_fft cudaMalloc failed!\n");
		return 1;
    }

	cudaStatus=cudaMalloc((void**)&d_pc_fft,sizeof(cufftComplex )*M*L);//设备段回波脉压频域
	if (cudaStatus != cudaSuccess) {
        printf( "d_pc_fft cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_pc,sizeof(cufftComplex )*M*L);//设备段回波脉压时域
	if (cudaStatus != cudaSuccess) {
        printf( "d_pc cudaMalloc failed!\n");
		return 1;
    }

	cudaStatus=cudaMalloc((void**)&d_ht,sizeof(cufftComplex )*L);//设备段回波脉压时域系数显存
	
	if (cudaStatus != cudaSuccess) {
        printf( "d_ht cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_ht_fft,sizeof(cufftComplex )*L);//设备段回波脉压频域系数显存
	if (cudaStatus != cudaSuccess) {
        printf( "d_ht_fft cudaMalloc failed!\n");
		return 1;
    }
	
	///初始赋值///////
	///主机端脉压时域系数
	for(int i=0;i<=L-1;++i)///
	{
		h_ht[i].x=cos(2*PI*(mu/2)*t[i]*t[i]);
		h_ht[i].y=sin(2*PI*(mu/2)*t[i]*t[i]);
	}
	//传递函数翻转共轭
	/*cufftComplex t_h;//交换用的临时变量
	for (int i = 0; i < L / 2; ++i)//传递函数
	{
		if (L - i - 1 >= 0)
		{
			t_h.x = h_ht[i].x;
			h_ht[i].x = h_ht[L - i - 1].x;
			h_ht[L - i - 1].x = t_h.x;
			t_h.y = -1 * h_ht[i].y;//共轭加负号
			h_ht[i].y = -1 * h_ht[L - i - 1].y;
			h_ht[L - i - 1].y = t_h.y;
		}
	}*/
   //for(int i=0;i<=L-1;++i)printf("h_ht[%d]=%0.8f+%0.8fi\n",i+1,h_ht[i].x,h_ht[i].y);
	///回波赋值
	for(int j=0;j<M;++j)//延迟计算
	{
			Vr=Vr+a*Tr*j;
		    delt_t[j]=2*(R_start+Vr*Tr*j)/C;
	}
	for(int i=0;i<M;++i)//慢时间
	{
		for(int j=0;j<L; ++j)//快时间
		{
			h_echo[j+i*L].x=cos(2*PI*(mu/2*(t[j]+delt_t[i])*(t[j]+delt_t[i])+fc*delt_t[i]));
			h_echo[j+i*L].y=-sin(2*PI*(mu/2*(t[j]+delt_t[i])*(t[j]+delt_t[i])+fc*delt_t[i]));
			//printf("h_echo[%d][%d]=%0.8f+%0.8fi  ",i+1,j+1,h_echo[j+i*L].x,h_echo[j+i*L].y);
			//if((j+1)%3==0) printf("\n");
		}
		//printf("\n");
	}
	//for(int i=L;i<=2*L-1;++i)printf("h_echo[%d]=%0.8f+%0.8fi\n",i+1,h_echo[i].x,h_echo[i].y);
	//主机端数据拷贝到设备段
	cudaStatus=cudaMemcpy(d_echo,h_echo,sizeof(cufftComplex )*M*L,cudaMemcpyHostToDevice);//回拨数据主机端到设备段
	if (cudaStatus != cudaSuccess) {
        printf("h_echo->d_echo cudaMemcpy failed!\n");
		return 1;
    }
	cudaStatus=cudaMemcpy(d_ht,h_ht,sizeof(cufftComplex )*L,cudaMemcpyHostToDevice);// 脉压时域系数,数据主机端到设备段
    if (cudaStatus != cudaSuccess) {
        printf("h_ht->d_ht cudaMemcpy failed!\n");
		return 1;
    }
	//对回波，脉压系数做fft变换到频域，在频域进行脉压
	//设置fft计划
	LARGE_INTEGER b_pc1,e_pc1,b_pc2,e_pc2;//脉压计时
	float time_pc1,time_pc2;
	cufftHandle plan_ML,plan_L;//做L点M批的fft，做L点1批的fft计划
	cufftResult Result_fft_ML,Result_fft_L;//检查执行结果，返回值类型还不是cudaError_t真奇了怪了
	cufftPlan1d(&plan_ML,L,CUFFT_C2C,M);//设置回波fft计划
	QueryPerformanceCounter(&b_pc1);
	Result_fft_ML=cufftExecC2C(plan_ML,d_echo,d_echo_fft,CUFFT_FORWARD);//对回波进行FFT
	//此FFT是0到fs千万别和matlab中fftshift后的数据比对
	/*if (Result_fft_ML != cudaSuccess)
	{
		printf("回波fft错误!\n");	
		return 1;
	}*/
	/*cudaStatus=cudaMemcpy(h_echo_fft,d_echo,sizeof(cufftComplex )*M*L,cudaMemcpyDeviceToHost);//回拨频域数据设备段到主机端
	if (cudaStatus != cudaSuccess) {
        printf("d_echo_fft->h_echo_fft cudaMemcpy failed!\n");
		return 1;
    }*/
	///
	/*for(int i=0;i<M;++i)//慢时间
	{
		for(int j=0;j<L; ++j)//快时间
		{
			printf("h_echo_fft[%d][%d]=%0.8f+%0.8fi  ",i+1,j+1,h_echo_fft[j+i*L].x,h_echo_fft[j+i*L].y);
			if((j+1)%2==0) printf("\n");
		}
		printf("\n");
	}*/
	//for(int i=0;i<=L-1;++i)printf("h_echo_fft[%d]=%0.8f+%0.8fi\n",i+1,h_echo_fft[i].x,h_echo_fft[i].y);
	cufftPlan1d(&plan_L,L,CUFFT_C2C,1);
	Result_fft_L=cufftExecC2C(plan_L,d_ht,d_ht_fft,CUFFT_FORWARD);//对脉压系数进行FFT CUFFT_INVERSE CUFFT_FORWARD
	if (Result_fft_L != cudaSuccess)
	{
		printf("脉压系数fft错误!\n");	
		return 1;
	}
	/*cudaStatus=cudaMemcpy(h_ht_fft,d_ht_fft,sizeof(cufftComplex )*L,cudaMemcpyDeviceToHost);//脉压系数频域数据设备段到主机端
	if (cudaStatus != cudaSuccess) {
        printf("d_ht_fft->h_ht_fft cudaMemcpy failed!\n");
		return 1;
    }
	for(int j=0;j<L; ++j)
	{
			printf("h_ht_fft[%d]=%0.8f+%0.8fi \n ",j+1,h_ht_fft[j].x,h_ht_fft[j].y);
	}*/
	//进行脉压/////////////////////////////////////////////
	//GPU点乘
	dim3 blcok(1,M);//每个块是一个距离单元
	dim3 threadPerBlock(L,1);//每个块中的一个线程是一个脉冲
	MulVector <<<blcok,threadPerBlock>>>(d_ht_fft, d_echo_fft, d_pc_fft, L*M);//做点乘，不知道为什么会有红线
	//检查点乘是否正确/////////////////////
	/*cudaStatus=cudaMemcpy(h_pc_fft,d_pc_fft,sizeof(cufftComplex )*M*L,cudaMemcpyDeviceToHost);//脉压频域数据设备段到主机端
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

	Result_fft_ML=cufftExecC2C(plan_ML,d_pc_fft,d_pc,CUFFT_INVERSE);//脉压IFFT
	if (Result_fft_ML != cudaSuccess)
	{
		printf("回波fft错误!\n");	
		return 1;
	}
	ChuVector << <blcok, threadPerBlock >> >(d_pc, M*L,L);//做点除
	QueryPerformanceCounter(&e_pc1);
	time_pc1=(float)(e_pc1.QuadPart-b_pc1.QuadPart)/(float)fp_cpu.QuadPart;
	QueryPerformanceCounter(&b_pc2);
	cudaStatus=cudaMemcpy(h_pc,d_pc,sizeof(cufftComplex )*M*L,cudaMemcpyDeviceToHost);//脉压频域数据设备段到主机端
	QueryPerformanceCounter(&e_pc2);
	time_pc2=(float)(e_pc2.QuadPart-b_pc2.QuadPart)/(float)fp_cpu.QuadPart;
	printf("脉压用时:%f s\n",time_pc1);
	printf("数据传输用时:%f s\n",time_pc2);
	printf("共用时:%f s\n\n",time_pc1+time_pc2);
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
	/////////脉压数据输出////////////////////////////////////
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
	//////////////////////////////////////////脉压完成////////////////////////////////////////////////
	///////////////////////////////////////////RFT///////////////////////////////////////////////////////
	//设置搜索步长
	float Vb = lamda/2/Tr;//盲速
	//printf("Vb=%f\n",Vb);
	float Vi = Vb/M; //速度搜索步长,可调
	//printf("Vi=%f\n",Vi);
	int    SP =400;//1024*6/M;   //搜索的盲速区间个数//
	float V_offset =SP/2*Vb;
	printf("SP=%d\n",SP);
	printf("速度区间[%f,%f]\n",-SP/2*Vb,SP/2*Vb);

	/////根据SP将所有速度、对应的多普勒存储好，在计算对应的DFT系数再存储好//////////////
	float *h_fd, *h_offset, *d_offset;
	h_fd=(float *)malloc(sizeof(float )*M*SP);//主机端fd存储
	h_offset=(float *)malloc(sizeof(float )*M*SP*M);//主机端没个速度下,每个脉冲重复周期的偏移量。
	cudaStatus=cudaMalloc((void**)&d_offset,sizeof(cuComplex )*M*SP*M);//设备端没个速度下,每个脉冲重复周期的偏移量。
	if (cudaStatus != cudaSuccess) {
        printf("d_offset cudaMalloc failed!\n");
		return 1;
    }
	////生成fd,偏移量//////////
	for(int i=0;i<M*SP;++i)
	{
		h_fd[i] = 2*Vi*i/lamda;
		for(int j=0; j<M; ++j)
		{
			h_offset[i*M+j]=floor(Vi*i*j*Tr/delt_R+0.5);
		}
		
	}
	cudaStatus=cudaMemcpy(d_offset,h_offset,sizeof(float)*M*SP*M,cudaMemcpyHostToDevice);//传递数据主机端偏移量到设备段
		if (cudaStatus != cudaSuccess) {
        printf("h_offset->d_offset  cudaMemcpy failed!\n");
		return 1;
    }
	///////////////生成DFT系数///////////////
	cuComplex *h_DFT,*d_DFT;
	h_DFT=(cuComplex *)malloc(sizeof(cuComplex)*M*SP*M);//主机端DFT存储
	cudaStatus=cudaMalloc((void**)&d_DFT,sizeof(cuComplex )*M*SP*M);//设备端DFT存储
	for(int i=0; i<M*SP; ++i)
	{
		for(int j=0; j<M; ++j)
		{
			h_DFT[i*M+j].x=cos(2*PI*h_fd[i]*j*Tr);
			h_DFT[i*M+j].y=sin(2*PI*h_fd[i]*j*Tr);
		}	
	}
	cudaStatus=cudaMemcpy(d_DFT,h_DFT,sizeof(cufftComplex)*M*SP*M,cudaMemcpyHostToDevice);//传递数据主机端DFT系数到设备段
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
	cufftComplex *Gv_cpu=(cufftComplex*)malloc(sizeof(cufftComplex)*SP*M*L);//V*L,速度-距离单元二维
	LARGE_INTEGER b_cpu,e_cpu;//计时
	float time_cpu;
	float fd_cpu;
	float *V=(float*)malloc(sizeof(float)*SP*M);
	int    index_real;
	QueryPerformanceCounter(&b_cpu);
	for(int i=0; i<SP*M; ++i)//速度
	{
		V[i]=Vi*i-V_offset;
		fd_cpu=2*V[i]/lamda;//多普勒
		//printf("V=%f\t",V[i]);
		//printf("fd_cpu=%f\t",fd_cpu);
		for(int j=0; j<L; ++j)//距离单元
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
					//Gv_cpu[j+i*L].y+=h_pc[Strat_R+ti*L].x * Sum_y + h_pc[Strat_R+ti*L].y * Sum_x;//i;//不能这么累加，数据量大啦就废了
				}	
			}
			Gv_cpu[j+i*L]=Sum_cpu;
		}	
	}

	QueryPerformanceCounter(&e_cpu);
	time_cpu=(float(e_cpu.QuadPart-b_cpu.QuadPart))/(float)(fp_cpu.QuadPart);
	printf("CPU-RFT运算时间:%f s\n\n",time_cpu);
	fd_cpu=2*Vi*(SP*M-1)/lamda;*/
	//for(int i=0; i<M;++i)printf("exp=%f+%fi\n",cos(2*PI*fd_cpu*i*Tr),sin(2*PI*fd_cpu*i*Tr));////fd
	//for(int i=0; i<SP*M;++i)printf("V=%0.8f\n",V[i]);

	//for(int i=0; i<SP*M;++i)printf("Maxoffset=%0.8f\n",0 -int( V[i]*Tr*M/delt_R));
	/*for(int i=0; i<1; ++i)//SP*M
	{
		for(int j=0; j<L; ++j)printf("Gv_cpu[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,Gv_cpu[j+i*L].x,Gv_cpu[j+i*L].y);
	}*/
	////////主机端RFT文本输出/////////////////////////////////////////
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
	////////////////////////////////清除，只保留脉压后的时域和频域结果
    cudaFree(d_echo);
    cudaFree(d_echo_fft);
	cudaFree(d_ht);
	cudaFree(d_ht_fft);
	/////////////////////////////////////////////////GPU_RFT开始//////////////////////////////////////////////
	//RFT结果声明和开辟空间
	cufftComplex *h_Gv, *d_Gv;//主机端，设备段标准RFT结果声明
	h_Gv=(cufftComplex*)malloc(sizeof(cufftComplex)*SP*M*L);
	cudaStatus=cudaMalloc((void**)&d_Gv,sizeof(cufftComplex)*SP*M*L);
	if (cudaStatus != cudaSuccess)
	{
		printf("d_Gv cudaMalloc fail!\n");
		return 1;
	}
	float DataQ=(float)sizeof(cufftComplex)*SP*M*L/1024/1024;//RFT结果数据量
	printf("数据量%fMB\n\n",DataQ);
	dim3 block_s(M,SP);//M,SP
	dim3 threadPerBlock_s(L,1);
	
	LARGE_INTEGER b1,b2,e1,e2;//开始时间，结束计算时间，结束传输时间
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
	printf("GPU-RFT运算时间event计时:%0.8f s\n",time_RFT/1000);
	cudaEventDestroy(RFT_start);
	cudaEventDestroy(RFT_end);
	//QueryPerformanceCounter(&e1);
	//time1=(float)(e1.QuadPart-b1.QuadPart)/(float)fp_cpu.QuadPart;
	//printf("GPU-RFT运算时间:%0.8f s\n",time1);

	//要先释放掉否则传输很慢 
	cudaFree(d_pc);//
	cudaFree(d_pc_fft);
	////
	//QueryPerformanceCounter(&b2);
	cudaEvent_t RFT_t_start,RFT_t_end;//传输
	float time_RFT_trans;//传输时间
	cudaEventCreate(&RFT_t_start);
	cudaEventCreate(&RFT_t_end);
	cudaEventRecord(RFT_t_start,0);
	cudaStatus=cudaMemcpy(h_Gv,d_Gv,sizeof(cufftComplex)*SP*M*L,cudaMemcpyDeviceToHost);// RFT结果,数据主机端到设备段
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
	printf("RFT数据传输时间event计时:%0.8f s\n",time_RFT_trans/1000);
	printf("RFT共耗时event计时:%0.8f s\n",(time_RFT_trans+time_RFT)/1000);
	//time2=(float)(e2.QuadPart-b2.QuadPart)/(float)fp_cpu.QuadPart;
	//printf("数据传输时间:%0.8f s\n",time2);
	//printf("共耗时:%0.8f s\n",(time2+time1));
	float TransSpeed =DataQ/time_RFT_trans;
	//int TransSpeed =DataQ/time2;
	printf("传输速度为%.2fMB s\n\n",TransSpeed*1000);

//	float Speedup=time_cpu/(time_RFT_trans+time_RFT)*1000;
	//printf("加速比:%f\n",Speedup);
	/*for(int i=10*M-2; i<10*M-1; ++i)//SP*M
	{
		for(int j=0; j<L; ++j)printf("h_Gv[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,h_Gv[j+i*L].x,h_Gv[j+i*L].y);
	}*/
	/////清除显存///////////////////////////////////////////////////////////
	//cudaFree(d_Gv);
	/////////////////////////////////////////////////////////////////////////////文本输出/////
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

