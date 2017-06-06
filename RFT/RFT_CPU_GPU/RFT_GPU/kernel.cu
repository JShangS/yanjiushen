
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
#define M 512//脉冲积累数
///kernel 函数
//点乘
__global__ void MulVector(cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c,int size)
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
//x.*ww
__global__ void MulVector_xw(cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c,int size)
{
	int id_by = blockIdx.y;//所在的块号纵轴,第几个块,第几个快时间距离单元
	int id_bx = blockIdx.x;//块的横轴表示第几个模糊因子
	int id_tx = threadIdx.x;//所在块的线程号横轴，第几个距离单元
	int threadPerBlock_x = blockDim.x;//块内尺寸的横轴长度，即2*M
	int threadPerBlock_y = blockDim.y;//块内尺寸的纵轴长度，即1
	int index_x=id_tx+id_by*threadPerBlock_x;//x要相乘的索引值
	int index_y=index_x;//y结果要存入的索引值
	int index_w=id_tx+size-1+id_by*threadPerBlock_x;//w要相乘的索引值
	if (id_tx < size)
	{
		c[index_y].x = a[index_x].x * b[index_w].x - a[index_x].y * b[index_w].y;//blockDim.x;//index_ht;//
		c[index_y].y = a[index_x].x * b[index_w].y + a[index_x].y * b[index_w].x;//blockDim.y;//index_real;//
	}
}
//点除
__global__ void ChuVector(cufftDoubleComplex *a, int size,int L)
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
__global__ void RFT(cufftDoubleComplex *pc, cufftDoubleComplex *Gv, double Vi, int L, int SP, double Tr , double delt_R, double lamda,cufftDoubleComplex *d_DFT, double *d_offset)
{
	cufftDoubleComplex Sum={0,0};
	int id_bx = blockIdx.x;//所在的块号横轴,第几个块
	int id_by = blockIdx.y;//所在的块号纵轴,第几个块
	int id_tx = threadIdx.x;//所在块的线程号横轴,初始距离单元
	int id_ty = threadIdx.y;//所在块的线程号纵轴,块内第几个速度
	int threadPerBlock_x = blockDim.x;//块内尺寸的横轴长,初始距离
	int threadPerBlock_y = blockDim.y;//块内尺寸的纵轴长度,一个块可以搜的速度个数
	int BlockPerGid_x= gridDim.x;//grid中的块坐标横轴
	int BlockPerGid_y= gridDim.y;//grid中的块坐标横轴
	//__shared__ cufftDoubleComplex Pc_share[M];
	//真正的索引值=      在块内                      +           跨过的块数        *           块的大小
	int index_real =(id_tx + id_ty*threadPerBlock_x) + (id_bx + id_by*BlockPerGid_x) * (threadPerBlock_x*threadPerBlock_y);
	int V_index = id_ty + threadPerBlock_y * (id_bx + id_by*BlockPerGid_x) ;//速度索引
	double V=Vi*V_index;//该线程搜索的速度
	int Strat_R=id_tx;//初始距离单元
	int maxoffset=floor(V*M*Tr/delt_R+0.5);//最大距离走动
	double Sum_x=0;
	double Sum_y=0;
	int   offset=0;
	if (Strat_R-maxoffset>=0 && index_real<SP*M*L && Strat_R<L)//没超出最大搜索距离就计算
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
//转置
__global__ void ChangeVector(cufftDoubleComplex *a,cufftDoubleComplex *b,int L)//矩阵转置
{
	__shared__ cufftDoubleComplex share_block[600];//存交换用的
	//cufftDoubleComplex temp;//换用的临时变量，是在寄存器中吗？
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
	double fc=100e6;//载频
	double B=4e6;//带宽
	double Tao=128e-6;//带宽
	double Fs=2*B;//采样频率
	double Ts=1/Fs;
	//printf("Ts=%0.8f\n",Ts);
	double mu=B/Tao;//调频率
	double C=3e8;//光速
	double delt_R=C/(2*Fs);//距离分辨率
	double R_start=500*delt_R;//初始距离
	double lamda=C/fc;//波长
	double PRF=500;//脉冲重复周期
	double Tr=1/PRF;//脉冲重复时间
	double Vr=1200;//初始速度
	printf("Vr=%f\n",Vr);
	double a=0;//加速度
	int    L=Tao*Fs;//距离单元数，距离采样点数
	printf("L=%d\n",L);
	int    N=MTPB/L;//块中线程纵轴维数；
	printf("N=%d\n",N);
	double *t=(double*)malloc(sizeof(double)*L);//快时间
	double *delt_t=(double*)malloc(sizeof(double)*M);//延迟时间


	for(int i=0;i<L;++i)
	{
		t[i]=-Tao/2+Ts*i;	
		//t[i]=Ts*i;	
	   // printf("t[%d]=%0.8f\n",i+1,t[i]);
	}
	
	///变量声明/////////////////////////////////////////////
	
	cufftDoubleComplex  *h_echo,*d_echo; //主机端,设备段回波时域
	cufftDoubleComplex  *h_echo_fft,*d_echo_fft; //主机端,设备段回波频域

	cufftDoubleComplex  *h_ht,*d_ht;    //主机端,设备段脉压时域系数
	cufftDoubleComplex  *h_ht_fft,*d_ht_fft;    //主机端,设备段脉压时域系数

	cufftDoubleComplex  *h_pc,*d_pc;   //主机端,设备脉压时域结果
	cufftDoubleComplex  *h_pc_fft,*d_pc_fft;//主机端,设备脉压频域结果

	//开辟内存，显存////////////////////////
	cudaError_t cudaStatus;//状态记录
    h_echo=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*M*L);//主机端回波时域
	h_echo_fft=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*M*L);//主机端回波频域

	h_ht=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*L);//主机端脉压时域系数
	h_ht_fft=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*L);//主机端脉压频域系数

	h_pc=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*M*L);//主机端脉压时域结果
	h_pc_fft=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex )*M*L);//主机端脉压频域结果

    cudaStatus=cudaMalloc((void**)&d_echo,sizeof(cufftDoubleComplex )*M*L);//设备段回波时域开辟显存	    
    if (cudaStatus != cudaSuccess) {
        printf( "d_echo cudaMalloc failed!\n %d \n",cudaStatus);
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_echo_fft,sizeof(cufftDoubleComplex )*M*L);//设备段回波频域开辟显存
	if (cudaStatus != cudaSuccess) {
        printf( "d_echo_fft cudaMalloc failed!\n");
		return 1;
    }

	cudaStatus=cudaMalloc((void**)&d_pc_fft,sizeof(cufftDoubleComplex )*M*L);//设备段回波脉压频域
	if (cudaStatus != cudaSuccess) {
        printf( "d_pc_fft cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_pc,sizeof(cufftDoubleComplex )*M*L);//设备段回波脉压时域
	if (cudaStatus != cudaSuccess) {
        printf( "d_pc cudaMalloc failed!\n");
		return 1;
    }

	cudaStatus=cudaMalloc((void**)&d_ht,sizeof(cufftDoubleComplex )*L);//设备段回波脉压时域系数显存
	
	if (cudaStatus != cudaSuccess) {
        printf( "d_ht cudaMalloc failed! \n");
		return 1;
    }
	cudaStatus=cudaMalloc((void**)&d_ht_fft,sizeof(cufftDoubleComplex )*L);//设备段回波脉压频域系数显存
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
	/*cufftDoubleComplex t_h;//交换用的临时变量
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
	cudaStatus=cudaMemcpy(d_echo,h_echo,sizeof(cufftDoubleComplex )*M*L,cudaMemcpyHostToDevice);//回拨数据主机端到设备段
	if (cudaStatus != cudaSuccess) {
        printf("h_echo->d_echo cudaMemcpy failed!\n");
		return 1;
    }
	cudaStatus=cudaMemcpy(d_ht,h_ht,sizeof(cufftDoubleComplex )*L,cudaMemcpyHostToDevice);// 脉压时域系数,数据主机端到设备段
    if (cudaStatus != cudaSuccess) {
        printf("h_ht->d_ht cudaMemcpy failed!\n");
		return 1;
    }
	//对回波，脉压系数做fft变换到频域，在频域进行脉压
	//设置fft计划
	LARGE_INTEGER b_pc1,e_pc1,b_pc2,e_pc2;//脉压计时
	double time_pc1,time_pc2;
	
	cufftHandle plan_ML,plan_L;//做L点M批的fft，做L点1批的fft计划
	cufftResult Result_fft_ML,Result_fft_L;//检查执行结果，返回值类型还不是cudaError_t真奇了怪了
	cufftPlan1d(&plan_ML,L,CUFFT_Z2Z,M);//设置回波fft计划
	QueryPerformanceCounter(&b_pc1);
	Result_fft_ML=cufftExecZ2Z(plan_ML,d_echo,d_echo_fft,CUFFT_FORWARD);//对回波进行FFT
	//此FFT是0到fs千万别和matlab中fftshift后的数据比对
	/*if (Result_fft_ML != cudaSuccess)
	{
		printf("回波fft错误!\n");	
		return 1;
	}*/
	/*cudaStatus=cudaMemcpy(h_echo_fft,d_echo,sizeof(cufftDoubleComplex )*M*L,cudaMemcpyDeviceToHost);//回拨频域数据设备段到主机端
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
	cufftPlan1d(&plan_L,L,CUFFT_Z2Z,1);
	Result_fft_L=cufftExecZ2Z(plan_L,d_ht,d_ht_fft,CUFFT_FORWARD);//对脉压系数进行FFT CUFFT_INVERSE CUFFT_FORWARD
	if (Result_fft_L != cudaSuccess)
	{
		printf("脉压系数fft错误!\n");	
		return 1;
	}
	/*cudaStatus=cudaMemcpy(h_ht_fft,d_ht_fft,sizeof(cufftDoubleComplex )*L,cudaMemcpyDeviceToHost);//脉压系数频域数据设备段到主机端
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
	/*cudaStatus=cudaMemcpy(h_pc_fft,d_pc_fft,sizeof(cufftDoubleComplex )*M*L,cudaMemcpyDeviceToHost);//脉压频域数据设备段到主机端
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

	Result_fft_ML=cufftExecZ2Z(plan_ML,d_pc_fft,d_pc,CUFFT_INVERSE);//脉压IFFT
	if (Result_fft_ML != cudaSuccess)
	{
		printf("回波fft错误!\n");	
		return 1;
	}
	ChuVector << <blcok, threadPerBlock >> >(d_pc, M*L,L);//做点除
	QueryPerformanceCounter(&e_pc1);
	time_pc1=(double)(e_pc1.QuadPart-b_pc1.QuadPart)/(double)fp_cpu.QuadPart;
	QueryPerformanceCounter(&b_pc2);
	cudaStatus=cudaMemcpy(h_pc,d_pc,sizeof(cufftDoubleComplex )*M*L,cudaMemcpyDeviceToHost);//脉压频域数据设备段到主机端
	QueryPerformanceCounter(&e_pc2);
	time_pc2=(double)(e_pc2.QuadPart-b_pc2.QuadPart)/(double)fp_cpu.QuadPart;
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
	//////////////////////////////////////////脉压完成////////////////////////////////////////////////
	///////////////////////////////////////////RFT///////////////////////////////////////////////////////
	//设置搜索步长
	double Vb = lamda/2/Tr;//盲速
	//printf("Vb=%f\n",Vb);
	double Vi = Vb/M; //速度搜索步长
	//printf("Vi=%f\n",Vi);
	int    SP =8;//1024*6/M;   //搜索的盲速区间个数//，数据量一大就黑屏，GTX650显卡，破显卡还是破程序？
	printf("SP=%d\n",SP);


	/////根据SP将所有速度、对应的多普勒存储好，在计算对应的DFT系数再存储好//////////////
	double *h_fd, *h_offset, *d_offset;
	h_fd=(double *)malloc(sizeof(double )*M*SP);//主机端fd存储
	h_offset=(double *)malloc(sizeof(double )*M*SP*M);//主机端没个速度下,每个脉冲重复周期的偏移量。
	cudaStatus=cudaMalloc((void**)&d_offset,sizeof(cuDoubleComplex )*M*SP*M);//设备端没个速度下,每个脉冲重复周期的偏移量。
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
	cudaStatus=cudaMemcpy(d_offset,h_offset,sizeof(double)*M*SP*M,cudaMemcpyHostToDevice);//传递数据主机端偏移量到设备段
		if (cudaStatus != cudaSuccess) {
        printf("h_offset->d_offset  cudaMemcpy failed!\n");
		return 1;
    }
	///////////////生成DFT系数///////////////
	cuDoubleComplex *h_DFT,*d_DFT;
	h_DFT=(cuDoubleComplex *)malloc(sizeof(cuDoubleComplex)*M*SP*M);//主机端DFT存储
	cudaStatus=cudaMalloc((void**)&d_DFT,sizeof(cuDoubleComplex )*M*SP*M);//设备端DFT存储
	for(int i=0; i<M*SP; ++i)
	{
		for(int j=0; j<M; ++j)
		{
			h_DFT[i*M+j].x=cos(2*PI*h_fd[i]*j*Tr);
			h_DFT[i*M+j].y=sin(2*PI*h_fd[i]*j*Tr);
		}	
	}
	cudaStatus=cudaMemcpy(d_DFT,h_DFT,sizeof(cufftDoubleComplex)*M*SP*M,cudaMemcpyHostToDevice);//传递数据主机端DFT系数到设备段
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
	cufftDoubleComplex *Gv_cpu=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*SP*M*L);//V*L,速度-距离单元二维
	LARGE_INTEGER b_cpu,e_cpu;//计时
	double time_cpu;
	double fd_cpu;
	double *V=(double*)malloc(sizeof(double)*SP*M);
	int    index_real;
	QueryPerformanceCounter(&b_cpu);
	for(int i=0; i<SP*M; ++i)//速度
	{
		V[i]=Vi*i;
		fd_cpu=2*V[i]/lamda;//多普勒
		//printf("V=%f\t",V[i]);
		//printf("fd_cpu=%f\t",fd_cpu);
		for(int j=0; j<L; ++j)//距离单元
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
					//Gv_cpu[j+i*L].y+=h_pc[Strat_R+ti*L].x * Sum_y + h_pc[Strat_R+ti*L].y * Sum_x;//i;//不能这么累加，数据量大啦就废了
				}	
			}
			Gv_cpu[j+i*L]=Sum_cpu;
		}	
	}

	QueryPerformanceCounter(&e_cpu);
	time_cpu=(double(e_cpu.QuadPart-b_cpu.QuadPart))/(double)(fp_cpu.QuadPart);
	printf("CPU-RFT运算时间:%f s\n\n",time_cpu);
	fd_cpu=2*Vi*(SP*M-1)/lamda;*/
	//for(int i=0; i<M;++i)printf("exp=%f+%fi\n",cos(2*PI*fd_cpu*i*Tr),sin(2*PI*fd_cpu*i*Tr));////fd
	//for(int i=0; i<SP*M;++i)printf("V=%0.8f\n",V[i]);

	//for(int i=0; i<SP*M;++i)printf("Maxoffset=%0.8f\n",0 -int( V[i]*Tr*M/delt_R));
	/*for(int i=0; i<1; ++i)//SP*M
	{
		for(int j=0; j<L; ++j)printf("Gv_cpu[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,Gv_cpu[j+i*L].x,Gv_cpu[j+i*L].y);
	}
	////////主机端RFT文本输出/////////////////////////////////////////
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
	////////////////////////////////清除，只保留脉压后的时域和频域结果
    cudaFree(d_echo);
    cudaFree(d_echo_fft);
	cudaFree(d_ht);
	cudaFree(d_ht_fft);
	/////////////////////////////////////////////CZT_RFT_GPU开始/////////////////////////////////////////////////////
	//系数声明和赋值////
	/*//////参见matlab中czt的实现
	double fai=lamda*B/L/C;//快时间频域分辨率
	double *fai_a=(double*)malloc(sizeof(double)*L);//1-fai的结果
	//cufftDoubleComplex *w=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*L);//不同快时间频率点处的基础czt系数
	for(int i=0;i<L;++i)//同一个快时间处用一个fai_a
	{
		fai_a[i]=1.0f-fai*(i+1.0f);
		//w=exp(-1j*2*pi*fai_a/L);
		//w[i].x=cos(2*PI*fai_a[i]/M);
		//w[i].y=-sin(2*PI*fai_a[i]/M);
		//printf("fai_a[%d]=%0.8f\n",i+1,fai_a[i]);
	}
	double *kk=(double*)malloc(sizeof(double)*2*M-1);//幂次
	double *kk2=(double*)malloc(sizeof(double)*2*M-1);//幂次平方
	for(int i=0;i<2*M;++i)
	{
		kk[i]=-1.0f*M+i+1;
		kk2[i]=kk[i]*kk[i]/2.0f;
		//printf("kk2[%d]=%0.8f\n",i+1,kk2[i]);
	}
	cufftDoubleComplex *h_ww=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*L*2*M);//不同快时间频率点处的czt系数
	cufftDoubleComplex *h_v=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*L*2*M);//不同快时间频率点处的czt 1./系数
	///ww共L行，2*M列，其中每一行对应一个快时间频点，对同一个快时间频点处的M个脉冲回波做CZT变换
	for(int i=0;i<L;++i)
	{
		//w=exp(-1j*2*pi*fai_a/L);
		//ww=w.^kk2
		//v=1./ww;
		for(int j=0; j<2*M; ++j)///出来结果和matlba是转置的关系 
		{
			//int indx_a=int(i/2/M);//同一个快时间频点
			//int indx_kk2=i%(2*M);//以2*M循环
			//printf("indx_a=%d,indx_kk2=%d\n",indx_a,indx_kk2);
			h_ww[i*2*M+j].x=cos(2*PI*fai_a[i]*kk2[j]/M);
			h_ww[i*2*M+j].y=-sin(2*PI*fai_a[i]*kk2[j]/M);
			h_v[i*2*M+j].x=cos(2*PI*fai_a[i]*kk2[j]/M);
			h_v[i*2*M+j].y=sin(2*PI*fai_a[i]*kk2[j]/M);
		}

		
	}
	
	cufftDoubleComplex *d_ww;
	
	cudaStatus=cudaMalloc((void**)&d_ww,sizeof(cufftDoubleComplex)*2*M*L);//在GPU端开辟一块将信号进行扩维将原来的M*L转置加扩维成L*M*2，将脉压频域赋值给d_x
	if (cudaStatus != cudaSuccess) {
        printf("d_ww cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMemcpy(d_ww,h_ww,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyHostToDevice);//传递数据到x中作czt
	if (cudaStatus != cudaSuccess) {
        printf("h_ww->d_ww cudaMemcpy failed!\n");
		return 1;
    }
	/////换一下位置，把后面M个w值和前面M个w值位置互换//////////////
	cufftDoubleComplex *d_ww_change;
	cudaStatus=cudaMalloc((void**)&d_ww_change,sizeof(cufftDoubleComplex)*2*M*L);//在GPU端开辟一块将信号进行扩维将原来的M*L转置加扩维成L*M*2，将脉压频域赋值给d_x
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
	//检查ww是否正确,没错2016/10/11///////////////////
	/*for(int i=0;i<1;++i)//快时间，行
	{
		for(int j=255;j<275;++j)//慢时间，列
		{
			printf("v[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,h_v[i*2*M+j].x,h_v[i*2*M+j].y);
		}
		
	}*/

	/*cufftDoubleComplex *d_v;
	cudaStatus=cudaMalloc((void**)&d_v,sizeof(cufftDoubleComplex)*2*M*L);//在GPU端开辟一块将信号进行扩维将原来的M*L转置加扩维成L*M*2，将脉压频域赋值给d_x
	if (cudaStatus != cudaSuccess) {
        printf("d_v cudaMalloc failed!\n");
		return 1;
    }
	cudaStatus=cudaMemcpy(d_v,h_v,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyHostToDevice);//传递数据到x中作czt
	if (cudaStatus != cudaSuccess) {
        printf("h_v->d_v cudaMemcpy failed!\n");
		return 1;
    }
	//扩展x做FFT前的准备
	cufftDoubleComplex *zeros_x=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*1*M*L);//一个M行L列的0矩阵，给d_x后半段赋值
	for(int i = 0; i < 1 *M * L; ++i)
	{
		zeros_x[i].x = 0;
		zeros_x[i].y = 0;
		//if(i%100==0) printf("%0.8f+%0.8fi\n",zeros_x[i].x,zeros_x[i].y);

	}
	cufftDoubleComplex *d_x;//pc_fft的扩展补零
	cudaStatus=cudaMalloc((void**)&d_x,sizeof(cufftDoubleComplex)*2*M*L);//在GPU端开辟一块将信号进行扩维将原来的M*L转置加扩维成L*M*2，将脉压频域赋值给d_x
	if (cudaStatus != cudaSuccess) {
        printf("d_x cudaMalloc failed!\n");
		return 1;
    }
	//把频域脉压结果补上
	//把0补齐
	cudaStatus=cudaMemcpy(&d_x[0],d_pc_fft,sizeof(cufftDoubleComplex)*M*L,cudaMemcpyDeviceToDevice);//传递数据到x中作czt
	if (cudaStatus != cudaSuccess) {
        printf("d_pc_fft->d_x cudaMemcpy failed!\n");
		return 1;
    }//把0补在x下面
	cudaStatus=cudaMemcpy(&d_x[M*L],zeros_x,sizeof(cufftDoubleComplex)*1*M*L,cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
        printf("zeros_x->d_x cudaMemcpy failed!\n");
		return 1;
    }*/
    
	
	///////////转置//////////////////////////////////////开始CZT//////////////////////////
	/*cufftDoubleComplex *d_x_change;//转置后的d_x
	cudaMalloc((void**)&d_x_change,sizeof(cufftDoubleComplex)*L*2*M);//
	dim3 blcok_chang(2*M, 1);
	dim3 threadPerBlock_chang(L, 1);
	////////////计时变量声明////////////////////
	LARGE_INTEGER b_czt_change,e_czt_change,b_czt,e_czt,b_diancheng,e_diancheng,b_fft,e_fft;//计时,chang:转置
	LARGE_INTEGER b_dianchu,e_dianchu, b_ifft, e_ifft;
	double time_czt_change,time_czt,time_diancheng,time_fft,time_ifft,time_dianchu;
	cudaEvent_t stop_e_change,start_e_change;
	float time_e_change;
	cudaEventCreate(&start_e_change);
	cudaEventCreate(&stop_e_change);
	//QueryPerformanceCounter(&b_czt_change);
	cudaEventRecord(start_e_change,0);
	ChangeVector << <blcok_chang, threadPerBlock_chang >> >(d_x, d_x_change,L);//矩阵转置
	cudaEventRecord(stop_e_change,0);
	cudaEventSynchronize(stop_e_change);
	cudaEventElapsedTime(&time_e_change,start_e_change,stop_e_change);

	cudaEventDestroy(stop_e_change);
	cudaEventDestroy(start_e_change);
	//QueryPerformanceCounter(&e_czt_change);
	//time_czt_change =(double(e_czt_change.QuadPart-b_czt_change.QuadPart))/(double)(fp_cpu.QuadPart);*/

	///检查转置是否成功////////////////////////
	/*cufftDoubleComplex *h_x_change=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex)*2*M*L);//
	cudaMemcpy(h_x_change,d_x_change,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyDeviceToHost);
	for(int i=0;i<1;++i)//快时间，行
	{
		for(int j=245;j<265;++j)//慢时间，列
		{
			printf("h_x_change[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,h_x_change[i*L+j].x,h_x_change[i*L+j].y);
		}
		
	}*/
	////////////////////////////////开始利用FFT实现CZT/////////////////////////////
	//x(:,0:M-1).*ww(:,M:2*M-1)
	////清除一波//////////////////给GPU中需要计算用到参数开辟空间/////////////
	/*cudaFree(d_x);
	cufftDoubleComplex *d_y;//x.*ww结果
	cudaStatus=cudaMalloc((void**)&d_y,sizeof(cufftDoubleComplex)*2*M*L);//在GPU端开辟一块将信号进行扩维将原来的M*L转置加扩维成L*M*2，将脉压频域赋值给d_x
	//printf("一组系数需要的存储空间为%dMB\n",sizeof(cufftDoubleComplex)*2*M*L/1024/1024);
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
	//检查x.*w是否正确
	cufftDoubleComplex *h_y=(cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex)*2*M*L);//
	cudaStatus=cudaMemcpy(h_y,d_y,sizeof(cufftDoubleComplex)*2*M*L,cudaMemcpyDeviceToHost);//
	if (cudaStatus != cudaSuccess) {
        printf("d_y->h_y cudaMemcpy failed!\n");
		return 1;
    }*/
	/*for(int i=0;i<1;++i)//快时间，行
	{
		for(int j=245;j<265;++j)//慢时间，列
		{
			printf("h_y[%d][%d]=%0.8f+%0.8fi\n",i+1,j+1,h_y[i*L+j].x,h_y[i*L+j].y);
		}
	}*/

	/*cufftDoubleComplex *d_fy,*d_fv,*d_ify,*d_sp,*d_sp_change;//x.*ww的fft结果
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
	////////fft的Plan和参数准备
	cufftHandle Plan_fy;//对y做fft
	cufftPlan1d(&Plan_fy,2*M,CUFFT_Z2Z,L);
	

	cufftHandle Plan_ify;//对y做ifft
	cufftPlan1d(&Plan_ify,2*M,CUFFT_Z2Z,L);
	cufftHandle Plan_ifsp;//对sp做ifft
	cufftPlan1d(&Plan_ifsp,L,CUFFT_Z2Z,M);
	dim3 blcok_czt(1,L);//每个块是一个距离单元
	dim3 threadPerBlock_czt(2*M,1);//每个块中的一个线程是一个脉冲
	dim3 blcok_sp(1,M);//每个块是快时间
	dim3 threadPerBlock_sp(L,1);//每个块中的一个线程是一个脉冲

	

	cudaEvent_t start_event_diancheng,stop_event_diancheng,start_e_czt,stop_e_czt;
	cudaEvent_t start_e_fft,stop_e_fft;
	float time_e_diancheng,time_e_czt;
	float time_e_fft;
	cudaEventCreate(&start_e_czt);
	cudaEventCreate(&stop_e_czt);

	cudaEventCreate(&start_e_fft);
	cudaEventCreate(&stop_e_fft);
	//QueryPerformanceCounter(&b_czt);//czt计时开始
	//QueryPerformanceCounter(&b_diancheng);//点乘计时开始
	cudaEventRecord(start_e_czt,0);
	//for (int i=0;i<1;++i)
	//{
		MulVector_xw << <blcok_xw, threadPerBlock_xw >> >(d_x,d_ww,d_y,M);//点乘x.*w
		//QueryPerformanceCounter(&e_diancheng);//点乘计时结束
		cufftExecZ2Z(Plan_fy,d_v,d_fv,CUFFT_FORWARD);
		//QueryPerformanceCounter(&b_fft);//fft计时开始
		cudaEventRecord(start_e_fft,0);
		cufftExecZ2Z(Plan_fy,d_y,d_fy,CUFFT_FORWARD);//fft
		cudaEventRecord(stop_e_fft,0);
		
		//cudaEventSynchronize(stop_e_fft);
		//QueryPerformanceCounter(&e_fft);//fft计时结束

		MulVector <<<blcok_czt,threadPerBlock_czt>>>(d_fy, d_fv, d_fy, 2*L*M);//做点乘，不知道为什么会有红线

		//QueryPerformanceCounter(&b_ifft);//ifft计时开始
		cufftExecZ2Z(Plan_fy,d_fy,d_ify,CUFFT_INVERSE);//ifft做完还要做要做除数为2*M的除法
		//QueryPerformanceCounter(&e_ifft);//ifft计时结束

		//QueryPerformanceCounter(&b_dianchu);//点除计时开始
		//ChuVector << <blcok_czt, threadPerBlock_czt >> >(d_ify, 2*M*L,2*M);//做点除
		//QueryPerformanceCounter(&e_dianchu);//点除计时结束

		MulVector_xw << <blcok_xw, threadPerBlock_xw >> >(d_ify,d_ww,d_sp,M);//ify(:,M:2*M-1).*w(:,M:2*M-1)

		ChangeVector<< <blcok_chang, threadPerBlock_chang >> >(d_sp,d_sp_change,M);//转置

		cufftExecZ2Z(Plan_ifsp,d_sp_change,d_sp_change,CUFFT_INVERSE);
		cufftDoubleComplex *h_sp;//czt最终结果
		h_sp=(cufftDoubleComplex * )malloc(sizeof(cufftDoubleComplex)*M*L);
		cudaStatus=cudaMemcpy(h_sp,d_sp_change,sizeof(cufftDoubleComplex)*M*L,cudaMemcpyDeviceToHost);//
		if (cudaStatus != cudaSuccess) {
			printf("d_sp_change->h_sp cudaMemcpy failed!\n");
			return 1;
		 }
		//////文本输出/////////////////////////////////
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
		//ChuVector << <blcok_sp, threadPerBlock_sp >> >(d_sp_change, M*L,L);//做点除
	//}
	
	cudaEventRecord(stop_e_czt,0);
	cudaEventSynchronize(stop_e_czt);
	cudaEventElapsedTime(&time_e_czt,start_e_czt,stop_e_czt);
	cudaEventElapsedTime(&time_e_fft,start_e_fft,stop_e_fft);
	cudaEventDestroy(start_e_czt);
	cudaEventDestroy(stop_e_czt);
	QueryPerformanceCounter(&e_czt);//czt计时结束

	//time_czt = (double(e_czt.QuadPart-b_czt.QuadPart))/(double)(fp_cpu.QuadPart);
	//time_diancheng = (double(e_diancheng.QuadPart-b_diancheng.QuadPart))/(double)(fp_cpu.QuadPart);
	//time_dianchu = (double(e_dianchu.QuadPart-b_dianchu.QuadPart))/(double)(fp_cpu.QuadPart);
	//time_ifft = (double(e_ifft.QuadPart-b_ifft.QuadPart))/(double)(fp_cpu.QuadPart);
	//time_fft = (double(e_fft.QuadPart-b_fft.QuadPart))/(double)(fp_cpu.QuadPart);
	//printf("一次转置时间=%0.8fs,CZT时间=%0.8fs\n",time_czt_change,time_czt);
	//printf("一次点乘时间=%0.8fs\n",time_diancheng);
	//printf("一次点除时间=%0.8fs\n",time_dianchu);
	//printf("一次fft时间=%0.8fs\n",time_fft);
	//printf("一次ifft时间=%0.8fs\n",time_ifft);
	//printf("在GPU内共耗时=%0.8fs\n",time_czt_change+time_czt);
	printf("一次fft时间event计时=%0.8fs\n",time_e_fft/1000);
	printf("一次转置时间event计时=%0.8fs\n",time_e_change/1000);
	printf("在GPU内共耗时event计时=%0.8fs\n",time_e_change/1000+time_e_czt/1000);
	////////////////////////////////////////////CZT_RFT_GPU结束//////////////////////////////////////////////////
	/////////////////////////清除CZT中使用的系数////////////////////////////////////////////////////////////////
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
	/////////////////////////////////////////////////GPU_RFT开始//////////////////////////////////////////////
	//RFT结果声明和开辟空间
	cufftDoubleComplex *h_Gv, *d_Gv;//主机端，设备段标准RFT结果声明
	h_Gv=(cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*SP*M*L);
	cudaStatus=cudaMalloc((void**)&d_Gv,sizeof(cufftDoubleComplex)*SP*M*L);
	if (cudaStatus != cudaSuccess)
	{
		printf("d_Gv cudaMalloc fail!\n");
		return 1;
	}
	double DataQ=(double)sizeof(cufftDoubleComplex)*SP*M*L/1024/1024;//RFT结果数据量
	printf("数据量%fMB\n\n",DataQ);
	dim3 block_s(M,SP);
	dim3 threadPerBlock_s(L,1);
	
	LARGE_INTEGER b1,b2,e1,e2;//开始时间，结束计算时间，结束传输时间
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
	printf("GPU-RFT运算时间event计时:%0.8f s\n",time_RFT/1000);
	cudaEventDestroy(RFT_start);
	cudaEventDestroy(RFT_end);
	//QueryPerformanceCounter(&e1);
	//time1=(double)(e1.QuadPart-b1.QuadPart)/(double)fp_cpu.QuadPart;
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
	cudaStatus=cudaMemcpy(h_Gv,d_Gv,sizeof(cufftDoubleComplex)*SP*M*L,cudaMemcpyDeviceToHost);// RFT结果,数据主机端到设备段
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
	//time2=(double)(e2.QuadPart-b2.QuadPart)/(double)fp_cpu.QuadPart;
	//printf("数据传输时间:%0.8f s\n",time2);
	//printf("共耗时:%0.8f s\n",(time2+time1));
	double TransSpeed =DataQ/time_RFT_trans;
	//int TransSpeed =DataQ/time2;
	printf("传输速度为%.2fMB s\n\n",TransSpeed*1000);

	//double Speedup=time_cpu/(time_RFT_trans+time_RFT)*1000;
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

