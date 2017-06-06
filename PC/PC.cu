#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cufft.h>
#include <stdio.h>
#include <stdlib.h>
#include "PC.h"
#include <mex.h>
#pragma comment(lib,"cufft.lib")
__global__ void DianCheng(cufftDoubleComplex *a,  cufftDoubleComplex *b, cufftDoubleComplex *c,int M, int L)//点乘的GPU函数
{
      int tx = threadIdx.y;
      int by = blockIdx.x;
      int i=by*L+tx;
      if (i<=M*L)
    {
       c[i].x = a[tx].x * b[i].x-a[tx].y*b[i].y;
       c[i].y = a[tx].x * b[i].y+a[tx].y*b[i].x;
    }   
}
void pc(double *realh, double *imagh, double *realecho, double *imagecho, 
        double *realpc, double * imagpc, int M, int L )
{
    printf("HI1~~~~\n");
    //把Matlab中生成的数据给到c中，从c中再给GPU///////////////////////////////
    double2 *h_h;
    double2 *h_echo;
    double2 *h_pc;
    h_h=(double2*)malloc(sizeof(double2)*M*L);
    h_echo=(double2*)malloc(sizeof(double2)*L);
    h_pc=(double2*)malloc(sizeof(double2)*M*L);
    for(int i=0; i<M*L;++i)
    {
        h_h[i].x=realh[i];
        h_h[i].y=imagh[i];
        h_echo[i].x=realecho[i];
        h_echo[i].y=imagecho[i];
    }
    printf("HI~~~~\n");
    /*for(int i=0; i<L;++i)
    {
        //printf("h[%d]=%0.8f+%0.8fi\n",i+1,realh[i],imagh[i]);
        printf("h[%d]=%0.8f+%0.8fi\n",i+1,h_h[i].x,h_h[i].y);
    }*/

    //声明GPU数组/////////////////////////////////////
    cufftDoubleComplex *dev_h;//做完共轭翻转的传递函数
    cufftDoubleComplex *dev_echo;//回波
    cufftDoubleComplex *dev_pc;//脉压时域结果
    cufftDoubleComplex *dev_pcfft;//脉压频域结果
    //开辟显存////////////////////////////////////
    cudaMalloc(&dev_h, M * L * sizeof(cufftDoubleComplex));
    cudaMalloc(&dev_echo, M * L * sizeof(cufftDoubleComplex));
    cudaMalloc(&dev_pc, M * L * sizeof(cufftDoubleComplex));
    cudaMalloc(&dev_pcfft, M * L * sizeof(cufftDoubleComplex));
    //将内存数据拷贝到显存 //////////////////////////////////
    cudaMemcpy(dev_h, h_h, M * L * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_echo, h_echo, M * L * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_pc, h_pc, M * L * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    //设置fft计划
    cufftHandle plan;
    cufftPlan1d(&plan,L,CUFFT_Z2Z,M);
    cufftExecZ2Z(plan,dev_h,dev_h,CUFFT_FORWARD);//传递函数频域
    cufftExecZ2Z(plan,dev_echo,dev_echo,CUFFT_FORWARD);
    //点乘，频移脉压
    dim3 BlockPerGrid(1,M);//网格包含的block的维度
    dim3 threadPerBlock(1,L);//每个block包含的线程维度
    DianCheng<<<BlockPerGrid,threadPerBlock>>>(dev_h,dev_echo,dev_pcfft,M,L);
    //逆傅里叶时域脉压结果
    cufftExecZ2Z(plan,dev_pcfft,dev_pc,1);//ifft
    cudaMemcpy(h_pc, dev_pc, M * L * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
    for(int i=0; i<M*L;++i)
    {
        realpc[i]=h_pc[i].x/200.0;
        imagpc[i]=h_pc[i].y/200.0;
    }
    //realpc=&(h_pc->x);
    //imagpc=&(h_pc->y);
    //释放资源
	free(h_h);
	free(h_echo);
	free(h_pc);
	cudaFree(dev_h);
	cudaFree(dev_echo);
	cudaFree(dev_pc);
	cudaFree(dev_pcfft);
}

