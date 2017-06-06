#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cufft.h>
#include <stdio.h>
#include <stdlib.h>
#include "PC.h"
#include <mex.h>
#pragma comment(lib,"cufft.lib")
__global__ void DianCheng(cufftDoubleComplex *a,  cufftDoubleComplex *b, cufftDoubleComplex *c,int M, int L)//��˵�GPU����
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
    //��Matlab�����ɵ����ݸ���c�У���c���ٸ�GPU///////////////////////////////
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

    //����GPU����/////////////////////////////////////
    cufftDoubleComplex *dev_h;//���깲�ת�Ĵ��ݺ���
    cufftDoubleComplex *dev_echo;//�ز�
    cufftDoubleComplex *dev_pc;//��ѹʱ����
    cufftDoubleComplex *dev_pcfft;//��ѹƵ����
    //�����Դ�////////////////////////////////////
    cudaMalloc(&dev_h, M * L * sizeof(cufftDoubleComplex));
    cudaMalloc(&dev_echo, M * L * sizeof(cufftDoubleComplex));
    cudaMalloc(&dev_pc, M * L * sizeof(cufftDoubleComplex));
    cudaMalloc(&dev_pcfft, M * L * sizeof(cufftDoubleComplex));
    //���ڴ����ݿ������Դ� //////////////////////////////////
    cudaMemcpy(dev_h, h_h, M * L * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_echo, h_echo, M * L * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_pc, h_pc, M * L * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    //����fft�ƻ�
    cufftHandle plan;
    cufftPlan1d(&plan,L,CUFFT_Z2Z,M);
    cufftExecZ2Z(plan,dev_h,dev_h,CUFFT_FORWARD);//���ݺ���Ƶ��
    cufftExecZ2Z(plan,dev_echo,dev_echo,CUFFT_FORWARD);
    //��ˣ�Ƶ����ѹ
    dim3 BlockPerGrid(1,M);//���������block��ά��
    dim3 threadPerBlock(1,L);//ÿ��block�������߳�ά��
    DianCheng<<<BlockPerGrid,threadPerBlock>>>(dev_h,dev_echo,dev_pcfft,M,L);
    //�渵��Ҷʱ����ѹ���
    cufftExecZ2Z(plan,dev_pcfft,dev_pc,1);//ifft
    cudaMemcpy(h_pc, dev_pc, M * L * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
    for(int i=0; i<M*L;++i)
    {
        realpc[i]=h_pc[i].x/200.0;
        imagpc[i]=h_pc[i].y/200.0;
    }
    //realpc=&(h_pc->x);
    //imagpc=&(h_pc->y);
    //�ͷ���Դ
	free(h_h);
	free(h_echo);
	free(h_pc);
	cudaFree(dev_h);
	cudaFree(dev_echo);
	cudaFree(dev_pc);
	cudaFree(dev_pcfft);
}

