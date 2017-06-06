#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cufft.h>
#include <stdio.h>
#include <stdlib.h>
__global__ void DianCheng(cufftDoubleComplex *a,  cufftDoubleComplex *b, cufftDoubleComplex *c,int M, int L)//µã³ËµÄGPUº¯Êý
{
      int tx = threadIdx.x;
      int by = blockIdx.y;
      int i=by*L+tx;
      if (i<=M*L)
    {
       c[i].x = a[i].x * b[i].x-a[i].y*b[i].y;
       c[i].y = a[i].x * b[i].y+a[i].y*b[i].x;
    }   
}