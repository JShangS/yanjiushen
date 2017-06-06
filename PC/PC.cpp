#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include "PC.h"
//#pragma comment(lib,"cufft.lib")
void mexFunction(int nlhs,mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    if (nrhs!=4)
    {
      //mexErrMsgTxt("参数输入不够");
      mexPrintf(" 参数输入为:%d \n 应该为:4,\n 要有传递函数实部、虚部，回波实部、虚部\n",nrhs);   
    }
    mexPrintf(" 参数输出为:%d\n",nlhs);
    if (nlhs!=2)   
        mexErrMsgTxt("参数输出不够,要有2个：脉压时域的实部、虚部");
    int L=mxGetM(prhs[2]);//获得行数
    int M=mxGetN(prhs[2]);//获得列数
    mexPrintf("L=%d,M=%d\n",L,M); 
    double *realh=(double*)mxGetPr(prhs[0]);//获取输入参数
   // mexPrintf("地址:realh[%d]=%d,realh[%d]=%d,realh[%d]=%d\n",0,&realh[0],1,&realh[1],2,&realh[2]); 
    double *imagh=(double*)mxGetPr(prhs[1]);
    double *realecho=(double*)mxGetPr(prhs[2]);
    double *imagecho=(double*)mxGetPr(prhs[3]);
     /*for(int i=0;i<M*L;++i){
         mexPrintf("%.8f,",realecho[i]);
         if (i==20) mexPrintf("\n");
     }*/
    plhs[0]=mxCreateNumericMatrix(M,L,mxDOUBLE_CLASS,mxREAL);//开辟输出参数空间
    plhs[1]=mxCreateNumericMatrix(M,L,mxDOUBLE_CLASS,mxREAL);//开辟输出参数空间
    double *realpc=(double*)mxGetPr(plhs[0]);//获得输出结果
    double *imagpc=(double*)mxGetPr(plhs[1]);
    pc(realh,imagh,realecho,imagecho,realpc,imagpc,M,L);
}