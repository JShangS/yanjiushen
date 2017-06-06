#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include "PC.h"
//#pragma comment(lib,"cufft.lib")
void mexFunction(int nlhs,mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    if (nrhs!=4)
    {
      //mexErrMsgTxt("�������벻��");
      mexPrintf(" ��������Ϊ:%d \n Ӧ��Ϊ:4,\n Ҫ�д��ݺ���ʵ�����鲿���ز�ʵ�����鲿\n",nrhs);   
    }
    mexPrintf(" �������Ϊ:%d\n",nlhs);
    if (nlhs!=2)   
        mexErrMsgTxt("�����������,Ҫ��2������ѹʱ���ʵ�����鲿");
    int L=mxGetM(prhs[2]);//�������
    int M=mxGetN(prhs[2]);//�������
    mexPrintf("L=%d,M=%d\n",L,M); 
    double *realh=(double*)mxGetPr(prhs[0]);//��ȡ�������
   // mexPrintf("��ַ:realh[%d]=%d,realh[%d]=%d,realh[%d]=%d\n",0,&realh[0],1,&realh[1],2,&realh[2]); 
    double *imagh=(double*)mxGetPr(prhs[1]);
    double *realecho=(double*)mxGetPr(prhs[2]);
    double *imagecho=(double*)mxGetPr(prhs[3]);
     /*for(int i=0;i<M*L;++i){
         mexPrintf("%.8f,",realecho[i]);
         if (i==20) mexPrintf("\n");
     }*/
    plhs[0]=mxCreateNumericMatrix(M,L,mxDOUBLE_CLASS,mxREAL);//������������ռ�
    plhs[1]=mxCreateNumericMatrix(M,L,mxDOUBLE_CLASS,mxREAL);//������������ռ�
    double *realpc=(double*)mxGetPr(plhs[0]);//���������
    double *imagpc=(double*)mxGetPr(plhs[1]);
    pc(realh,imagh,realecho,imagecho,realpc,imagpc,M,L);
}