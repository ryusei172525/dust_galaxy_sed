#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

void readfile(int,double agbdata[][4],double sndata[][4]);
void linear(int,double agbdata[][4],double sndata[][4]);

int main()
{
  
  int j;
  double agbdata[15][4],sndata[7][4];

  //loop for metallicity
  for(j=5;j<=1000;j++){
      
    readfile(j,agbdata,sndata);
    linear(j,agbdata,sndata);

  }
  
}

void readfile(int j,double agbdata[][4],double sndata[][4])
{
  int i,k;
  FILE *file1,*file2;
  double z;
  char label[80],prename[256],postname[80];
 
  z = 0.001*j;
  
  sprintf(label,"%3.3lf",z);
  strcpy(prename,"newyield/AGBmetal_interpolate/metal_interpolate/AGB_z");
  strcpy(postname,".dat");
  
  strcat(prename,label);
  strcat(prename,postname);
  
 if ((file1 = fopen(prename,"r")) == NULL){
    printf("file open is impossible\n");
    exit(1);
  }
 for(k=0;k<15;k++){
   for(i=0;i<4;i++){
     fscanf(file1,"%lf",&agbdata[k][i]);
     if(feof(file1)){
       break;  
     }//ファイルの末尾まで読み込み
   }
 }

 strcpy(prename,"newyield/SNmetal_kobayashi_interpolate/metal_interpolate/SN_z");
 strcpy(postname,".dat");
 
 strcat(prename,label);
 strcat(prename,postname);
 
 if ((file2 = fopen(prename,"r")) == NULL){
    printf("file open is impossible\n");
    exit(1);
  }
 for(k=0;k<7;k++){
   for(i=0;i<4;i++){
     fscanf(file2,"%lf",&sndata[k][i]);
     if(feof(file2)){
       break;  
     }//ファイルの末尾まで読み込み
   }
 }

 fclose(file1);
 fclose(file2);

}


void linear(int j,double agbdata[][4],double sndata[][4])
{
  int i,k,l;
  double smass,smass1,smass2,logm1,logm2,logm;
  double zyield1,zyield2,logz1,logz2,logz,z;
  double cyield1,cyield2,logc1,logc2,logc,c;
  double syield1,syield2,logs1,logs2,logs,s;
  FILE *file;
  double zz;
  char label[80],prename[256],postname[80];
  
  zz = 0.001*j;
  
  sprintf(label,"%3.3lf",zz);
  strcpy(prename,"newyield/totalmetal/totalmetalmass_z");
  strcpy(postname,".dat");
  
  strcat(prename,label);
  strcat(prename,postname);
  
  if ((file = fopen(prename,"w")) == NULL){
    printf("file open is impossible\n");
    exit(1);
  }
  
  
  //AGB part     	
  
  for(i=0;i<14;i++){//AGB stellar mass - 0.01Msunまで
    
    smass1 = agbdata[i][0];
    smass2 = agbdata[i+1][0];	 
    zyield1 = agbdata[i][1];
    zyield2 = agbdata[i+1][1];
    cyield1 = agbdata[i][2];
    cyield2 = agbdata[i+1][2];
    syield1 = agbdata[i][3];
    syield2 = agbdata[i+1][3];	  
    
    smass = agbdata[i][0];      
    while(agbdata[i+1][0] - smass > 0.005){//doubleはintと異なり少数に値がある可能性がある。それを回避するため差をとった。
      
      //more heavier than He
      logz1=log10(zyield1);
      logz2=log10(zyield2);
      logm1=log10(smass1);
      logm2=log10(smass2);
      logm=log10(smass);
      
      //Carbon
      logc1=log10(cyield1);
      logc2=log10(cyield2);
      
      //Silicon
      logs1=log10(syield1);
      logs2=log10(syield2);
      
      logz=(logz2-logz1)/(logm2-logm1)*logm+(logm2*logz1-logm1*logz2)/(logm2-logm1);
      logc=(logc2-logc1)/(logm2-logm1)*logm+(logm2*logc1-logm1*logc2)/(logm2-logm1);
      logs=(logs2-logs1)/(logm2-logm1)*logm+(logm2*logs1-logm1*logs2)/(logm2-logm1);
      
      z=pow(10.0,logz);
      c=pow(10.0,logc);
      s=pow(10.0,logs);
      
      fprintf(file,"%5e    %5e    %5e    %5e\n",smass,z,c,s);
            
      smass += 0.01;//0.01Msun刻み
      
    }
  }      
  
	//AGB and SN part

  smass1 = agbdata[14][0];
  smass2 = sndata[0][0];	 
  zyield1 = agbdata[14][1];
  zyield2 = sndata[0][1];
  cyield1 = agbdata[14][2];
  cyield2 = sndata[0][2];
  syield1 = agbdata[14][3];
  syield2 = sndata[0][3];	  
  
  smass = agbdata[14][0];      
  while(sndata[0][0] - smass > 0.005){    
    
    //more heavier than He
    logz1=log10(zyield1);
    logz2=log10(zyield2);
    logm1=log10(smass1);
    logm2=log10(smass2);
    logm=log10(smass);
      
    //Carbon
    logc1=log10(cyield1);
    logc2=log10(cyield2);
      
    //Silicon
    logs1=log10(syield1);
    logs2=log10(syield2);
      
    logz=(logz2-logz1)/(logm2-logm1)*logm+(logm2*logz1-logm1*logz2)/(logm2-logm1);
    logc=(logc2-logc1)/(logm2-logm1)*logm+(logm2*logc1-logm1*logc2)/(logm2-logm1);
    logs=(logs2-logs1)/(logm2-logm1)*logm+(logm2*logs1-logm1*logs2)/(logm2-logm1);
    
    z=pow(10.0,logz);
    c=pow(10.0,logc);
    s=pow(10.0,logs);
    
    fprintf(file,"%5e    %5e    %5e    %5e\n",smass,z,c,s);
    
    smass += 0.01;//0.01Msun刻み
    
  }
      
  //SN part    
  
  for(i=0;i<6;i++){//SN stellar mass - 0.01Msunまで
    
    smass1 = sndata[i][0];
    smass2 = sndata[i+1][0];	 
    zyield1 = sndata[i][1];
    zyield2 = sndata[i+1][1];
    cyield1 = sndata[i][2];
    cyield2 = sndata[i+1][2];
    syield1 = sndata[i][3];
    syield2 = sndata[i+1][3];	  
    
    smass = sndata[i][0];      
    while(sndata[i+1][0] - smass > 0.005){    
      
      //more heavier than He
      logz1=log10(zyield1);
      logz2=log10(zyield2);
      logm1=log10(smass1);
      logm2=log10(smass2);
      logm=log10(smass);
      
      //Carbon
      logc1=log10(cyield1);
      logc2=log10(cyield2);
      
      //Silicon
      logs1=log10(syield1);
      logs2=log10(syield2);
      
      logz=(logz2-logz1)/(logm2-logm1)*logm+(logm2*logz1-logm1*logz2)/(logm2-logm1);
      logc=(logc2-logc1)/(logm2-logm1)*logm+(logm2*logc1-logm1*logc2)/(logm2-logm1);
      logs=(logs2-logs1)/(logm2-logm1)*logm+(logm2*logs1-logm1*logs2)/(logm2-logm1);
      
      z=pow(10.0,logz);
      c=pow(10.0,logc);
      s=pow(10.0,logs);
      
      fprintf(file,"%5e    %5e    %5e    %5e\n",smass,z,c,s);
            
      smass += 0.01;//0.01Msun刻み
      
    }
  }      

  //for plot at stellar mass 40.0Msun
  
  fprintf(file,"%5e    %5e    %5e    %5e\n",sndata[6][0],sndata[6][1],sndata[6][2],sndata[6][3]);
  
  fclose(file);
  
}

