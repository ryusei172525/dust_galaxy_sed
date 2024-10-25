#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//to read the file of Kobayashi et al. (2006)
typedef struct Data{
  double metallicity[85];
  char set[85][10];
  double yield[7][85];
}Metal;

void readfile(Metal *,char *);
void linear(Metal *);

int main()
{
  int i;
  Metal sndata[4];
  char filename[256]; 

  readfile(&sndata[0],"SNmetal_kobayashi_z0.0.dat");
  readfile(&sndata[1],"SNmetal_kobayashi_z0.001.dat");
  readfile(&sndata[2],"SNmetal_kobayashi_z0.004.dat");
  readfile(&sndata[3],"SNmetal_kobayashi_z0.02.dat");
  //printf("%5e\n",sndata[3].mass[14][0]);
  
  linear(sndata);
   
}

void readfile(Metal *p,char *filename)
{
  int i,j;
  FILE *file;

 if ((file = fopen(filename,"r")) == NULL){
    printf("file open is impossible\n");
    exit(1);
  }
 for(j=0;j<85;j++){
   fscanf(file,"%lf",&p->metallicity[j]);
   fscanf(file,"%s",p->set[j]);
   for(i=0;i<7;i++){
     fscanf(file,"%lf",&p->yield[i][j]);
     if(feof(file)){
       break;  
     }//ファイルの末尾まで読み込み
   }
 }
 fclose(file);
}

void linear(Metal *sndata)
{
  int i,j,k,l;
  double metal1,metal2,logm1,logm2,metal,logm;
  double zyield1,zyield2,logz1,logz2,logz,z;
  double cyield1,cyield2,logc1,logc2,logc,c;
  double syield1,syield2,logs1,logs2,logs,s;
  FILE *file;
  
 
    for(l=0;l<7;l++){//stellar mass
      for(i=0;i<3;i++){//metallicity

	zyield1 = zyield2 = cyield1 = cyield2 = syield1 = syield2 = 0.0;

	for(j=6;j<85;j++){//Heよりも重い元素の足し算
	  zyield1 += (sndata+i)->yield[l][j];
	  zyield2 += (sndata+i+1)->yield[l][j];
	}
	
	for(j=11;j<13;j++){//carbonの足し算
	  cyield1 += (sndata+i)->yield[l][j];
	  cyield2 += (sndata+i+1)->yield[l][j];
	}
	
	for(j=27;j<30;j++){//siliconの足し算
	  syield1 += (sndata+i)->yield[l][j];
	  syield2 += (sndata+i+1)->yield[l][j];
	}
	
	metal1=(sndata+i)->metallicity[0]/0.02;//solar metallicityに規格化
	metal2=(sndata+i+1)->metallicity[0]/0.02;      
	
	if(metal1 == 0.0) metal1 = 1.e-8;

	metal = metal1;      
	while(metal < metal2){    
	  
	  //more heavier than He
	  logz1=log10(zyield1);
	  logz2=log10(zyield2);
	  logm1=log10(metal1);
	  logm2=log10(metal2);
	  logm=log10(metal);
	  
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

	  char label[80];
	  sprintf(label,"%3.3lf",metal);
	  //printf("%s\n",c);
	  
	  char prename[256],postname[80];
	  strcpy(prename,"newyield/SNmetal_kobayashi_interpolate/metal_interpolate/SN_z");
	  strcpy(postname,".dat");
	  
	  strcat(prename,label);
	  strcat(prename,postname);
	  
	  if ((file = fopen(prename,"a+")) == NULL){
	    printf("file open is impossible\n");
	    exit(1);
	  }
	  
	  fprintf(file,"%5e    %5e    %5e    %5e\n",(sndata)->yield[l][0],z,c,s);
	  fclose(file);
	  
	  metal += 0.001;//0.001Zsun刻み
	  
	}
      }
    }

    //for plot at metallicity 1.0Zsun
    for(l=0;l<7;l++){//stellar mass
      
      zyield1 = cyield1 = syield1 = 0.0;
      
      for(j=6;j<85;j++){//Heよりも重い元素の足し算
	zyield1 += (sndata+3)->yield[l][j];
      }
      
      for(j=11;j<13;j++){//carbonの足し算
	cyield1 += (sndata+3)->yield[l][j];
      }
      
      for(j=27;j<30;j++){//siliconの足し算
	syield1 += (sndata+3)->yield[l][j];
      }
      
      metal = 1.000;
      
      char label[80];
      sprintf(label,"%3.3lf",metal);
      //printf("%s\n",c);
	  
      char prename[256],postname[80];
      strcpy(prename,"newyield/SNmetal_kobayashi_interpolate/metal_interpolate/SN_z");
      strcpy(postname,".dat");
      
      strcat(prename,label);
      strcat(prename,postname);
      
      if ((file = fopen(prename,"a+")) == NULL){
	printf("file open is impossible\n");
	exit(1);
      }
      
      fprintf(file,"%5e    %5e    %5e    %5e\n",(sndata)->yield[l][0],zyield1,cyield1,syield1);
      fclose(file);
      
    }
    
}

