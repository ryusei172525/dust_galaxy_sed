#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//to read the file of Karakas (2010)
typedef struct Data{
  double mass[15][77];
  char set[15][77][2];
  double metallicity[15][77];
  double M_final[15][77];
  char element[15][77][2];
  int atomic[15][77];
  double yield[15][77];
  double lostwind[15][77];
  double initwind[15][77];
  double averabun[15][77];
  double initfra[15][77];
  double pridfac[15][77];
}Metal;

void readfile(Metal *,char *);
void linear(Metal *);

int main()
{
  int i;
  Metal agbdata[4];
  char filename[256]; 

  readfile(&agbdata[0],"AGBmetaldata_z0.0001.dat");
  readfile(&agbdata[1],"AGBmetaldata_z0.004.dat");
  readfile(&agbdata[2],"AGBmetaldata_z0.008.dat");
  readfile(&agbdata[3],"AGBmetaldata_z0.02.dat");
  //printf("%5e\n",agbdata[3].mass[14][0]);
  
  linear(agbdata);
   
}

void readfile(Metal *p,char *filename)
{
  int i,j;
  FILE *file;

 if ((file = fopen(filename,"r")) == NULL){
    printf("file open is impossible\n");
    exit(1);
  }
 for(j=0;j<15;j++){
   for(i=0;i<77;i++){
     fscanf(file,"%lf",&p->mass[j][i]);
     fscanf(file,"%s",p->set[j][i]);
     fscanf(file,"%lf",&p->metallicity[j][i]);
     fscanf(file,"%lf",&p->M_final[j][i]);
     fscanf(file,"%s",p->element[j][i]);
     fscanf(file,"%d",&p->atomic[j][i]);
     fscanf(file,"%lf",&p->yield[j][i]);
     fscanf(file,"%lf",&p->lostwind[j][i]);
     fscanf(file,"%lf",&p->initwind[j][i]);
     fscanf(file,"%lf",&p->averabun[j][i]);
     fscanf(file,"%lf",&p->initfra[j][i]);
     fscanf(file,"%lf",&p->pridfac[j][i]);
     
     if(feof(file)){
       break;  
     }//ファイルの末尾まで読み込み
   }
 }
 fclose(file);
}

void linear(Metal *agbdata)
{
  int i,j,k,l;
  double metal1,metal2,logm1,logm2,metal,logm;
  double zyield1,zyield2,logz1,logz2,logz,z;
  double cyield1,cyield2,logc1,logc2,logc,c;
  double syield1,syield2,logs1,logs2,logs,s;
  FILE *file;
  
 
    for(l=0;l<15;l++){//stellar mass
      for(i=0;i<3;i++){//metallicity

	zyield1 = zyield2 = cyield1 = cyield2 = syield1 = syield2 = 0.0;

	for(j=6;j<77;j++){//Heよりも重い元素の足し算
	  zyield1 += (agbdata+i)->lostwind[l][j];
	  zyield2 += (agbdata+i+1)->lostwind[l][j];
	}
	
	for(j=9;j<12;j++){//carbonの足し算
	  cyield1 += (agbdata+i)->lostwind[l][j];
	  cyield2 += (agbdata+i+1)->lostwind[l][j];
	}

	for(j=44;j<51;j++){//siliconの足し算
	  syield1 += (agbdata+i)->lostwind[l][j];
	  syield2 += (agbdata+i+1)->lostwind[l][j];
	}
	
	metal1=(agbdata+i)->metallicity[0][0]/0.02;//solar metallicityに規格化
	metal2=(agbdata+i+1)->metallicity[0][0]/0.02;      
	
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
	  strcpy(prename,"newyield/AGBmetal_interpolate/metal_interpolate/AGB_z");
	  strcpy(postname,".dat");
	  
	  strcat(prename,label);
	  strcat(prename,postname);
	  
	  if ((file = fopen(prename,"a+")) == NULL){
	    printf("file open is impossible\n");
	    exit(1);
	  }
	  
	  fprintf(file,"%5e    %5e    %5e    %5e\n",(agbdata)->mass[l][0],z,c,s);
	  fclose(file);
	  
	  metal += 0.001;//0.001Zsun刻み
	  
	}
      }
    }

    //for plot at metallicity 1.0Zsun
     for(l=0;l<15;l++){//stellar mass
      
	zyield1 = cyield1 = syield1 = 0.0;

	for(j=6;j<77;j++){//Heよりも重い元素の足し算
	  zyield1 += (agbdata+3)->lostwind[l][j];
	}
	
	for(j=9;j<12;j++){//carbonの足し算
	  cyield1 += (agbdata+3)->lostwind[l][j];
	}

	for(j=44;j<51;j++){//siliconの足し算
	  syield1 += (agbdata+3)->lostwind[l][j];
	}

	metal = 1.000;
	
	  char label[80];
	  sprintf(label,"%3.3lf",metal);
	  //printf("%s\n",c);
	  
	  char prename[256],postname[80];
	  strcpy(prename,"newyield/AGBmetal_interpolate/metal_interpolate/AGB_z");
	  strcpy(postname,".dat");
	  
	  strcat(prename,label);
	  strcat(prename,postname);
	  
	  if ((file = fopen(prename,"a+")) == NULL){
	    printf("file open is impossible\n");
	    exit(1);
	  }
	  
	  fprintf(file,"%5e    %5e    %5e    %5e\n",(agbdata)->mass[l][0],zyield1,cyield1,syield1);
	  fclose(file);
	  
	}
      

}

