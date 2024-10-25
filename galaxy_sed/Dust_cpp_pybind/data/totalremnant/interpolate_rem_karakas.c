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
  FILE *file;
  
 
    for(l=0;l<15;l++){//stellar mass
      for(i=0;i<3;i++){//metallicity

	metal1=(agbdata+i)->metallicity[0][0]/0.02;//solar metallicityに規格化
	metal2=(agbdata+i+1)->metallicity[0][0]/0.02;      

	zyield1 = (agbdata+i)->M_final[l][0];
	zyield2 = (agbdata+i+1)->M_final[l][0];
	
	metal = metal1;      
	while(metal < metal2){    
	  
	  logz1=log10(zyield1);
	  logz2=log10(zyield2);
	  logm1=log10(metal1);
	  logm2=log10(metal2);
	  logm=log10(metal);
	  
	  logz=(logz2-logz1)/(logm2-logm1)*logm+(logm2*logz1-logm1*logz2)/(logm2-logm1);
	
	  z=pow(10.0,logz);
	
	  char label[80];
	  sprintf(label,"%3.3lf",metal);
	  
	  char prename[256],postname[80];
	  strcpy(prename,"newyield/AGBrem_interpolate/AGB_z");
	  strcpy(postname,".dat");
	  
	  strcat(prename,label);
	  strcat(prename,postname);
	  
	  if ((file = fopen(prename,"a+")) == NULL){
	    printf("file open is impossible\n");
	    exit(1);
	  }
	  
	  fprintf(file,"%5e    %5e\n",(agbdata)->mass[l][0],z);
	  fclose(file);
	  
	  metal += 0.001;//0.001Zsun刻み
	  
	}
      }
    }

    //for plot at metallicity 1.0Zsun
     for(l=0;l<15;l++){//stellar mass
      
       metal = 1.000;
       zyield1 = (agbdata+3)->M_final[l][0];

	  char label[80];
	  sprintf(label,"%3.3lf",metal);
	  //printf("%s\n",c);
	  
	  char prename[256],postname[80];
	  strcpy(prename,"newyield/AGBrem_interpolate/AGB_z");
	  strcpy(postname,".dat");
	  
	  strcat(prename,label);
	  strcat(prename,postname);
	  
	  if ((file = fopen(prename,"a+")) == NULL){
	    printf("file open is impossible\n");
	    exit(1);
	  }
	  
	  fprintf(file,"%5e    %5e\n",(agbdata)->mass[l][0],zyield1);
	  fclose(file);
	  
	}
      

}

