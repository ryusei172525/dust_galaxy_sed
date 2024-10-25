/*calculation of dust size distribution for AGB stars (lognormal)*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int main()
{
  int i;
  double aa0,aa1;
  double mu,sigma;
  double mode,ave,var;

  /*parameter set: Winters et al. 1997; grain size distribution*/
  /*mu = log(109.49);
  sigma = log10(2.0);
  /********************************/


  /*parameter set: Yasuda & Kozasa 2012; mass spectrum*/
  //mu = log(309.49);
  //sigma = log10(4.0);
    /********************************/

  /*parameter set: Yasuda & Kozasa 2012; mass spectrum*/
  //resize on 24/2/2012, file name is agb_dsdis_yasudakozasa_resize.dat
  //mu = log(559.49*1.e-7);
  //sigma = log10(4.0);
  /********************************/

 /*parameter set: Yasuda & Kozasa 2012; mass spectrum*/
  //resize on 24/2/2012, file name is agb_dsdis_yasudakozasa_resize_carbon.dat
  //mu = log(109.49*1.e-7);
  //sigma = log10(3.0);
  /********************************/

 /*parameter set: Yasuda & Kozasa 2012; mass spectrum*/
  //resize on 24/2/2012, file name is agb_dsdis_ventura.dat
  //extreme case
  //mu = log(65.49*1.e-7);
  //sigma = log10(5.0);
  /********************************/

  /*parameter set: Yasuda & Kozasa 2012; mass spectrum, peak value changes*/
  //resize on 24/2/2012, file name is agb_dsdis_yasudakozasa_resize_carbon.dat
  //mu = log(309.49*1.e-7);
  //sigma = log10(3.0);
  /********************************/

  /*parameter set: */
  //test size distribution for study about initial stellar GSD, 28/8/2013, file name is stellar_testGSD_**um.dat, ** is the peak size.
  mu = log(1009.49*1.e-7);
  sigma = log10(3.0);
  /********************************/



  mode=exp(mu)/exp(pow(sigma,2.0));
  ave=exp(mu)*exp(pow(sigma,2.0)/2.0);
  var=(exp(pow(sigma,2.0))-1.0)*pow(exp(mu),2.0)*exp(pow(sigma,2.0));

  /*printf("%5e    %5e    %5e\n",mode,ave,var);*/

  //******以下サイズ分布関数を作成(Winters et al. 1997)*******
  /*for(i=1;i<=19;i++){
    printf("%5e     %5e\n",pow(10.0,-8.0+0.1*i),0.0);
      }
  for(i = 0;i <= 20;i++){
    //aa0 = 10.0 + 990.0*(1.0*i/100.0);//サイズの単位はnm
    
    aa0 = pow(10.0,1.0+0.1*i);

    aa1 = 1.0/(aa0*sigma*sqrt(2.0*M_PI))*exp(-pow(log(aa0)-mu,2.0)/2.0/sigma/sigma);
   
    printf("%5e     %5e\n",aa0*pow(10.0,-7.0),aa1);

 }  
  for(i=1;i<=9;i++){
    printf("%5e     %5e\n",pow(10.0,-4.0+0.1*i),0.0);
    }*/

//******以下サイズ分布関数を作成(Yasuda & Kozasa 2011)*******
  /* for(i=1;i<=19;i++){
    printf("%5e     %5e\n",pow(10.0,-8.0+0.1*i),0.0);
      }
  for(i = 0;i <= 26;i++){
    //aa0 = 10.0 + 990.0*(1.0*i/100.0);//サイズの単位はnm
    
    //aa0 = pow(10.0,1.0+0.1*i);
    aa0 = pow(10.0,-6.0+0.1*i);
    aa1 = 1.0/(aa0*sigma*sqrt(2.0*M_PI))*exp(-pow(log(aa0)-mu,2.0)/2.0/sigma/sigma)/pow(aa0,4.0);
    //aa1 = 1.0/(aa0*sigma*sqrt(2.0*M_PI))*exp(-pow(log(aa0)-mu,2.0)/2.0/sigma/sigma);
    //printf("%5e     %5e\n",aa0*pow(10.0,-7.0),aa1);
    printf("%5e     %5e\n",aa0,aa1);
 }  
  for(i=1;i<=3;i++){
    printf("%5e     %5e\n",pow(10.0,-3.4+0.1*i),0.0);
    }*/

  for(i = 1;i <= 49;i++){
    
    aa0 = pow(10.0,-8.0+0.1*i);
    aa1 = 1.0/(aa0*sigma*sqrt(2.0*M_PI))*exp(-pow(log(aa0)-mu,2.0)/2.0/sigma/sigma)/pow(aa0,4.0);
    //aa1 = 1.0/(aa0*sigma*sqrt(2.0*M_PI))*exp(-pow(log(aa0)-mu,2.0)/2.0/sigma/sigma);
    printf("%5e     %5e\n",aa0,aa1);
    }  



}

