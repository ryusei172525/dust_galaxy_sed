g1=read_file('30n1_linear.dat',/dbl,/str)
;最後の,/dbl,/strは読み込む値をdoubleとして読み込んでいる(これがない場合はfloatとして読み込んでいる)

red=[0,1,1.0,0,0,0,1,1,0.8,0.2,0.5,1.0,0.5]
grn=[0,0,0.5,1,0,1,1,0,0.8,0.5,0.0,0.5,0.5]
blu=[0,0,0.0,0,1,1,1,1,0.2,0.5,0.9,0.2,0.5]

tvlct,red*255,grn*255,blu*255

set_plot,'ps'
;place='/home/asano/dust_evolution/dust_destruction/'
place='/home/ryosuke/dust_evolution/dust_destruction/size_distribution/shockeddust/linear/num_den1/'

device,bits=8,filename=place+'massdis_sn30.ps'
device,xsize=16,ysize=15,xoffset=3,yoffset=8
device,/color

xmin=1.e-8
xmax=1.e-3
;ymin=1.e+10
;ymax=1.e+25
ymin=1.e-10
ymax=1.e+1

  plot,g1(0,*),g1(1,*),title='!6',/xstyle,/xlog,/ystyle,/ylog,$
   xtitle=textoidl('!6grain size [cm]'),$
    ;ytitle=textoidl('!6grain size distribution(relative)'),$
    ytitle=textoidl('!6a^{4}f(a)'),$
       /nodata,xthick=2.0, ythick=2.0,chars=1.8,charthick=2,$
    xrange=[xmin,xmax],$
    yrange=[ymin,ymax]

  oplot,g1(0,*),g1(1,*)*g1(0,*)*g1(0,*)*g1(0,*)*g1(0,*),linestyle=0,symsize=1.5,thick=4,color=1
  oplot,g1(0,*),g1(2,*)*g1(0,*)*g1(0,*)*g1(0,*)*g1(0,*),linestyle=1,symsize=1.5,thick=4,color=2
  oplot,g1(0,*),g1(3,*)*g1(0,*)*g1(0,*)*g1(0,*)*g1(0,*),linestyle=2,symsize=1.5,thick=4,color=3
  oplot,g1(0,*),g1(4,*)*g1(0,*)*g1(0,*)*g1(0,*)*g1(0,*),linestyle=3,symsize=1.5,thick=4,color=4
  oplot,g1(0,*),g1(5,*)*g1(0,*)*g1(0,*)*g1(0,*)*g1(0,*),linestyle=4,symsize=1.5,thick=4,color=5
  oplot,g1(0,*),g1(6,*)*g1(0,*)*g1(0,*)*g1(0,*)*g1(0,*),linestyle=5,symsize=1.5,thick=4,color=7
  oplot,g1(0,*),g1(7,*)*g1(0,*)*g1(0,*)*g1(0,*)*g1(0,*),linestyle=6,symsize=1.5,thick=4,color=8
  oplot,g1(0,*),g1(8,*)*g1(0,*)*g1(0,*)*g1(0,*)*g1(0,*),linestyle=7,symsize=1.5,thick=4,color=9
  oplot,g1(0,*),g1(9,*)*g1(0,*)*g1(0,*)*g1(0,*)*g1(0,*),linestyle=8,symsize=1.5,thick=4,color=10

                                ;oplot,g1(0,0:48),g1(11,4900:4948),linestyle=1,symsize=1.5,thick=4,color=4;for 1000 Myr
                                ;oplot,g1(0,0:48),g1(11,49000:49048),linestyle=2,symsize=1.5,thick=4,color=7;for 10000 Myr


;  oplot,g3(0,*)*0.1,g3(3,*),linestyle=1,symsize=1.5,thick=4,color=4
;  oplot,g4(0,*),g4(3,*),linestyle=1,symsize=1.5,thick=4,color=7
;  oplot,g3(0,*)*0.1,g3(3,*),linestyle=1,symsize=1.5,thick=4,color=4
;g2は時間スケールが1Myr刻みであるが、それではplotされない部分もあるの
;で、g3でその部分をカバーしている(0.1Myr刻み)。


  plots,[0.00002*xmax,0.00007*xmax],[1.0e+10*ymin,1.0e+10*ymin],$
        linestyle=0,symsize=1.5,thick=4,color=1
  xyouts,0.00008*xmax,7.e+9*ymin,'!6C',charthick=2,chars=1.5,color=1

  plots,[0.00002*xmax,0.00007*xmax],[1.0e+9*ymin,1.0e+9*ymin],$
        linestyle=1,symsize=1.5,thick=4,color=2
  xyouts,0.00008*xmax,7.e+8*ymin,'!6Si',charthick=2,chars=1.5,color=2

  plots,[0.00002*xmax,0.00007*xmax],[1.0e+8*ymin,1.0e+8*ymin],$
        linestyle=2,symsize=1.5,thick=4,color=3
  xyouts,0.00008*xmax,7.e+7*ymin,'!6Fe',charthick=2,chars=1.5,color=3

  plots,[0.0002*xmax,0.0007*xmax],[1.0e+10*ymin,1.0e+10*ymin],$
        linestyle=3,symsize=1.5,thick=4,color=4
  xyouts,0.0008*xmax,7.e+9*ymin,'!6FeS',charthick=2,chars=1.5,color=4

  plots,[0.0002*xmax,0.0007*xmax],[1.0e+9*ymin,1.0e+9*ymin],$
        linestyle=4,symsize=1.5,thick=4,color=5
  xyouts,0.0008*xmax,7.e+8*ymin,textoidl('!6Al_{2}O_{3}'),charthick=2,chars=1.5,color=5

  plots,[0.0002*xmax,0.0007*xmax],[1.0e+8*ymin,1.0e+8*ymin],$
        linestyle=5,symsize=1.5,thick=4,color=7
  xyouts,0.0008*xmax,7.e+7*ymin,textoidl('!6MgSiO_{3}'),charthick=2,chars=1.5,color=7

  plots,[0.008*xmax,0.04*xmax],[1.0e+10*ymin,1.0e+10*ymin],$
        linestyle=6,symsize=1.5,thick=4,color=8
  xyouts,0.05*xmax,7.e+9*ymin,textoidl('!6Mg_{2}SiO_{4}'),charthick=2,chars=1.5,color=8

  plots,[0.008*xmax,0.04*xmax],[1.0e+9*ymin,1.0e+9*ymin],$
        linestyle=7,symsize=1.5,thick=4,color=9
  xyouts,0.05*xmax,7.e+8*ymin,textoidl('!6SiO_{2}'),charthick=2,chars=1.5,color=9

  plots,[0.008*xmax,0.04*xmax],[1.0e+8*ymin,1.0e+8*ymin],$
        linestyle=8,symsize=1.5,thick=4,color=10
  xyouts,0.05*xmax,7.e+7*ymin,textoidl('!6MgO'),charthick=2,chars=1.5,color=10

                                ;RHO_CARSIL densities of dust species
                                ;   1       C  2.2631e+00    u
                                ;   2      Si  2.3244e+00    u
                                ;   3      Fe  7.8913e+00    u
                                ;   4     FeS  4.8358e+00    u
                                ;   5   Al2O3  3.9824e+00   m/u
                                ;   6  MgSiO3  3.1813e+00   m/u
                                ;   7 Mg2SiO4  3.2018e+00   m/u
                                ;   8    SiO2  2.6428e+00   m/u
                                ;   9     MgO  3.5633e+00    u
                                ;  10   Fe3O4  5.2073e+00    m
                                ; *******************************/

  device,/close
  set_plot,'x'


end

