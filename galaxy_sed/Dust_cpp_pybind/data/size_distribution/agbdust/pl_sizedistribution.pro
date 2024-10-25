g1=read_file('agb_dsdis_new.dat',/dbl,/str)
;最後の,/dbl,/strは読み込む値をdoubleとして読み込んでいる(これがない場合はfloatとして読み込んでいる)

red=[0,1,1.0,0,0,0,1,1,0.8,0.2,0.5,1.0,0.5]
grn=[0,0,0.5,1,0,1,1,0,0.8,0.5,0.0,0.5,0.5]
blu=[0,0,0.0,0,1,1,1,1,0.2,0.5,0.9,0.2,0.5]

tvlct,red*255,grn*255,blu*255

set_plot,'ps'
;place='/home/asano/dust_evolution/dust_destruction/'
place='/home/asano/dust_evolution/dust_destruction/size_distribution/agbdust/'

device,bits=8,filename=place+'gradis_agb.ps'
device,xsize=16,ysize=15,xoffset=3,yoffset=8
device,/color

xmin=1.e-8
xmax=1.e-3
ymin=1.e-14
ymax=1.e

  plot,g1(0,*),g1(1,*),title='!6',/xstyle,/xlog,/ystyle,/ylog,$
   xtitle=textoidl('!6grain size [cm]'),$
    ytitle=textoidl('!6grain size distribution(relative)'),$
    ;ytitle=textoidl('!6dust mass[M_{'+sunsymbol()+'}]'),$
       /nodata,xthick=2.0, ythick=2.0,chars=1.8,charthick=2,$
    xrange=[xmin,xmax],$
    yrange=[ymin,ymax]

  oplot,g1(0,*),g1(1,*),linestyle=0,symsize=1.5,thick=4,color=1;for 100Myr
  ;oplot,g1(0,0:48),g1(11,4900:4948),linestyle=1,symsize=1.5,thick=4,color=4;for 1000 Myr
  ;oplot,g1(0,0:48),g1(11,49000:49048),linestyle=2,symsize=1.5,thick=4,color=7;for 10000 Myr
  

;  oplot,g3(0,*)*0.1,g3(3,*),linestyle=1,symsize=1.5,thick=4,color=4
;  oplot,g4(0,*),g4(3,*),linestyle=1,symsize=1.5,thick=4,color=7
;  oplot,g3(0,*)*0.1,g3(3,*),linestyle=1,symsize=1.5,thick=4,color=4
;g2は時間スケールが1Myr刻みであるが、それではplotされない部分もあるの
;で、g3でその部分をカバーしている(0.1Myr刻み)。  


;plots,[0.0002*xmax,0.0004*xmax],[1.0e+6*ymin,1.0e+6*ymin],$
      ;linestyle=0,symsize=1.5,thick=4,color=1
  ;xyouts,0.03*xmax,1.e+4*ymin,'!6AGB',charthick=2,chars=1.5,color=1

  ;plots,[0.0002*xmax,0.0004*xmax],[1.7e+5*ymin,1.7e+5*ymin],$
      ;linestyle=0,symsize=1.5,thick=4,color=4
  ;xyouts,0.002*xmax,3.e+2*ymin,'!6SN',charthick=2,chars=1.5,color=4
 
device,/close
set_plot,'x'


end

