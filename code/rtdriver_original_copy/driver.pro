pro driver
close, /all
!p.multi=[0,2,4]
   set_plot, 'PS'
   device, /inch, /color, /portrait, ysize = 9.5, xsize = 7.75, $
   yoffset = 0.75, xoffset = 0.5, file='driver.ps'
   loadct, 39
;
;  ymax/min for plots
;
ypmax = 1.0e-6
ypmin = 10.0
;
xp  = fltarr(2)
yp  = fltarr(2)
num = 0
openr, 1, 'driver_data'
readf,1,num
;
albedo  = 0.0
cosz    = 0.0
solar   = 0.0
tautot  = 0.0
tautotc = 0.0
conrnu  = 0.0
firmax  = 0.0
firmin  = 0.0
fvmax   = 0.0
fvmin   = 0.0
fdmax   = 0.0
fdmin   = 0.0
firup   = fltarr(num)
firdn   = fltarr(num)
fvup    = fltarr(num)
fvdn    = fltarr(num)
fir     = fltarr(num)
fv      = fltarr(num)
pres    = fltarr(num)
heatv   = fltarr(num)
heatir  = fltarr(num)
heatt   = fltarr(num)
temp    = fltarr(num)
pres2   = fltarr(num-1)
;
; Upward fluxes
;
readf,1, firmax, firmin, fvmax, fvmin
for n=0,num-1 do begin
  readf,1, p1, f1, f2
  pres(n)  = p1
  fvup(n)  = f1
  firup(n) = f2
endfor
;
if (fvmax le 0.1) then begin
  fvmin = -0.5
  fvmax = 0.5
endif
plot, fvup, pres, xrange=[fvmin,fvmax], xs=1, $
/ylog, yrange = [ypmin,ypmax], $
symsize=0.6, psym = -4, $
title = 'UPWARD VISIBLE FLUX', xtitle = 'FLUX (W/M!U2!N)',     $
ytitle = 'PRESSURE (mbar)'
plot, firup, pres, xrange=[firmin,firmax], xs=1, $
/ylog, yrange = [ypmin,ypmax], $
symsize=0.6, psym = -1, $
title = 'UPWARD IR FLUX', xtitle = 'FLUX (W/M!U2!N)',     $
ytitle = 'PRESSURE (mbar)'
;
; Downward fluxes
;
readf,1, firmax, firmin, fvmax, fvmin
for n=0,num-1 do begin
  readf,1, p1, f1, f2
  pres(n)  = p1
  fvdn(n)  = f1
  firdn(n) = f2
endfor
if (fvmax le 0.1) then begin
  fvmin = -0.5
  fvmax = 0.5
endif
;
plot, fvdn, pres, xrange=[fvmin,fvmax], xs=1, $
/ylog, yrange = [ypmin,ypmax], $
symsize=0.6, psym = -4, $
title = 'DOWNWARD VISIBLE FLUX', xtitle = 'FLUX (W/M!U2!N)',     $
ytitle = 'PRESSURE (mbar)'
plot, firdn, pres, xrange=[firmin,firmax], xs=1, $
/ylog, yrange = [ypmin,ypmax], $
symsize=0.6, psym = -1, $
title = 'DOWNWARD IR FLUX', xtitle = 'FLUX (W/M!U2!N)',     $
ytitle = 'PRESSURE (mbar)'
;
; Net fluxes and other original DRIVER plots
;
readf,1, num
readf,1, cosz, albedo, solar, tautot, conrnu, tautotc
readf,1, firmax, firmin, fvmax, fvmin
for n=0,num-1 do begin
  readf, 1, p1,f1,f2
  pres(n) = p1
  fv(n)   = f1
  fir(n)  = f2
endfor
;
;  heating rates
;
readf,1,hmin, hmax
for n=0,num-1 do begin
  readf,1, p1, hv, hi, ht
  pres(n)   = p1
  heatv(n)  = hv
  heatir(n) = hi
  heatt(n)  = ht
endfor
;
;  temperature profile
for n=0,num-1 do begin
  readf, 1, p1,t1
  pres(n)  = p1
  temp(n)  = t1
endfor
readf, 1, gtemp, PSF
;      
if (abs(fvmax) le 0.1) then begin
  fvmin = -0.5
  fvmax = 0.5
endif
dx=0.1*(fvmax-fvmin)
fvmax = fvmax + dx
x2 = fvmin+2.5*dx
plot, fv, pres, xrange=[fvmin,fvmax], xs=1, $
/ylog, yrange = [ypmin,ypmax], $
symsize=0.6, psym = -4, $
title = 'Net Solar Flux', xtitle = 'FLUX (W/M!U2!N)', $
ytitle = 'PRESSURE (mbar)'
rheader=string(format='("SOLAR FLUX AT MARS = "(f5.1)," W/m!U2!N")',solar)
xyouts, x2, 1e-5, rheader, size=0.6
rheader=string(format='("COSZ = "(f4.2))',cosz)
xyouts, x2, 3.2e-5, rheader, size=0.6
rheader=string(format='("ALBEDO = "(f5.3))',albedo)
xyouts, x2, 1e-4, rheader, size=0.6
rheader=string(format='("Psf = "(f6.1)," mbar")',psf)
xyouts, x2, 3.2e-4, rheader, size=0.6
rheader=string(format='("CONRNU = "(f6.3))',conrnu)
xyouts, x2, 1.0e-3, rheader, size=0.6
rheader=string(format='("!4s!3!Ld!N = "(f7.2))',tautot)
xyouts, x2, 3.2e-3, rheader, size=0.6
rheader=string(format='("!4s!3!Lc!N = "(f7.2))',tautotc)
xyouts, x2, 1.0e-2, rheader, size=0.6
;
dx=0.1*(firmax-firmin)
firmin = firmin - 2.0*dx
plot, fir, pres, xrange=[firmin,firmax], xs=1, $
/ylog, yrange = [ypmin,ypmax], $
symsize=0.6, psym = -1, $
title = 'Net IR Flux', xtitle = 'FLUX (W/M!U2!N)', ytitle = 'PRESSURE (mbar)'
;
;  heating rates
;
hmin = -500.0
hmax =  500.0
plot, heatt, pres, xs=1, xrange=[hmin,hmax], $
/ylog, yrange = [ypmin,10*ypmax], title="Heating Rate (k/sol)", $
xtitle="Heating (K/sol)", ytitle="Pressure (mbar)"
oplot, heatv, pres, color=50
oplot, heatir, pres, color=240
;
xpp=[0.0,0.0]
ypp=[ypmin,ypmax]
oplot, xpp, ypp
;
dx=0.1*(hmax-hmin)
x2 = hmin+dx
xyouts,x2,1.0e-2,"VISIBLE", color=50, size=0.8
xyouts,x2,1.0e-1,"IR",color=240, size=0.8
xyouts,x2,1.0e00,"TOTAL", size=0.8
;
xyouts, hmin, 500.0, systime(0)+'        driver.pro', size=0.6
;
;  temperature profile
;
plot, temp, pres, xrange=[100,300], xs=1,  $
/ylog, yrange = [ypmin,ypmax], $
 title = "Temperature Profile", xtitle="Temperature (K)", $
ytitle = 'PRESSURE (mbar)'
usx=[-0.5, 0.0, 0.5, -0.5]
usy=[1, 0, 1, 1]
USERSYM, usx, usy, /fill
xxp=fltarr(1)
yyp=fltarr(1)
xxp(0) = gtemp
yyp(0) = ypmin
oplot, xxp, yyp, psym=-8, symsize=2
;
end
