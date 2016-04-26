; Save margins
x_margin = !x.margin & y_margin = !y.margin & z_margin = !z.margin

; Read files
filename1 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.0001_GC_upwelling_output.nc'
filename2 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.01_GC_upwelling_output.nc'
filename3 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.10_GC_upwelling_output.nc'
filename4 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.20_GC_upwelling_output.nc'
filename5 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.30_GC_upwelling_output.nc'
filename6 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.40_GC_upwelling_output.nc'
filename7 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.50_GC_upwelling_output.nc'
filename8 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.60_GC_upwelling_output.nc'
filename9 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.70_GC_upwelling_output.nc'
filename10 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.80_GC_upwelling_output.nc'
filename11 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.90_GC_upwelling_output.nc'
filename12 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_1.00_GC_upwelling_output.nc'

obj   = Obj_New('NCDF_DATA', filename1)
data1 = obj -> ReadFile(filename1)
obj   = Obj_New('NCDF_DATA', filename2)
data2 = obj -> ReadFile(filename2)
obj   = Obj_New('NCDF_DATA', filename3)
data3 = obj -> ReadFile(filename3)
obj   = Obj_New('NCDF_DATA', filename4)
data4 = obj -> ReadFile(filename4)
obj   = Obj_New('NCDF_DATA', filename5)
data5 = obj -> ReadFile(filename5)
obj   = Obj_New('NCDF_DATA', filename6)
data6 = obj -> ReadFile(filename6)
obj   = Obj_New('NCDF_DATA', filename7)
data7 = obj -> ReadFile(filename7)
obj   = Obj_New('NCDF_DATA', filename8)
data8 = obj -> ReadFile(filename8)
obj   = Obj_New('NCDF_DATA', filename9)
data9 = obj -> ReadFile(filename9)
obj    = Obj_New('NCDF_DATA', filename10)
data10 = obj -> ReadFile(filename10)
obj    = Obj_New('NCDF_DATA', filename11)
data11 = obj -> ReadFile(filename11)
obj    = Obj_New('NCDF_DATA', filename12)
data12 = obj -> ReadFile(filename12)

; Create layer pressure and altitude
zmid = reverse(zmid(data1.zs.data))
pmid = reverse(pmid(data1.ps.data))

; Save SZA, VZA, and RAA
sza = data1.solarzenithangle.data
vza = data1.viewingzenithangle.data
raa = data1.relativeazimuthangle.data

sw1 = data1.scatweights.data & in1 = data1.radiance.data
sw2 = data2.scatweights.data & in2 = data2.radiance.data
sw3 = data3.scatweights.data & in3 = data3.radiance.data
sw4 = data4.scatweights.data & in4 = data4.radiance.data
sw5 = data5.scatweights.data & in5 = data5.radiance.data
sw6 = data6.scatweights.data & in6 = data6.radiance.data
sw7 = data7.scatweights.data & in7 = data7.radiance.data
sw8 = data8.scatweights.data & in8 = data8.radiance.data
sw9 = data9.scatweights.data & in9 = data9.radiance.data
sw10 = data10.scatweights.data & in10 = data10.radiance.data
sw11 = data11.scatweights.data & in11 = data11.radiance.data
sw12 = data12.scatweights.data & in12 = data12.radiance.data

nsza = data1._dimensions.nsza
nvza = data1._dimensions.nvza
nraa = data1._dimensions.naza
nlay = data1._dimensions.nlayer

sw1m  = FLTARR(nlay,nsza,nvza,nraa) & in1m = FLTARR(nsza,nvza,nraa)
sw2m  = FLTARR(nlay,nsza,nvza,nraa) & in2m = FLTARR(nsza,nvza,nraa)
sw3m  = FLTARR(nlay,nsza,nvza,nraa) & in3m = FLTARR(nsza,nvza,nraa)
sw4m  = FLTARR(nlay,nsza,nvza,nraa) & in4m = FLTARR(nsza,nvza,nraa)
sw5m  = FLTARR(nlay,nsza,nvza,nraa) & in5m = FLTARR(nsza,nvza,nraa)
sw6m  = FLTARR(nlay,nsza,nvza,nraa) & in6m = FLTARR(nsza,nvza,nraa)
sw7m  = FLTARR(nlay,nsza,nvza,nraa) & in7m = FLTARR(nsza,nvza,nraa)
sw8m  = FLTARR(nlay,nsza,nvza,nraa) & in8m = FLTARR(nsza,nvza,nraa)
sw9m  = FLTARR(nlay,nsza,nvza,nraa) & in9m = FLTARR(nsza,nvza,nraa)
sw10m = FLTARR(nlay,nsza,nvza,nraa) & in10m = FLTARR(nsza,nvza,nraa)
sw11m = FLTARR(nlay,nsza,nvza,nraa) & in11m = FLTARR(nsza,nvza,nraa)
sw12m = FLTARR(nlay,nsza,nvza,nraa) & in12m = FLTARR(nsza,nvza,nraa)

For isza = 0, nsza-1 do begin
For ivza = 0, nvza-1 do begin
For iraa = 0, nraa-1 do begin

   in1m[isza,ivza,iraa] = in1[iraa+ivza*nraa+isza*nraa*nvza]
   in2m[isza,ivza,iraa] = in2[iraa+ivza*nraa+isza*nraa*nvza]
   in3m[isza,ivza,iraa] = in3[iraa+ivza*nraa+isza*nraa*nvza]
   in4m[isza,ivza,iraa] = in4[iraa+ivza*nraa+isza*nraa*nvza]
   in5m[isza,ivza,iraa] = in5[iraa+ivza*nraa+isza*nraa*nvza]
   in6m[isza,ivza,iraa] = in6[iraa+ivza*nraa+isza*nraa*nvza]
   in7m[isza,ivza,iraa] = in7[iraa+ivza*nraa+isza*nraa*nvza]
   in8m[isza,ivza,iraa] = in8[iraa+ivza*nraa+isza*nraa*nvza]
   in9m[isza,ivza,iraa] = in9[iraa+ivza*nraa+isza*nraa*nvza]
   in10m[isza,ivza,iraa] = in10[iraa+ivza*nraa+isza*nraa*nvza]
   in11m[isza,ivza,iraa] = in11[iraa+ivza*nraa+isza*nraa*nvza]
   in12m[isza,ivza,iraa] = in12[iraa+ivza*nraa+isza*nraa*nvza]

For ilay = 0, nlay-1 do begin 
 
   sw1m[ilay,isza,ivza,iraa] = sw1[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw2m[ilay,isza,ivza,iraa] = sw2[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw3m[ilay,isza,ivza,iraa] = sw3[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw4m[ilay,isza,ivza,iraa] = sw4[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw5m[ilay,isza,ivza,iraa] = sw5[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw6m[ilay,isza,ivza,iraa] = sw6[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw7m[ilay,isza,ivza,iraa] = sw7[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw8m[ilay,isza,ivza,iraa] = sw8[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw9m[ilay,isza,ivza,iraa] = sw9[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw10m[ilay,isza,ivza,iraa] = sw10[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw11m[ilay,isza,ivza,iraa] = sw11[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw12m[ilay,isza,ivza,iraa] = sw12[ilay,iraa+ivza*nraa+isza*nraa*nvza]


Endfor
Endfor
Endfor
Endfor

; Plot radiance SZA dependence
cgplot, cos(sza*!pi/180.0), in2m[*,0,0], color = 1, yrange = [0.0, MAX(in12m)]
For ivza = 0, nvza-1 do begin
   cgplot, cos(sza*!pi/180.0), in2m[*,ivza,0], color = 1+ivza, /Overplot
   cgplot, cos(sza*!pi/180.0), in10m[*,ivza,0], color = 1+ivza, /Overplot
Endfor

; Plot radiance VZA dependence
cgplot, cos(vza*!pi/180.0), in2m[0,*,0], color = 1, yrange = [0.0, MAX(in12m)]
For isza = 0, nsza-1 do begin
   cgplot, cos(vza*!pi/180.0), in2m[isza,*,0], color = 1+isza, /Overplot
   cgplot, cos(vza*!pi/180.0), in10m[isza,*,0], color = 1+isza, /Overplot
Endfor

diff = (sw1m-sw2m)
cgplot, diff[*,0,0,0], zmid, color = 1, /Nodata, xrange = [min(diff),max(diff)]
For isza = 0, nsza-1 do begin
For ivza = 0, nvza-1 do begin
For iraa = 0, nraa-1 do begin

   cgplot, diff[*,isza,ivza,iraa], zmid, color = 1, /Overplot

Endfor
Endfor
Endfor

END
