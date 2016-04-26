; Save margins
x_margin = !x.margin & y_margin = !y.margin & z_margin = !z.margin

; Read files
filename1 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_1030_0.00_1000_0.01_GC_upwelling_output.nc'
filename2 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0950_0.00_1000_0.01_GC_upwelling_output.nc'
filename3 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0900_0.00_1000_0.01_GC_upwelling_output.nc'
filename4 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0850_0.00_1000_0.01_GC_upwelling_output.nc'
filename5 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0800_0.00_1000_0.01_GC_upwelling_output.nc'
filename6 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0750_0.00_1000_0.01_GC_upwelling_output.nc'
filename7 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0700_0.00_1000_0.01_GC_upwelling_output.nc'
filename8 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0650_0.00_1000_0.01_GC_upwelling_output.nc'
filename9 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0600_0.00_1000_0.01_GC_upwelling_output.nc'
filename10 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0550_0.00_1000_0.01_GC_upwelling_output.nc'
filename11 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0500_0.00_1000_0.01_GC_upwelling_output.nc'
filename12 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_l325_0450_0.00_1000_0.01_GC_upwelling_output.nc'

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

surf  = [1030.0, 950.0, 900.0, 850.0, 800.0, 750.0, 700.0, 650.0, 600.0, 550.0, 500.0, 450.0]
; Plot radiance albedo dependence
rad  = [in1m[0,0,0],in2m[0,0,0],in3m[0,0,0],in4m[0,0,0],in5m[0,0,0],in6m[0,0,0],$
        in7m[0,0,0],in8m[0,0,0],in9m[0,0,0],in10m[0,0,0],in11m[0,0,0],in12m[0,0,0]]
cgplot, alog(surf), rad, color = 1, yrange = [0.0, MAX(in12m[*,*,*])]

rad  = [in1m[5,5,5],in2m[5,5,5],in3m[5,5,5],in4m[5,5,5],in5m[5,5,5],in6m[5,5,5],$
        in7m[5,5,5],in8m[5,5,5],in9m[5,5,5],in10m[5,5,5],in11m[5,5,5],in12m[5,5,5]]
cgplot, alog(surf), rad, color = 3, /Overplot

; Plot scattering albedo dependence
sca  = [sw1m[45,0,0,0],sw2m[45,0,0,0],sw3m[45,0,0,0],sw4m[45,0,0,0],sw5m[45,0,0,0],sw6m[45,0,0,0],$
        sw7m[45,0,0,0],sw8m[45,0,0,0],sw9m[45,0,0,0],sw10m[45,0,0,0],sw11m[45,0,0,0],sw12m[45,0,0,0]]
cgplot, alog(surf), sca, color = 1, yrange = [0.0, MAX(sw12m[*,*,*,*])]

sca  = [sw1m[20,5,5,5],sw2m[20,5,5,5],sw3m[20,5,5,5],sw4m[20,5,5,5],sw5m[20,5,5,5],sw6m[20,5,5,5],$
        sw7m[20,5,5,5],sw8m[20,5,5,5],sw9m[20,5,5,5],sw10m[20,5,5,5],sw11m[20,5,5,5],sw12m[20,5,5,5]]
cgplot, alog(surf), sca, color = 3, /Overplot

END
