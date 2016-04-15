; Save margins
x_margin = !x.margin & y_margin = !y.margin & z_margin = !z.margin

; Read files
filename1 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_h125_1.00_1000_0.80_GC_upwelling_output.nc'
filename2 = '/data/tempo2/ggonzale/GEOCAPE-TOOL/H2O/CALCULATIONS/H2O_h125_1.00_0200_0.80_GC_upwelling_output.nc'

obj   = Obj_New('NCDF_DATA', filename1)
data1 = obj -> ReadFile(filename1)

obj   = Obj_New('NCDF_DATA', filename2)
data2 = obj -> ReadFile(filename2)

; Create layer pressure and altitude
zmid = reverse(zmid(data1.zs.data))
pmid = reverse(pmid(data1.ps.data))

; Save SZA, VZA, and RAA
sza = data1.solarzenithangle.data
vza = data1.viewingzenithangle.data
raa = data1.relativeazimuthangle.data

sw1 = data1.scatweights.data
sw2 = data2.scatweights.data

nsza = data1._dimensions.nsza
nvza = data1._dimensions.nvza
nraa = data1._dimensions.naza
nlay = data1._dimensions.nlayer
sw1m = FLTARR(nlay,nsza,nvza,nraa)
sw2m = FLTARR(nlay,nsza,nvza,nraa)

For isza = 0, nsza-1 do begin
For ivza = 0, nvza-1 do begin
For iraa = 0, nraa-1 do begin
For ilay = 0, nlay-1 do begin 
 
   sw1m[ilay,isza,ivza,iraa] = sw1[ilay,iraa+ivza*nraa+isza*nraa*nvza]
   sw2m[ilay,isza,ivza,iraa] = sw2[ilay,iraa+ivza*nraa+isza*nraa*nvza]

Endfor
Endfor
Endfor
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
