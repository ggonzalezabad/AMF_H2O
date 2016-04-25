PRO fit_toms_clouds

; We are going to try to fit the radiances to the expression in the OMI algorithm 
; technical note:
; I = I0(sza, vza) + I1(sza, vza)cos(RAA) + I2(sza, vza) cos(2(RAA)) + ((RIr(sza, vza))/(1-RSb))
; where R is the albedo. There are therefore five fitting parameters:
; I0, I1, I2, Ir and Sb and two variables RAA and R.
; Wavelenghts below 316 nm with high SZA and VZA angles don't fit
; well to the expression!!

; To prevent the procedure to use more than 100% CPU in one node
;CPU, TPOOL_NTHREADS = 1
; To know if a regular IDL license is used
;IF lmgr(/runtime) THEN message, /info, 'runtime mode is on' ELSE message, 'not in runtime mode'
; Folder to output the results of the fitting
DIR_FITTING = '../fitting_results_cloud'

; Variables to take account of the ozone profiles, lower level pressure, solar zenith angles, albedos,
; viewing zenith angles, relative azimuth angles.
; Ozone profiles
NOZO = 26
TOMS =['h125', 'h275', 'h425', 'h575', 'l325', 'l475', 'm225', 'm375', 'm525', 'h175', $
       'h325', 'h475', 'l225', 'l375', 'm125', 'm275', 'm425', 'm575', 'h225', 'h375', $
       'h525', 'l275', 'l425', 'm175', 'm325', 'm475']

; Lower level pressure (surface)
NPRE = 12
SURF =['1030', '0950', '0900', '0850', '0800', '0750', '0700', '0650', '0600', '0550', '0500', '0450']

; Lower level pressure (Cloud pressure)
NCLD = 12
CLDP =['1000', '0950', '0900', '0850', '0800', '0750', '0700', '0600', '0500', '0400', '0300', '0200']

; Read one file to get the grid (these should not change for different
; viewing geometries, albedo and cloud properties)
filename = '../../CALCULATIONS/H2O_h125_1030_1.00_1000_0.80_GC_upwelling_output.nc'
obj  = Obj_New('NCDF_DATA', filename)
data = Obj -> ReadFile(filename)
; Read number of sza, vza, raa, layers and wavelengths
nsza = data._dimensions.nsza
nvza = data._dimensions.nvza
nraa = data._dimensions.naza
nlev = data._dimensions.nlevel
nlay = data._dimensions.nlayer
nwav = data._dimensions.nw

; Read SZA, VZA, and RAA
sza = data.solarzenithangle.data              ; [degrees]
vza = data.viewingzenithangle.data            ; [degrees]
raa = data.relativeazimuthangle.data          ; [degrees]
wav = data.wavelength.data                    ; [nm]

; Loop over different files. I need to take into account the different
; inputs. The dependence with albedo for clouds is dropped since we
; consider clouds as Lambertian surface with a reflectance of 0.8.
For iozo = 0, nozo-1 do begin
   print, nozo-iozo
   For ipre = 0, npre-1 do begin
      read_prof = 1
      For icld = 0, ncld-1 do begin
         ; If the cloud is below the surface then skip this loop
         IF ( FLOAT(CLDP[icld]) GT FLOAT(SURF[ipre]) ) THEN CONTINUE

         radiance = FLTARR(nsza,nvza,nraa)
         scatteri = FLTARR(nsza,nvza,nraa,nlay)
         filename = '../../CALCULATIONS/H2O_'+TOMS[iozo]+'_'+ $
                    SURF[ipre]+'_1.00_'+CLDP[icld]+'_0.80_GC_upwelling_output.nc'
         print, filename
   ; Open and read file
         obj  = Obj_New('NCDF_DATA', filename)
         data = Obj -> ReadFile(filename)
   ; Read atmospheric profiles (aircolumn, gascolumn, temperature)
         IF (read_prof) THEN BEGIN
   ; Read pressure, altitude and wavelength
            zlev = data.zs.data    ; Level [km]
            plev = data.ps.data    ; Level [hPa]
   ; Create layer pressure and altitude
            zmid = reverse(zmid(data.zs.data)) ; Layer [km]
            pmid = reverse(pmid(data.ps.data)) ; Layer [hPa]
   ; Read air column, gas column and temperature
            airc = data.aircol.data ; Layer [molecules/cm-2]
            gasc = data.gascol.data ; Layer [molecules/cm-2]
            temp = data.ts.data     ; Level [K]
   ; Don't read profiles again
            read_prof = 0
         ENDIF
      ; Fill up radiance, jacobians, and scattering weights matrices
         For isza = 0, nsza-1 do begin
            For ivza = 0, nvza-1 do begin
               For iraa = 0, nraa-1 do begin
                  radiance[isza,ivza,iraa] = data.radiance.data[iraa+ivza*nraa+isza*nraa*nvza]
                  For ilay = 0, nlay-1 do begin 
                     scatteri[isza,ivza,iraa,ilay] = data.scatweights.data[ilay,iraa+ivza*nraa+isza*nraa*nvza]
                  Endfor
               Endfor
            Endfor
         Endfor

         radiance_bis = radiance
   ; Now that I have all the data I need I can carry on with the parameterization.
   ; I'm going to use these two again and again
         COSRAA  = COS(RAA*!dtor)
         COS2RAA = COS(2*RAA*!dtor)
         PRINT, STRCOMPRESS('OZONE '+TOMS[iozo]+', NWAV '+STRING(NWAV)+', NLAYERS '+STRING(NLAY)+', NLEVELS '+$
                            STRING(NLEV)+' SURFACE LEVEL '+SURF[ipre]+', CLOUD PRESSURE '+CLDP[icld])
   ; Variable to keep the calculated parameters.
         I0        = FLTARR(NSZA,NVZA)
         I1        = FLTARR(NSZA,NVZA)
         I2        = FLTARR(NSZA,NVZA)
         Ir        = 0
         Sb        = 0

   ; Starting the fitting of the curve. Using CURVEFIT.
   ; The X values are RAA and Y are ALB. First I want to split the VZA and RAA variables.
   ;-------------------------------------------------------------------------
   ; First the fitting for the IO, I1 and I2 using ALBEDO = 0 calculations.
         FOR II = 0, NSZA-1 DO BEGIN
            FOR JJ = 0, NVZA-1 DO BEGIN

               AI0I1I2 = [0.004,-0.0002,0.0003]
               EXPRAD = REFORM(RADIANCE[II,JJ,*])
               test = FLTARR(nraa,2)
               test[*,0] = raa
               test[*,1] = 0.80
               Weights = 1/REFORM(RADIANCE[II,JJ,*])
               yfit = CURVEFIT(raa, EXPRAD, Weights, AI0I1I2, FUNCTION_NAME = 'TOMRAD_I0I1I2', $
                               ITER = ITER, ITMAX = 30, CHISQ = CHISQ, YERROR = YERROR, STATUS = STATUS, /DOUBLE)

               I0[II,JJ]          = AI0I1I2[0]
               I1[II,JJ]          = AI0I1I2[1]
               I2[II,JJ]          = AI0I1I2[2]

            ENDFOR
         ENDFOR

   ; Saving the results for the radiance fit
         TOMS_CLI = TOMS[iozo]
         SAVE, nsza, nvza, nraa, nlev, nlay, nwav, zlev, plev,   $
               zmid, pmid, sza, vza, raa, wav, airc, gasc, temp, $
               I0, I1, I2, Ir, Sb, TOMS_CLI, $
               FILENAME = DIR_FITTING+'/Radiance_fit_'+TOMS[iozo]+'_'+SURF[ipre]+'_cloud_'+CLDP[icld]+'.dat'
   
         For isza = 0, nsza-1 do begin
            For ivza = 0, nvza-1 do begin
               For iraa = 0, nraa-1 do begin
                  radiance_bis[isza,ivza,iraa] = $
                     I0[isza,ivza] + $
                     I1[isza,ivza] * cosRAA[iraa]  + $
                     I2[isza,ivza] * cos2RAA[iraa] + $
                     ( 0.8 * Ir ) / ( 1 - ( 0.8 * Sb ) ) 
               Endfor
            Endfor
         Endfor

         test_rad = (radiance-radiance_bis)/radiance*100.0

   ; Fitting of the Jacobian
         dI0  = DBLARR(NSZA,NVZA,NLAY)
         dI1  = DBLARR(NSZA,NVZA,NLAY)
         dI2  = DBLARR(NSZA,NVZA,NLAY)

         scatteri_bis = scatteri
         NCOUNT   = 0.0
         ACCEPTED = 0.0
         ACCEPTED = 0.0
         NSTATUS  = 0.0
         For isza = 0, 8 do begin  ;nsza-1 do begin
            For ivza = 0, 7 do begin ;0, nvza-1 do begin
               For ilay = 0, 46 do begin ;0, 46 do begin ;nlay-1 do begin
                  AI0I1I2               = [0.5D, 0.05D, 0.5D0]
                  EXPJAC                = REFORM(SCATTERI[isza,ivza,*,ilay])
                  IF TOTAL(EXPJAC) LT 0.0001 THEN CONTINUE
                  Weights               = 1.0d/EXPJAC
                  FITA = [1,1,1]
                  yfit = CURVEFIT(RAA, EXPJAC, Weights, AI0I1I2, FUNCTION_NAME = 'TOMRAD_dI0I1I2', $
                                  ITER=ITER, ITMAX = 30, CHISQ = CHISQ, YERROR = YERROR,           $
                                  STATUS = STATUS, FITA = FITA, /DOUBLE)
                  
                  dI0[isza,ivza,ilay] = AI0I1I2[0] 
                  dI1[isza,ivza,ilay] = AI0I1I2[1] 
                  dI2[isza,ivza,ilay] = AI0I1I2[2] 
                  scatteri_bis[isza,ivza,*,ilay] = dI0[isza,ivza,ilay] + $
                                                   dI1[isza,ivza,ilay] * COSRAA + $
                                                   dI2[isza,ivza,ilay] * COS2RAA
                  
               endfor
            endfor
         endfor
; Saving the results of the Scattering fittings. dI0, dI1, dI2
         SAVE, nsza, nvza, nraa, nlev, nlay, nwav, zlev, plev, $
               zmid, pmid, sza, vza, raa, wav, airc, gasc, temp, $
               dI0, dI1, dI2, TOMS_CLI, $
               FILENAME = DIR_FITTING+'/Scattering_fit_'+TOMS[iozo]+'_'+SURF[ipre]+'_cloud_'+CLDP[icld]+'.dat'   
         test_scat = (scatteri-scatteri_bis)/scatteri*100.0         
         print, max(test_rad), min(test_rad), max(test_scat), min(test_scat)

      Endfor ; End cloud pressure loop
   Endfor ; End surface pessure loop
Endfor ; End TOMS loop
; Just in case a file is left open
CLOSE, /ALL
  
END


PRO TOMRAD_I0I1I2, RAA, P, YMOD, PDER

  X = RAA

  ; The RAA dependence
  COSRAA  = COS(X*!dtor)
  COS2RAA = COS(2*X*!dtor)
    
  ; The function
  YMOD = P[0] + P[1] * COSRAA + P[2] * COS2RAA
  
  ; The derivatives P[0], P[1] and P[2]  
  IF N_PARAMS() GE 3 THEN Pder = [[REPLICATE(1,N_ELEMENTS(X))], [COSRAA], [COS2RAA]]
   
END

PRO TOMRAD_dI0I1I2, RAA, P, YMOD, PDER

  X = RAA[*,0]

  ; The RAA dependence
  COSRAA  = COS(X*!dtor)
  COS2RAA = COS(2*X*!dtor)
    
  ; The function
  YMOD = P[0] + $
         P[1] * COSRAA  + $
         P[2] * COS2RAA
  
  ; The derivatives P[0], P[1] and P[2]  
  IF N_PARAMS() GE 4 THEN Pder = [[REPLICATE(1,N_ELEMENTS(X))], $
                                  [COSRAA],  $
                                  [COS2RAA]]
   
END
