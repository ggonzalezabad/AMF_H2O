PRO create_h2o_h5_table

  FILENAME         = 'AMFs_lookup_h2o_442nm_test.he5'
  Author           = 'Gonzalo Gonzalez Abad ggonzalezabad@cfa.harvard.edu'
  Institution      = 'Smithsonian Astrophysical Observatory'
  Last_file_update = SYSTIME()
  ; Create file
  file_id = H5F_CREATE(filename)

  ; Fill in global attributes
  datatype  = H5T_IDL_CREATE(Author)
  dataspace = H5S_CREATE_SIMPLE(1)
  att_id  = H5A_CREATE(file_id, 'Author', datatype, dataspace)
  H5A_WRITE, att_id, Author
  H5S_CLOSE, dataspace
  H5T_CLOSE, datatype
  H5A_CLOSE, att_id
  datatype  = H5T_IDL_CREATE(Institution)
  dataspace = H5S_CREATE_SIMPLE(1)
  att_id  = H5A_CREATE(file_id, 'Institution', datatype, dataspace)
  H5A_WRITE, att_id, Institution
  H5S_CLOSE, dataspace
  H5T_CLOSE, datatype
  H5A_CLOSE, att_id
  datatype  = H5T_IDL_CREATE(Last_file_update)
  dataspace = H5S_CREATE_SIMPLE(1)
  att_id  = H5A_CREATE(file_id, 'Last file update', datatype, dataspace)
  H5A_WRITE, att_id, Last_file_update
  H5S_CLOSE, dataspace
  H5T_CLOSE, datatype
  H5A_CLOSE, att_id

  ; Create groups
  grp_id = H5G_CREATE(file_id, 'Grid')
  H5G_CLOSE, grp_id
  grp_id = H5G_CREATE(file_id, 'Intensity')
  H5G_CLOSE, grp_id
  grp_id = H5G_CREATE(file_id, 'Scattering Weights')
  H5G_CLOSE, grp_id
  grp_id = H5G_CREATE(file_id, 'Profiles')
  H5G_CLOSE, grp_id

  ; Fill in Grid group
  filename = '../fitting_results_clear/Radiance_fit_l375_1030_clear.dat'
  RESTORE, filename
  grp_id =H5G_OPEN(file_id,'Grid')
  ; Wavelength
  datatype  = H5T_IDL_CREATE(wav)
  dataspace = H5S_CREATE_SIMPLE(1)
  dat_id    = H5D_CREATE(grp_id, 'Wavelength', datatype, dataspace)
  H5D_WRITE, dat_id, wav
  unit   = '[nm]' & datatype = H5T_IDL_CREATE(unit)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; VZA
  datatype  = H5T_IDL_CREATE(vza)
  dataspace = H5S_CREATE_SIMPLE(SIZE(vza,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'VZA', datatype, dataspace)
  H5D_WRITE, dat_id, vza
  unit   = '[Degrees]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; SZA
  datatype  = H5T_IDL_CREATE(sza)
  dataspace = H5S_CREATE_SIMPLE(SIZE(sza,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'SZA', datatype, dataspace)
  H5D_WRITE, dat_id, sza
  unit   = '[Degrees]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Albedo
  albedo = [0.0001,0.01,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00]
  nalbedo = N_ELEMENTS(albedo)
  datatype  = H5T_IDL_CREATE(albedo)
  dataspace = H5S_CREATE_SIMPLE(SIZE(albedo,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Albedo', datatype, dataspace)
  H5D_WRITE, dat_id, albedo
  unit   = 'Unitless' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Cloud pressure
  cldp = [1000.0,0950.0,0900.0,0850.0,0800.0,0750.0,0700.0,0600.0,0500.0,0400.0,0300.0,0200.0]
  cldps = ['1000','0950','0900','0850','0800','0750','0700','0600','0500','0400','0300','0200']
  ncldp = N_ELEMENTS(cldp)
  datatype  = H5T_IDL_CREATE(cldp)
  dataspace = H5S_CREATE_SIMPLE(SIZE(cldp,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Cloud Pressure', datatype, dataspace)
  H5D_WRITE, dat_id, cldp
  unit   = '[hPa]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Surface pressure
  surf = [1030.0,0950.0,0900.0,0850.0,0800.0,0750.0,0700.0,0650.0,0600.0,0550.0,0500.0,0450.0]
  surfs = ['1030','0950','0900','0850','0800','0750','0700','0650','0600','0550','0500','0450']
  nsurf = N_ELEMENTS(surf)
  datatype  = H5T_IDL_CREATE(surf)
  dataspace = H5S_CREATE_SIMPLE(SIZE(surf,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Surface Pressure', datatype, dataspace)
  H5D_WRITE, dat_id, surf
  unit   = '[hPa]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; TOMS profile
  toms = ['h125', 'h175', 'h225', 'h275', 'h325', 'h375', 'h425', 'h475', 'h525', 'h575', $
          'l225', 'l275', 'l325', 'l375', 'l425', 'l475', $
          'm125', 'm175', 'm225', 'm275', 'm325', 'm375', 'm425', 'm475', 'm525', 'm575']
  ntoms = N_ELEMENTS(toms)
  datatype  = H5T_IDL_CREATE(toms)
  dataspace = H5S_CREATE_SIMPLE(SIZE(toms,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'TOMS', datatype, dataspace)
  H5D_WRITE, dat_id, toms
  unit   = '[String]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  H5G_CLOSE, grp_id

  ; Fill in profiles group
  grp_id =H5G_OPEN(file_id,'Profiles')
  ; Create variables to hold data
  altitude_level = FLTARR(NSURF,NLEV)
  altitude_layer = FLTARR(NSURF,NLAY)
  pressure_level = FLTARR(NSURF,NLEV)
  pressure_layer = FLTARR(NSURF,NLAY)
  temperat_level = FLTARR(NTOMS,NSURF,NLEV)
  aircolum_layer = FLTARR(NTOMS,NSURF,NLAY)
  ozonecol_layer = FLTARR(NTOMS,NSURF,NLAY)
  For isurf = 0, nsurf-1 do begin
     FOR itoms = 0, ntoms-1 do begin
        filename = '../fitting_results_clear/Radiance_fit_'+$
                   TOMS[itoms]+'_'+SURFS[isurf]+'_clear.dat'
        RESTORE, filename
        altitude_level[isurf,*] = zlev
        altitude_layer[isurf,*] = zmid
        pressure_level[isurf,*] = plev
        pressure_layer[isurf,*] = pmid
        temperat_level[itoms,isurf,*] = temp
        aircolum_layer[itoms,isurf,*] = airc
        ozonecol_layer[itoms,isurf,*] = gasc
     ENDFOR
  ENDFOR

  ; Pressure level
  datatype  = H5T_IDL_CREATE(pressure_level)
  dataspace = H5S_CREATE_SIMPLE(SIZE(pressure_level,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Pressure Level', datatype, dataspace)
  H5D_WRITE, dat_id, pressure_level
  unit   = '[hPa]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Pressure layer
  datatype  = H5T_IDL_CREATE(pressure_layer)
  dataspace = H5S_CREATE_SIMPLE(SIZE(pressure_layer,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Pressure Layer', datatype, dataspace)
  H5D_WRITE, dat_id, pressure_layer
  unit   = '[hPa]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Altitude level
  datatype  = H5T_IDL_CREATE(altitude_level)
  dataspace = H5S_CREATE_SIMPLE(SIZE(altitude_level,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Altitude Level', datatype, dataspace)
  H5D_WRITE, dat_id, altitude_level
  unit   = '[km]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Altitude layer
  datatype  = H5T_IDL_CREATE(altitude_layer)
  dataspace = H5S_CREATE_SIMPLE(SIZE(altitude_layer,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Altitude Layer', datatype, dataspace)
  H5D_WRITE, dat_id, altitude_layer
  unit   = '[km]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Temperature (level)
  datatype  = H5T_IDL_CREATE(temperat_level)
  dataspace = H5S_CREATE_SIMPLE(SIZE(temperat_level,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Temperature Level', datatype, dataspace)
  H5D_WRITE, dat_id, temperat_level
  unit   = '[K]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Air columns (layer)
  datatype  = H5T_IDL_CREATE(aircolum_layer)
  dataspace = H5S_CREATE_SIMPLE(SIZE(aircolum_layer,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Air Column Layer', datatype, dataspace)
  H5D_WRITE, dat_id, aircolum_layer
  unit   = '[mol/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Ozone columns (layer)
  datatype  = H5T_IDL_CREATE(ozonecol_layer)
  dataspace = H5S_CREATE_SIMPLE(SIZE(ozonecol_layer,/dimensions))
  dat_id    = H5D_CREATE(grp_id, 'Ozone Column Layer', datatype, dataspace)
  H5D_WRITE, dat_id, ozonecol_layer
  unit   = '[mol/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  H5G_CLOSE, grp_id

  ; Fill in intensity group
  grp_id = H5G_OPEN(file_id,'Intensity')
  ; Create clear sky group
  grp2_id = H5G_CREATE(grp_id, 'Clear Sky')
  H5G_CLOSE, grp2_id
  ; Create cloudy sky group
  grp2_id = H5G_CREATE(grp_id, 'Cloud Sky')
  H5G_CLOSE, grp2_id

  ; Create variables to hold data (clear sky)
  I0_all = FLTARR(NTOMS,NSURF,NVZA,NSZA)
  I1_all = FLTARR(NTOMS,NSURF,NVZA,NSZA)
  I2_all = FLTARR(NTOMS,NSURF,NVZA,NSZA)
  Ir_all = FLTARR(NTOMS,NSURF,NVZA,NSZA)
  Sb_all = FLTARR(NTOMS,NSURF)
  ; Fill in clear sky variables
  For itoms = 0, ntoms-1 do begin
     For isurf = 0, nsurf-1 do begin
        filename = '../fitting_results_clear/Radiance_fit_'+$
                   TOMS[itoms]+'_'+SURFS[isurf]+'_clear.dat'
        print, filename
        RESTORE, filename
        I0_all[itoms,isurf,*,*] = TRANSPOSE(I0)
        I1_all[itoms,isurf,*,*] = TRANSPOSE(I1)
        I2_all[itoms,isurf,*,*] = TRANSPOSE(I2)
        Ir_all[itoms,isurf,*,*] = TRANSPOSE(Ir)
        Sb_all[itoms,isurf]     = Sb
     Endfor
  Endfor
  ; Clear sky I0 (NSZA,NVZA,NTOMS)
  grp2_id = H5G_OPEN(grp_id, 'Clear Sky')
  datatype  = H5T_IDL_CREATE(I0_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(I0_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'I0', datatype, dataspace)
  H5D_WRITE, dat_id, I0_all
  unit   = '[W/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky I1 (NSZA,NVZA,NTOMS)
  datatype  = H5T_IDL_CREATE(I1_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(I1_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'I1', datatype, dataspace)
  H5D_WRITE, dat_id, I1_all
  unit   = '[W/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky I2 (NSZA,NVZA,NTOMS)
  datatype  = H5T_IDL_CREATE(I2_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(I2_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'I2', datatype, dataspace)
  H5D_WRITE, dat_id, I2_all
  unit   = '[W/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky Ir (NSZA,NVZA,NTOMS)
  datatype  = H5T_IDL_CREATE(Ir_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(Ir_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'Ir', datatype, dataspace)
  H5D_WRITE, dat_id, Ir_all
  unit   = '[W/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky Sb (NTOMS)
  datatype  = H5T_IDL_CREATE(Sb_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(Sb_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'Sb', datatype, dataspace)
  H5D_WRITE, dat_id, Sb_all
  unit   = '[Unitless]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  H5G_CLOSE,grp2_id

  ; Create variables to hold data (cloud sky)
  I0_all = FLTARR(NTOMS,NSURF,NCLDP,NVZA,NSZA)
  I1_all = FLTARR(NTOMS,NSURF,NCLDP,NVZA,NSZA)
  I2_all = FLTARR(NTOMS,NSURF,NCLDP,NVZA,NSZA)
  Ir_all = FLTARR(NTOMS,NSURF,NCLDP,NVZA,NSZA)
  Sb_all = FLTARR(NTOMS,NSURF,NCLDP)
  ; Fill in cloudy sky variables
  For itoms = 0, ntoms-1 do begin
     For isurf = 0, nsurf-1 do begin
        For ipre = 0, ncldp-1 do begin
           IF ( CLDP[ipre] GT SURF[isurf] ) THEN CONTINUE
           filename = '../fitting_results_cloud/Radiance_fit_'+$
                      TOMS[itoms]+'_'+SURFS[isurf]+'_cloud_'+CLDPS[ipre]+'.dat'
           print, filename
           RESTORE, filename
           I0_all[itoms,isurf,ipre,*,*] = TRANSPOSE(I0)
           I1_all[itoms,isurf,ipre,*,*] = TRANSPOSE(I1)
           I2_all[itoms,isurf,ipre,*,*] = TRANSPOSE(I2)
        Endfor
     Endfor
  Endfor
  ; Clear sky I0 (NSZA,NVZA,NTOMS)
  grp2_id = H5G_OPEN(grp_id, 'Cloud Sky')
  datatype  = H5T_IDL_CREATE(I0_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(I0_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'I0', datatype, dataspace)
  H5D_WRITE, dat_id, I0_all
  unit   = '[W/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky I1 (NSZA,NVZA,NTOMS)
  datatype  = H5T_IDL_CREATE(I1_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(I1_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'I1', datatype, dataspace)
  H5D_WRITE, dat_id, I1_all
  unit   = '[W/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky I2 (NSZA,NVZA,NTOMS)
  datatype  = H5T_IDL_CREATE(I2_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(I2_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'I2', datatype, dataspace)
  H5D_WRITE, dat_id, I2_all
  unit   = '[W/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky Ir (NSZA,NVZA,NTOMS)
  datatype  = H5T_IDL_CREATE(Ir_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(Ir_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'Ir', datatype, dataspace)
  H5D_WRITE, dat_id, Ir_all
  unit   = '[W/cm^2]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky Sb (NTOMS)
  datatype  = H5T_IDL_CREATE(Sb_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(Sb_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'Sb', datatype, dataspace)
  H5D_WRITE, dat_id, Sb_all
  unit   = '[Unitless]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  H5G_CLOSE,grp2_id
  H5G_CLOSE,grp_id

  ; Fill in Scattering Weights group
  grp_id = H5G_OPEN(file_id,'Scattering Weights')
  ; Create clear sky group
  grp2_id = H5G_CREATE(grp_id, 'Clear Sky')
  H5G_CLOSE, grp2_id
  ; Create cloudy sky group
  grp2_id = H5G_CREATE(grp_id, 'Cloud Sky')
  H5G_CLOSE, grp2_id

  ; Create variables to hold data (Clear sky)
  dI0_all = FLTARR(NTOMS,NSURF,NLAY,NALBEDO,NVZA,NSZA)
  dI1_all = FLTARR(NTOMS,NSURF,NLAY,NALBEDO,NVZA,NSZA)
  dI2_all = FLTARR(NTOMS,NSURF,NLAY,NALBEDO,NVZA,NSZA)
  ; Fill in cloud sky variables
  For itoms = 0, ntoms-1 do begin
     For isurf = 0, nsurf-1 do begin
        filename = '../fitting_results_clear/Scattering_fit_'+$
                   TOMS[itoms]+'_'+SURFS[isurf]+'_clear.dat'
        print, filename
        RESTORE, filename
        dI0_all[itoms,isurf,*,*,*,*] = TRANSPOSE(dI0)
        dI1_all[itoms,isurf,*,*,*,*] = TRANSPOSE(dI1)
        dI2_all[itoms,isurf,*,*,*,*] = TRANSPOSE(dI2)
     Endfor
  Endfor
  ; Clear sky dI0 (NSZA,NVZA,NALBEDO,NLAY,NTOMS)
  grp2_id = H5G_OPEN(grp_id, 'Clear Sky')
  datatype  = H5T_IDL_CREATE(dI0_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(dI0_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'dI0', datatype, dataspace)
  H5D_WRITE, dat_id, dI0_all
  unit   = '[Unitless]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky dI1 (NSZA,NVZA,NALBEDO,NLAY,NTOMS)
  datatype  = H5T_IDL_CREATE(dI1_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(dI1_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'dI1', datatype, dataspace)
  H5D_WRITE, dat_id, dI1_all
  unit   = '[Unitless]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; Clear sky I2 (NSZA,NVZA,NALBEDO,NLAY,NTOMS)
  datatype  = H5T_IDL_CREATE(dI2_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(dI2_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'dI2', datatype, dataspace)
  H5D_WRITE, dat_id, dI2_all
  unit   = '[Unitless]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  H5G_CLOSE, grp2_id

  ; Create variables to hold data (cloud sky)
  dI0_all = FLTARR(NTOMS,NSURF,NLAY,NCLDP,NVZA,NSZA)
  dI1_all = FLTARR(NTOMS,NSURF,NLAY,NCLDP,NVZA,NSZA)
  dI2_all = FLTARR(NTOMS,NSURF,NLAY,NCLDP,NVZA,NSZA)
  ; Fill in cloud sky variables
  For itoms = 0, ntoms-1 do begin
     For isurf = 0, nsurf-1 do begin
        For ipre = 0, ncldp-1 do begin
           IF ( CLDP[ipre] GT SURF[isurf] ) THEN CONTINUE
           filename = '../fitting_results_cloud/Scattering_fit_'+$
                      TOMS[itoms]+'_'+SURFS[isurf]+'_cloud_'+CLDPS[ipre]+'.dat'
           print, filename
           RESTORE, filename
           dI0_all[itoms,isurf,*,ipre,*,*] = TRANSPOSE(dI0)
           dI1_all[itoms,isurf,*,ipre,*,*] = TRANSPOSE(dI1)
           dI2_all[itoms,isurf,*,ipre,*,*] = TRANSPOSE(dI2)
        Endfor
     Endfor
  Endfor
  ; cloud sky dI0 (NSZA,NVZA,NCLDP,NLAY,NTOMS)
  grp2_id = H5G_OPEN(grp_id, 'Cloud Sky')
  datatype  = H5T_IDL_CREATE(dI0_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(dI0_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'dI0', datatype, dataspace)
  H5D_WRITE, dat_id, dI0_all
  unit   = '[Unitless]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; cloud sky dI1 (NSZA,NVZA,NCLDP,NLAY,NTOMS)
  datatype  = H5T_IDL_CREATE(dI1_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(dI1_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'dI1', datatype, dataspace)
  H5D_WRITE, dat_id, dI1_all
  unit   = '[Unitless]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  ; cloud sky I2 (NSZA,NVZA,NCLDP,NLAY,NTOMS)
  datatype  = H5T_IDL_CREATE(dI2_all)
  dataspace = H5S_CREATE_SIMPLE(SIZE(dI2_all,/dimensions))
  dat_id    = H5D_CREATE(grp2_id, 'dI2', datatype, dataspace)
  H5D_WRITE, dat_id, dI2_all
  unit   = '[Unitless]' & datatype = H5T_IDL_CREATE(unit) & dataspace = H5S_CREATE_SIMPLE(1)
  att_id = H5A_CREATE(dat_id, 'Unit', datatype, dataspace)
  H5A_WRITE, att_id, unit
  H5A_CLOSE, att_id
  H5D_CLOSE, dat_id
  H5G_CLOSE, grp2_id
  H5G_CLOSE,grp_id

  ; Close file
  H5F_CLOSE, file_id

END
