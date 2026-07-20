load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________


; ================================ read NCAR_LC ===================================
setfileoption("bin","ReadByteOrder","BigEndian")
NCAR_path   = "/g/data/w97/mm3972/scripts/wrf_scripts/make_LIS_landinfo/binary_files/NCAR_AU_Dynamic_land_cover_dataset.1gd4r"
; NCAR_path   = "/g/data/ua8/LIS/LIS_PARAMS/UMD/1KM/landcover_CABLE_NCAR.1gd4r"
; NCAR_path   = "./for_LIS_WRF/NCAR_AU_MODIS_land_cover_dataset.1gd4r"
NCAR_LC     = fbindirread (NCAR_path, 0, (/15000,36000/), "float")
NCAR_LC_int = NCAR_LC ;tointeger(NCAR_LC)
; NCAR_LC_int@_FillValue = 16
printVarSummary(NCAR_LC_int)
; NCAR_LC_int!0   = "lat"
; NCAR_LC_int!1   = "lon"
; NCAR_nlat       = floattointeger((89.995+59.995)/0.01+1)
; NCAR_nlon       = floattointeger((179.995+179.995)/0.01+1)
; NCAR_LC_int&lat = fspan(-59.995,89.995,  NCAR_nlat)
; NCAR_LC_int&lon = fspan(-179.995,179.995, NCAR_nlon)
; printVarSummary(NCAR_LC_int)

; NCAR_LC_int     = where(NCAR_LC_int .lt. 0, 0, NCAR_LC_int)
; NCAR_LC_int     = where(NCAR_LC_int .gt. 17, 0, NCAR_LC_int)


; ; ================================ write MODIS_LC binary ===================================
; setfileoption("bin","WriteByteOrder","BigEndian")
; file_out        = "./for_LIS_WRF/NCAR_AU_MODIS_land_cover_dataset_ocean_to_zero.1gd4r"
; system("/bin/rm -f " + file_out)
; fbindirwrite(file_out,NCAR_LC_int)
; print(max(NCAR_LC_int))
; print(min(NCAR_LC_int))

; ; __________________ check output __________________
; file_out_nc     = "./for_LIS_WRF/NCAR_AU_MODIS_land_cover_dataset_ocean_to_zero.nc"
; system("/bin/rm -f " + file_out_nc)    ; remove if exists
; fout            = addfile (file_out_nc, "c")     ; open output file

; setfileoption(fout,"DefineMode",True)

; ; create global attributes of the file
; fAtt               = True            ; assign file attributes
; fAtt@title         = "NCAR + MODIS dataset (over AU) landcover map"
; fAtt@source_file   = "/g/data/ua8/LIS/LIS_PARAMS/UMD/1KM/landcover_CABLE_NCAR.1gd4r and /g/data/w97/mm3972/data/MODIS/MODIS_landcover/var_landcover_only/MCD12Q1.061_landcover_clim_PFT_5km.nc"
; fAtt@Conventions   = "None"
; fAtt@creation_date = systemfunc ("date")
; fileattdef( fout, fAtt )            ; copy file attributes

; dimNames = (/"lat", "lon"/)
; dimSizes = (/ NCAR_nlat,  NCAR_nlon/)
; dimUnlim = (/False, False/)
; filedimdef( fout,dimNames,dimSizes,dimUnlim)

; filevardef( fout, "lat", "float", "lat")
; filevardef( fout, "lon", "float", "lon")
; filevardef( fout, "Landcover", "float", (/ "lat", "lon" /))
; setfileoption(fout,"DefineMode",False)

; fout->lat@long_name     = "latitude"
; fout->lat@standard_name = "latitude"
; fout->lat@axis          = "Y"
; fout->lat@units         = "degrees_North"
; fout->lat               = NCAR_LC_int&lat*1.

; fout->lon@long_name     = "longitude"
; fout->lon@standard_name = "longitude"
; fout->lon@axis          = "X"
; fout->lon@units         = "degrees_East"
; fout->lon               = NCAR_LC_int&lon*1.

; fout->Landcover@long_name     = "land cover type"
; fout->Landcover@units         = "-"
; fout->Landcover@_FillValue    = 0
; fout->Landcover               = NCAR_LC_int*1.

; ; ============================== Check output file ===============================

; setfileoption("bin","ReadByteOrder","BigEndian")
; NCAR_path = "./for_LIS_WRF/NCAR_AU_MODIS_land_cover_dataset_ocean_to_zero.1gd4r"
; NCAR_LC   = fbindirread (NCAR_path, 0, (/15000,36000/), "float") ;
; printVarSummary(NCAR_LC)
; plot      = gsn_csm_contour(wks,NCAR_LC,res)

end
