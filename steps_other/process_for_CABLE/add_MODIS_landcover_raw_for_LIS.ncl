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
NCAR_path   = "/g/data/ua8/LIS/LIS_PARAMS/UMD/1KM/landcover_CABLE_NCAR.1gd4r"
NCAR_LC     = fbindirread (NCAR_path, 0, (/15000,36000/), "float")
NCAR_LC_int = NCAR_LC ;tointeger(NCAR_LC)
; NCAR_LC_int@_FillValue = 16

NCAR_LC_int!0   = "lat"
NCAR_LC_int!1   = "lon"
NCAR_nlat       = floattointeger((89.995+59.995)/0.01+1)
NCAR_nlon       = floattointeger((179.995+179.995)/0.01+1)
NCAR_LC_int&lat = fspan(-59.995,89.995,  NCAR_nlat)
NCAR_LC_int&lon = fspan(-179.995,179.995, NCAR_nlon)
printVarSummary(NCAR_LC_int)


; ================================ read MODIS_LC ===================================
MODIS_in = "/g/data/w97/mm3972/data/MODIS/MODIS_landcover/var_landcover_only/MCD12Q1.061_landcover_clim_PFT_5km.nc"
f_lc     = addfile(MODIS_in, "r")
iveg     = f_lc->LC(:,:)*1.

i_num    = 0
do i = -44.,-10., 0.05
  j_num  = 0
  do j = 112., 154., 0.05
    printVarSummary(NCAR_LC_int)
    NCAR_LC_int({i-0.05:i+0.05},{j-0.05:j+0.05}) = iveg(i_num,j_num)
    j_num = j_num + 1
  end do
  i_num = i_num + 1 
end do

NCAR_LC_int@_FillValue = 16.

setfileoption("bin","WriteByteOrder","BigEndian")
file_out = "./binary_file/NCAR_AU_MODIS_land_cover_dataset.1gd4r"
system("/bin/rm -f " + file_out)
fbindirwrite(file_out,NCAR_LC_int)
print(max(NCAR_LC_int))
print(min(NCAR_LC_int))

; __________________ check output __________________
file_out_nc = "./binary_file/NCAR_AU_MODIS_land_cover_dataset.nc"
system("/bin/rm -f " + file_out_nc)    ; remove if exists
fout  = addfile (file_out_nc, "c")     ; open output file

setfileoption(fout,"DefineMode",True)

; create global attributes of the file
fAtt               = True            ; assign file attributes
fAtt@title         = "NCAR + MODIS dataset (over AU) landcover map"
fAtt@source_file   = "/g/data/ua8/LIS/LIS_PARAMS/UMD/1KM/landcover_CABLE_NCAR.1gd4r and /g/data/w97/mm3972/data/MODIS/MODIS_landcover/var_landcover_only/MCD12Q1.061_landcover_clim_PFT_5km.nc"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( fout, fAtt )            ; copy file attributes

dimNames = (/"lat", "lon"/)
dimSizes = (/ NCAR_nlat,  NCAR_nlon/)
dimUnlim = (/False, False/)
filedimdef( fout,dimNames,dimSizes,dimUnlim)

filevardef( fout, "lat", "float", "lat")
filevardef( fout, "lon", "float", "lon")
filevardef( fout, "Landcover", "float", (/ "lat", "lon" /))
setfileoption(fout,"DefineMode",False)

fout->lat@long_name     = "latitude"
fout->lat@standard_name = "latitude"
fout->lat@axis          = "Y"
fout->lat@units         = "degrees_North"
fout->lat               = NCAR_LC_int&lat*1.

fout->lon@long_name     = "longitude"
fout->lon@standard_name = "longitude"
fout->lon@axis          = "X"
fout->lon@units         = "degrees_East"
fout->lon               = NCAR_LC_int&lon*1.

fout->Landcover@long_name     = "land cover type"
fout->Landcover@units         = "-"
fout->Landcover@_FillValue    = 16
fout->Landcover               = NCAR_LC_int*1.

; ============================== Check output file ===============================

setfileoption("bin","ReadByteOrder","BigEndian")
NCAR_path = "./binary_file/NCAR_AU_MODIS_land_cover_dataset.1gd4r"
NCAR_LC   = fbindirread (NCAR_path, 0, (/14902,19161/), "float") ;
printVarSummary(NCAR_LC)
plot      = gsn_csm_contour(wks,NCAR_LC,res)

end
