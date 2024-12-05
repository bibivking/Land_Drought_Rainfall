;*******************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;*******************************************

; __________________________________ MMY _______________________________________
; Purpose: Regrid Offline CABLE output to LIS-CABLE resolution
; ______________________________________________________________________________


;______ From Mark Decker's code ______
setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 5000000000 ;
end setvalues
;_____________________________________

begin

; ================================== Options ===================================
fname_lis_in   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght1719_bdy_data/bdy_data/lis_input.d01.nc"
fname_lis_out  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght1719_bdy_data/bdy_data/lis_input_new.d01.nc"
fname_soil_in  = "/g/data/w97/mm3972/scripts/wrf_scripts/make_LIS_landinfo/nc_file/Openlandmap_soilcomposition_CORDEX_180E_depth_varying.nc"

system("rm " + fname_lis_out)
system("cp " + fname_lis_in + " " + fname_lis_out)

var_names = (/"SAND_VEC","CLAY_VEC","SILT_VEC","OC_VEC","BULK_DEN_VEC"/)

; _________________ set range of lat and lon ______________________
fout      = addfile(fname_lis_out, "w")

setfileoption(fout,"DefineMode",True)

; create global attributes of the file
fAtt               = True           ; assign file attributes
fAtt@modification  = "Add SAND_VEC, CLAY_VEC, SILT_VEC, OCSOIL_VEC & BULKDENSITY_VEC from "+fname_soil_in
fileattdef( fout, fAtt )            ; copy file attributes

print("point 1")

dimNames = (/"soil_depth"/)
dimSizes = (/ 6/)
dimUnlim = (/False/)
filedimdef( fout,dimNames,dimSizes,dimUnlim)

filevardef( fout, "SAND_VEC", "float", (/"soil_depth","north_south", "east_west"/) )
filevardef( fout, "CLAY_VEC", "float", (/"soil_depth","north_south", "east_west"/) )
filevardef( fout, "SILT_VEC", "float", (/"soil_depth","north_south", "east_west"/) )
filevardef( fout, "OCSOIL_VEC", "float", (/"soil_depth","north_south", "east_west"/) )
filevardef( fout, "BULKDENSITY_VEC", "float", (/"soil_depth","north_south", "east_west"/) )

setfileoption(fout,"DefineMode",False)
print("point 2")

fout->SAND_VEC@standard_name = "Sand fraction of each soil layer"
fout->SAND_VEC@units         = ""
fout->SAND_VEC@scale_factor  = 1.
fout->SAND_VEC@add_offset    = 0.
fout->SAND_VEC@missing_value = -9999.
fout->SAND_VEC@vmin          = 0.
fout->SAND_VEC@vmax          = 0.
fout->SAND_VEC@num_bins      = -1

fout->CLAY_VEC@standard_name = "Clay fraction of each soil layer"
fout->CLAY_VEC@units         = ""
fout->CLAY_VEC@scale_factor  = 1.
fout->CLAY_VEC@add_offset    = 0.
fout->CLAY_VEC@missing_value = -9999.
fout->CLAY_VEC@vmin          = 0.
fout->CLAY_VEC@vmax          = 0.
fout->CLAY_VEC@num_bins      = -1

fout->SILT_VEC@standard_name = "Silt fraction of each soil layer"
fout->SILT_VEC@units         = ""
fout->SILT_VEC@scale_factor  = 1.
fout->SILT_VEC@add_offset    = 0.
fout->SILT_VEC@missing_value = -9999.
fout->SILT_VEC@vmin          = 0.
fout->SILT_VEC@vmax          = 0.
fout->SILT_VEC@num_bins      = -1

fout->OCSOIL_VEC@standard_name = "Organic Carbon soil content of each soil layer"
fout->OCSOIL_VEC@units         = "-"
fout->OCSOIL_VEC@scale_factor  = 1.
fout->OCSOIL_VEC@add_offset    = 0.
fout->OCSOIL_VEC@missing_value = -9999.
fout->OCSOIL_VEC@vmin          = 0.
fout->OCSOIL_VEC@vmax          = 0.
fout->OCSOIL_VEC@num_bins      = -1

fout->BULKDENSITY_VEC@standard_name = "Soil bulk density of each soil layer"
fout->BULKDENSITY_VEC@units         = ""
fout->BULKDENSITY_VEC@scale_factor  = 1.
fout->BULKDENSITY_VEC@add_offset    = 0.
fout->BULKDENSITY_VEC@missing_value = -9999.
fout->BULKDENSITY_VEC@vmin          = 0.
fout->BULKDENSITY_VEC@vmax          = 0.
fout->BULKDENSITY_VEC@num_bins      = -1
print("point 3")

; ===================== Interpolation =====================
fsoil     = addfile(fname_soil_in, "r")
lat_in    = fsoil->lat(1500:3300)
lon_in    = fsoil->lon(4500:7300)

print(lat_in)
print(lon_in)

nlat      = dimsizes(lat_in)
nlon      = dimsizes(lon_in)
lat_in_2D = new((/nlat,nlon/), float)
lon_in_2D = new((/nlat,nlon/), float)

do j = 0, nlon-1
  lat_in_2D(:,j) = lat_in
end do

do i = 0, nlat-1
  lon_in_2D(i,:) = lon_in
end do

lat_in_1D  = ndtooned(lat_in_2D)
lon_in_1D  = ndtooned(lon_in_2D)

sand_vec   = fsoil->SAND_VEC(:,1500:3300,4500:7300)
clay_vec   = fsoil->CLAY_VEC(:,1500:3300,4500:7300)
silt_vec   = fsoil->SILT_VEC(:,1500:3300,4500:7300)
ocsoil_vec = fsoil->OC_VEC(:,1500:3300,4500:7300)
rhosoil_vec= fsoil->BULK_DEN_VEC(:,1500:3300,4500:7300)
delete(fsoil)
print("point 4")

; ----------- get lat/lon out info ------------
lat_out_2D = fout->lat
lon_out_2D = fout->lon
nlat_out   = dimsizes(lat_out_2D(:,0))
nlon_out   = dimsizes(lat_out_2D(0,:))
SAND_VEC   = new((/6,nlat_out,nlon_out/),float)
CLAY_VEC   = new((/6,nlat_out,nlon_out/),float)
SILT_VEC   = new((/6,nlat_out,nlon_out/),float)
OCSOIL_VEC = new((/6,nlat_out,nlon_out/),float)
BULKDENSITY_VEC = new((/6,nlat_out,nlon_out/),float)

print("point 5")

do lvl = 0,5
  print("lvl="+lvl)
	sand_tmp     = ndtooned(sand_vec(lvl,:,:))
	clay_tmp     = ndtooned(clay_vec(lvl,:,:))
	silt_tmp     = ndtooned(silt_vec(lvl,:,:))
	ocsoil_tmp   = ndtooned(ocsoil_vec(lvl,:,:))
	rhosoil_tmp  = ndtooned(rhosoil_vec(lvl,:,:))

 	sand     = triple2grid2d(lon_in_1D,lat_in_1D,sand_tmp,lon_out_2D,lat_out_2D, False)
	clay     = triple2grid2d(lon_in_1D,lat_in_1D,clay_tmp,lon_out_2D,lat_out_2D, False)
	silt     = triple2grid2d(lon_in_1D,lat_in_1D,silt_tmp,lon_out_2D,lat_out_2D, False)
	ocsoil   = triple2grid2d(lon_in_1D,lat_in_1D,ocsoil_tmp,lon_out_2D,lat_out_2D, False)
	rhosoil  = triple2grid2d(lon_in_1D,lat_in_1D,rhosoil_tmp,lon_out_2D,lat_out_2D, False)

	SAND_VEC(lvl,:,:)        = sand
	CLAY_VEC(lvl,:,:)        = clay
	SILT_VEC(lvl,:,:)        = silt
	OCSOIL_VEC(lvl,:,:)      = ocsoil
	BULKDENSITY_VEC(lvl,:,:) = rhosoil
end do

fout->SAND_VEC        = SAND_VEC
fout->CLAY_VEC        = CLAY_VEC
fout->SILT_VEC        = SILT_VEC
fout->OCSOIL_VEC      = OCSOIL_VEC
fout->BULKDENSITY_VEC = BULKDENSITY_VEC

end
