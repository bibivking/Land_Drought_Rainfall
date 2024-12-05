import pygrib

# Open the GRIB file
file_path = "/g/data/ua8/LIS/LIS_PARAMS/UMD/1KM/elev_GTOPO30.1gd4r"
grbs = pygrib.open(file_path)

# Read the first message (or iterate over all messages if needed)
grb = grbs.message(1)  # Access the first GRIB message
print(grb)

# Extract latitude and longitude grid
lats, lons = grb.latlons()

# Print information
print("Latitudes:", lats)
print("Longitudes:", lons)

# Close the GRIB file
grbs.close()
