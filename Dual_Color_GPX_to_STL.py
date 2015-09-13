## Code for the Dual Color GPX Converter

## Run from the command line as:
## python Dual_Color_GPX_to_STL.py gpx_filename.gpx uk track_width force_UTM_zone

## Step 1 : Convert the GPX file into a DXF (containing lines and arcs)
## Step 2 : Write an OpenSCAD script to:
##          Offset and Linear_Extrude the DXF into 3D form
##          Export the Difference between the DXF and the STL landscape as STL (Color 1)
## Step 3 : Write a second OpenSCAD script to:
##          Export the Intersection between the DXF and the STL landscape as STL (color 2)
## Step 4 : Run both scripts through OpenSCAD
## You can then use Slic3r to print the two STL files in dual color

## Download dxfwrite from: https://pypi.python.org/pypi/dxfwrite/
## Documentation at: https://pythonhosted.org/dxfwrite/

## utm gratefully plagiarised from: https://pypi.python.org/pypi/utm

## WGS84toOSGB36 gratefully plagiarised from:
## http://hannahfry.co.uk/2012/02/01/converting-latitude-and-longitude-to-british-national-grid/

import numpy
import sys
from scipy import *
import string
from dxfwrite import DXFEngine as dxf
import os

def WGS84toOSGB36(lat, lon):

    #First convert to radians
    #These are on the wrong ellipsoid currently: GRS80. (Denoted by _1)
    lat_1 = lat*pi/180
    lon_1 = lon*pi/180

    #Want to convert to the Airy 1830 ellipsoid, which has the following:
    a_1, b_1 =6378137.000, 6356752.3141 #The GSR80 semi-major and semi-minor axes used for WGS84(m)
    e2_1 = 1- (b_1*b_1)/(a_1*a_1)   #The eccentricity of the GRS80 ellipsoid
    nu_1 = a_1/sqrt(1-e2_1*sin(lat_1)**2)

    #First convert to cartesian from spherical polar coordinates
    H = 0 #Third spherical coord.
    x_1 = (nu_1 + H)*cos(lat_1)*cos(lon_1)
    y_1 = (nu_1+ H)*cos(lat_1)*sin(lon_1)
    z_1 = ((1-e2_1)*nu_1 +H)*sin(lat_1)

    #Perform Helmut transform (to go between GRS80 (_1) and Airy 1830 (_2))
    s = 20.4894*10**-6 #The scale factor -1
    tx, ty, tz = -446.448, 125.157, -542.060 #The translations along x,y,z axes respectively
    rxs,rys,rzs = -0.1502, -0.2470, -0.8421  #The rotations along x,y,z respectively, in seconds
    rx, ry, rz = rxs*pi/(180*3600.), rys*pi/(180*3600.), rzs*pi/(180*3600.) #In radians
    x_2 = tx + (1+s)*x_1 + (-rz)*y_1 + (ry)*z_1
    y_2 = ty + (rz)*x_1  + (1+s)*y_1 + (-rx)*z_1
    z_2 = tz + (-ry)*x_1 + (rx)*y_1 +  (1+s)*z_1

    #Back to spherical polar coordinates from cartesian
    #Need some of the characteristics of the new ellipsoid
    a, b = 6377563.396, 6356256.909 #The GSR80 semi-major and semi-minor axes used for WGS84(m)
    e2 = 1- (b*b)/(a*a)   #The eccentricity of the Airy 1830 ellipsoid
    p = sqrt(x_2**2 + y_2**2)

    #Lat is obtained by an iterative proceedure:
    lat = arctan2(z_2,(p*(1-e2))) #Initial value
    latold = 2*pi
    while abs(lat - latold)>10**-16:
        lat, latold = latold, lat
        nu = a/sqrt(1-e2*sin(latold)**2)
        lat = arctan2(z_2+e2*nu*sin(latold), p)

    #Lon and height are then pretty easy
    lon = arctan2(y_2,x_2)
    H = p/cos(lat) - nu

    #E, N are the British national grid coordinates - eastings and northings
    F0 = 0.9996012717                   #scale factor on the central meridian
    lat0 = 49*pi/180                    #Latitude of true origin (radians)
    lon0 = -2*pi/180                    #Longtitude of true origin and central meridian (radians)
    N0, E0 = -100000, 400000            #Northing & easting of true origin (m)
    n = (a-b)/(a+b)

    #meridional radius of curvature
    rho = a*F0*(1-e2)*(1-e2*sin(lat)**2)**(-1.5)
    eta2 = nu*F0/rho-1

    M1 = (1 + n + (5/4)*n**2 + (5/4)*n**3) * (lat-lat0)
    M2 = (3*n + 3*n**2 + (21/8)*n**3) * sin(lat-lat0) * cos(lat+lat0)
    M3 = ((15/8)*n**2 + (15/8)*n**3) * sin(2*(lat-lat0)) * cos(2*(lat+lat0))
    M4 = (35/24)*n**3 * sin(3*(lat-lat0)) * cos(3*(lat+lat0))

    #meridional arc
    M = b * F0 * (M1 - M2 + M3 - M4)          

    I = M + N0
    II = nu*F0*sin(lat)*cos(lat)/2
    III = nu*F0*sin(lat)*cos(lat)**3*(5- tan(lat)**2 + 9*eta2)/24
    IIIA = nu*F0*sin(lat)*cos(lat)**5*(61- 58*tan(lat)**2 + tan(lat)**4)/720
    IV = nu*F0*cos(lat)
    V = nu*F0*cos(lat)**3*(nu/rho - tan(lat)**2)/6
    VI = nu*F0*cos(lat)**5*(5 - 18* tan(lat)**2 + tan(lat)**4 + 14*eta2 - 58*eta2*tan(lat)**2)/120

    N = I + II*(lon-lon0)**2 + III*(lon- lon0)**4 + IIIA*(lon-lon0)**6
    E = E0 + IV*(lon-lon0) + V*(lon- lon0)**3 + VI*(lon- lon0)**5 

    #Job's a good'n.
    return E,N

# UTM

K0 = 0.9996

E = 0.00669438
E2 = E * E
E3 = E2 * E
E_P2 = E / (1.0 - E)

SQRT_E = math.sqrt(1 - E)
_E = (1 - SQRT_E) / (1 + SQRT_E)
_E2 = _E * _E
_E3 = _E2 * _E
_E4 = _E3 * _E
_E5 = _E3 * _E

M1 = (1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256)
M2 = (3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024)
M3 = (15 * E2 / 256 + 45 * E3 / 1024)
M4 = (35 * E3 / 3072)

P2 = (3. / 2 * _E - 27. / 32 * _E3 + 269. / 512 * _E5)
P3 = (21. / 16 * _E2 - 55. / 32 * _E4)
P4 = (151. / 96 * _E3 - 417. / 128 * _E5)
P5 = (1097. / 512 * _E4)

R = 6378137

ZONE_LETTERS = [
    (84, None), (72, 'X'), (64, 'W'), (56, 'V'), (48, 'U'), (40, 'T'),
    (32, 'S'), (24, 'R'), (16, 'Q'), (8, 'P'), (0, 'N'), (-8, 'M'), (-16, 'L'),
    (-24, 'K'), (-32, 'J'), (-40, 'H'), (-48, 'G'), (-56, 'F'), (-64, 'E'),
    (-72, 'D'), (-80, 'C')
]

def from_latlon(latitude, longitude, force_zone):
    if not -80.0 <= latitude <= 84.0:
        raise OutOfRangeError('latitude out of range (must be between 80 deg S and 84 deg N)')
    if not -180.0 <= longitude <= 180.0:
        raise OutOfRangeError('northing out of range (must be between 180 deg W and 180 deg E)')

    lat_rad = math.radians(latitude)
    lat_sin = math.sin(lat_rad)
    lat_cos = math.cos(lat_rad)

    lat_tan = lat_sin / lat_cos
    lat_tan2 = lat_tan * lat_tan
    lat_tan4 = lat_tan2 * lat_tan2

    lon_rad = math.radians(longitude)

    zone_number = latlon_to_zone_number(latitude, longitude)
    if force_zone != '': zone_number = force_zone
    central_lon = zone_number_to_central_longitude(zone_number)
    central_lon_rad = math.radians(central_lon)

    zone_letter = latitude_to_zone_letter(latitude)

    n = R / math.sqrt(1 - E * lat_sin**2)
    c = E_P2 * lat_cos**2

    a = lat_cos * (lon_rad - central_lon_rad)
    a2 = a * a
    a3 = a2 * a
    a4 = a3 * a
    a5 = a4 * a
    a6 = a5 * a

    m = R * (M1 * lat_rad -
             M2 * math.sin(2 * lat_rad) +
             M3 * math.sin(4 * lat_rad) -
             M4 * math.sin(6 * lat_rad))

    easting = K0 * n * (a +
                        a3 / 6 * (1 - lat_tan2 + c) +
                        a5 / 120 * (5 - 18 * lat_tan2 + lat_tan4 + 72 * c - 58 * E_P2)) + 500000

    northing = K0 * (m + n * lat_tan * (a2 / 2 +
                                        a4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c**2) +
                                        a6 / 720 * (61 - 58 * lat_tan2 + lat_tan4 + 600 * c - 330 * E_P2)))

    if latitude < 0:
        northing += 10000000

    return easting, northing, zone_number, zone_letter


def latitude_to_zone_letter(latitude):
    for lat_min, zone_letter in ZONE_LETTERS:
        if latitude >= lat_min:
            return zone_letter

    return None


def latlon_to_zone_number(latitude, longitude):
    if 56 <= latitude <= 64 and 3 <= longitude <= 12:
        return 32

    if 72 <= latitude <= 84 and longitude >= 0:
        if longitude <= 9:
            return 31
        elif longitude <= 21:
            return 33
        elif longitude <= 33:
            return 35
        elif longitude <= 42:
            return 37

    return int((longitude + 180) / 6) + 1


def zone_number_to_central_longitude(zone_number):
    return (zone_number - 1) * 6 - 180 + 3

if __name__ == '__main__':
    try:
        print 'Dual Color GPX'

        dtw = 2.0 # Default track width (mm)
        dth = 200. # Default track height (mm)

        script1 = 'Dual_Color_GPX_1.scad'        
        script2 = 'Dual_Color_GPX_2.scad'
        openscad = "\"C:\\Program Files\\OpenSCAD\\openscad\""
        params_file = 'STL_Params.npz'

        gpx_filename = ''
        stl_filename = ''
        UK = ''
        force_zone = ''
        track_width = ''
        track_height = ''
        east_offset = 0.
        north_offset = 0.

        print 'Reading parameters...'
        try:
            filedata = numpy.load(params_file,mmap_mode='r')
            stl_filename = str(filedata['stl_file'])
            min_east = float(filedata['min_east'])
            max_east = float(filedata['max_east'])
            min_north = float(filedata['min_north'])
            max_north = float(filedata['max_north'])
            scale = float(filedata['scale'])
            track_height = float(filedata['max_height']) + 10. # Add an extra 10mm to track_height as a precaution
        except:
            raise Exception('Invalid parameters file!')
            
        ## Everest
        #gpx_filename = 'Everest_1953_Route.gpx'
        #UK = 'n'
        #track_width = dtw
        
        east_offset = min_east
        north_offset = min_north
        dp = 20. * scale # Minimum spacing (m) to avoid duplicate points
        
        if gpx_filename == '':
            # Check if the gpx filename was passed in argv
            if len(sys.argv) > 1: gpx_filename = sys.argv[1]
            # Check if the UK flag was passed in argv 
            if len(sys.argv) > 2: UK = sys.argv[2]
            # Check if track width was passed in argv
            if len(sys.argv) > 3: track_width = sys.argv[3]
            # Check if a forced UTM zone was passed in argv
            if len(sys.argv) > 4: force_zone = sys.argv[4]
            
        if gpx_filename == '': gpx_filename = raw_input('Enter the gpx filename: ') # Get the gpx filename

        if UK == '': UK = raw_input('Is this UK data (y/n)? ') # Is it UK data?
        if UK == '' or UK == 'Y' or UK == 'y': UK = True
        else: UK = False

        if UK == False: # Force a UTM zone for non-UK data?
            if force_zone == '': force_zone = raw_input('If you want to force a UTM zone, enter it now: ')
            if force_zone != '':
                try:
                    force_zone = int(force_zone)
                except:
                    raise Exception('Invalid zone!')

        if track_width == '': track_width = raw_input('Enter the track width in mm (default '+str(dtw)+'): ')
        if track_width != '':
            try:
                track_width = float(track_width)
            except:
                raise Exception('Invalid track width!')
        else:
            track_width = dtw
                
        if track_height == '': track_height = raw_input('Enter the track height in mm (default '+str(dth)+'): ')
        if track_height != '':
            try:
                track_height = float(track_height)
            except:
                raise Exception('Invalid track height!')
        else:
            track_height = dth
                
        print 'Processing',gpx_filename

        dxffile = gpx_filename[:-4] + '.dxf'
        
        try:
            fi = open(gpx_filename,"r")
        except:
            raise Exception('Invalid gpx file!')

        east = []
        north = []

        for line in fi:
            if ('<trkpt' in line) or ('<rtept' in line): # identify route or track points
                nums = []
                for num in string.split(line,sep='"'): # separate off the lat and lon
                    try:
                        nums.append(float(num))
                    except ValueError:
                        pass
                if len(nums) == 2: # check we have (only) lat and lon
                    if UK:
                        e,n = WGS84toOSGB36(nums[0],nums[1]) # UK
                    else:
                        e,n,zn,zl = from_latlon(nums[0],nums[1],force_zone) # Rest of the World
                    if (len(east)) > 0:
                        # Check for duplicate point by checking we have moved by at least dp meters
                        move = ((((e - east_offset) * scale) - east[-1])**2 + (((n - north_offset) * scale) - north[-1])**2)**0.5
                        if move > dp:
                        #if (e < east[-1] - dp) or (e > east[-1] + dp) or (n < north[-1] - dp) or (n > north[-1] + dp):
                            east.append((e - east_offset) * scale)
                            north.append((n - north_offset) * scale)
                    else:
                        east.append((e - east_offset) * scale)
                        north.append((n - north_offset) * scale)
        fi.close()

        track_points = len(east)

        print 'Found',track_points,'non-duplicate points'
        if track_points < 2: raise Exception('Not enough track points!')

        # Scale the data
        #east = east - min_east
        #north = north - min_north
        #east = east * scale
        #north = north * scale

        print 'Generating DXF data...'

        drawing = dxf.drawing(dxffile)

        for l in range(track_points-1):

            # Calculate bearing from point 0 to point 1
            b01 = (numpy.pi / 2.) - numpy.arctan2(north[l+1] - north[l], east[l+1] - east[l])

            # Create two points on a line parallel to the one running from point 0 to point 1
            x1 = east[l] + ((track_width / 4.) * numpy.cos(b01))
            y1 = north[l] - ((track_width / 4.) * numpy.sin(b01))
            x2 = east[l+1] + ((track_width / 4.) * numpy.cos(b01))
            y2 = north[l+1] - ((track_width / 4.) * numpy.sin(b01))

            # First line
            drawing.add(dxf.line((x1,y1),(x2,y2),color=1))
            # First arc
            drawing.add(dxf.arc((track_width / 4.),(east[l+1],north[l+1]),numpy.degrees(0.-b01),numpy.degrees(numpy.pi-b01),color=1))
            
            # Now do it again creating a line on the other side of the track
            x3 = east[l+1] - ((track_width / 4.) * numpy.cos(b01))
            y3 = north[l+1] + ((track_width / 4.) * numpy.sin(b01))
            x4 = east[l] - ((track_width / 4.) * numpy.cos(b01))
            y4 = north[l] + ((track_width / 4.) * numpy.sin(b01))

            # Second line
            drawing.add(dxf.line((x3,y3),(x4,y4),color=1))
            # Second arc
            drawing.add(dxf.arc((track_width / 4.),(east[l],north[l]),numpy.degrees(numpy.pi-b01),numpy.degrees(0.-b01),color=1))
            # Repeat first line (possibly redundant?)
            drawing.add(dxf.line((x1,y1),(x2,y2),color=1))

        # Write data to dxf file
        print 'Writing DXF file',dxffile
        try: 
            drawing.save()
        except:
            raise Exception('Could not write dxf file!')

        # Write script files
        print 'Creating OpenSCAD scripts',script1,'&',script2
        try:
            fo1 = open(script1,"w")
        except:
            raise Exception('Could not open script file 1!')
        try:
            fo2 = open(script2,"w")
        except:
            raise Exception('Could not open script file 2!')

        filestr = 'difference ()\r\n'
        fo1.write(filestr)
        filestr = 'intersection ()\r\n'
        fo2.write(filestr)
        filestr = '{\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\timport(\"'+stl_filename+'\",convexity=10);\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\tlinear_extrude(height='+str(track_height)+',convexity=10)\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t{\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\tunion()\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t{\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t\toffset(r='+str(track_width/4.)+')\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t\t{\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t\t\timport(\"'+dxffile+'\",convexity=10);\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t\t}\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t\toffset(r=-'+str(track_width/4.)+')\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t\t{\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t\t\timport(\"'+dxffile+'\",convexity=10);\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t\t}\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t\t}\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '\t}\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        filestr = '}\r\n'
        fo1.write(filestr)
        fo2.write(filestr)
        fo1.close()
        fo2.close()

        # Run script 1 through OpenSCAD
        stlfile = "-o "+gpx_filename[:-4]+"_Color_1.stl"
        script = openscad+" "+stlfile+" "+script1
        print 'Executing',script1
        print 'Result (should be 0):',os.system(script)

        # Run script 2 through OpenSCAD
        stlfile = "-o "+gpx_filename[:-4]+"_Color_2.stl"
        script = openscad+" "+stlfile+" "+script2
        print 'Executing',script2
        print 'Result (should be 0):',os.system(script)

        print 'Complete!'

    except KeyboardInterrupt:
        print 'CTRL+C received...'
     
    finally:
        fi.close()
        fo1.close()
        fo2.close()
        print 'Bye!'

