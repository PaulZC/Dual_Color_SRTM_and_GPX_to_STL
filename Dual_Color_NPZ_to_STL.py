## Code for the Dual Color NPZ to STL Converter

## Converts a npz file into binary stl for dual color 3D printing

## The code is identical to SRTM_to_STL_1 \ NPZ_to_STL.py except it also creates a file called
## STL_Params.npz containing the name, position and scale of the SRTM STL. Dual_Color_GPX_to_STL.py
## uses that file to constrain the GPX route to the same position and scale

## Run from the command line as:
## python Dual_Color_NPZ_to_STL.py filename.npz printheight printwidth baseoffset basethickness
## printheight is the maximum print height in mm
## printwidth is the (maximum) print width (or height - whichever is larger) in mm
## baseoffset is the extra material in mm that is added to non-minimum data to help define coastlines etc. (set to zero as required)
## basethickness is the base layer thickness in mm (set to zero as required)

import numpy
import sys

if __name__ == '__main__':
    try:
        filename = ''
        print_height = 0.
        print_width = 0.
        base_offset = 0.
        base_thickness = 0.
        # Binary STL header must be 80 characters wide
        header = 'Generated from NASA SRTM data                                                   '
        params_file = 'STL_Params.npz'
        
        ## Everest
        ## N27E086_join_trim is roughly 9km square so scale to 90mm (1:100,000)
        #filename = 'N27E086_join_trim.npz'
        #print_height = 87.5 # 8748m / 100,000
        #print_width = 90.0
        #base_offset = 0.0
        #base_thickness = 0.0
        #header = 'Everest from NASA SRTM. XYZ Scale 1:~100,000.                                   '

        if filename == '':
            # Check if the npz filename was passed in argv
            if len(sys.argv) > 1: filename = sys.argv[1]
            # Check if the print height was passed in argv
            if len(sys.argv) > 2: print_height = float(sys.argv[2])
            # Check if the print width was passed in argv
            if len(sys.argv) > 3: print_width = float(sys.argv[3])
            # Check if the base offset was passed in argv
            if len(sys.argv) > 4: base_offset = float(sys.argv[4])
            # Check if the base thickness was passed in argv
            if len(sys.argv) > 5: base_thickness= float(sys.argv[5])

        if filename == '': filename = raw_input('Enter the npz filename: ') # Get the npz filename

        print 'Processing',filename

        outfile = filename[:-4] + '.stl' # Create the output file name
        
        try:
            # read data: width, height, max_hgt, east[], north[], hgt[]
            filedata = numpy.load(filename,mmap_mode='r')
            width = int(filedata['width'])
            height = int(filedata['height'])
            hgt_max = int(filedata['hgt_max'])
            east = numpy.ravel(filedata['east'])
            north = numpy.ravel(filedata['north'])
            hgt = numpy.ravel(filedata['hgt'])
            points = len(hgt) # Get the number of data points
            filedata.close()
        except:
            raise Exception('Invalid file!')

        if points != width * height:
            raise Exception('Invalid file size!')

        # Find the minimum height, ignoring invalid data points (-32768)
        hgt_min = 32767.
        for l in range(points):
            if (hgt[l] < hgt_min) and (hgt[l] > -32768.): hgt_min = hgt[l]

        hgt_max = hgt.max() # Find the max_height

        print 'Points:',points
        print 'Maximum height (before scaling):',hgt_max
        print 'Minimum height (before scaling):',hgt_min

        if print_height == 0.:
            try:
                print_height = float(raw_input('Enter the print height in mm: '))
            except:
                raise Exception('Invalid print height!')
            
        if print_width == 0.:
            try:
                print_width = float(raw_input('Enter the print width in mm: '))
            except:
                raise Exception('Invalid print width!')
            
        if base_offset == 0.:
            try:
                base_offset = float(raw_input('Enter the base offset (for coastline definition) in mm: '))
            except:
                raise Exception('Invalid base offset!')
            
        if base_thickness == 0.:
            try:
                base_thickness = float(raw_input('Enter the base layer thickness in mm: '))
            except:
                raise Exception('Invalid base thickness!')
            
        #Make all valid points zero height or greater
        if hgt_min < 0.:
            for l in range(points):
                if hgt[l] > -32768.: hgt[l] -= hgt_min

        #Scale height and add base_offset to all valid points to help define coastlines etc..
        #Zero invalid points
        for l in range(points):
            if hgt[l] > 0.:
                hgt[l] = base_offset + (print_height * hgt[l] / hgt_max)
            else:
                hgt[l] = 0.

        # Add base thickness
        hgt += base_thickness

        # Find extremes and save for use in Dual_Color_GPX_to_STL
        print 'Use the following values in Dual_Color_GPX_to_STL:'
        min_east = east.min()
        print 'Minimum East =',min_east
        max_east = east.max()
        print 'Maximum East =',max_east
        min_north = north.min()
        print 'Minimum North =',min_north
        max_north = north.max()
        print 'Maximum North =',max_north

        # Shift E,N origin to zero
        east -= min_east
        north -= min_north

        # Find the maximum E or N dimension
        max_east = east.max()
        max_north = north.max()
        maxd = max_north
        if max_east > max_north: maxd = max_east

        scale = print_width / maxd
        print 'Scale Factor =',scale

        print_max = hgt.max()
        print 'Maximum Print Height =',print_max

        print 'Saving to',params_file
        try:
            numpy.savez(params_file, stl_file=outfile, min_east=min_east, max_east=max_east, min_north=min_north, max_north=max_north, scale=scale, max_height = print_max)
        except:
            raise Exception('Could not save parameters!')

        # Scale to print width
        east *= scale
        north *= scale

        base = hgt * 0. # Create a zero base layer array for triangle generation

        # Duplicate coordinates arrays for triangle generation
        east = numpy.concatenate((east,east))
        north = numpy.concatenate((north,north))
        hgt = numpy.concatenate((hgt,base))

        print 'Generating triangles...'

        triangles = []
        #shades = [] # For VisCAM / SolidView color STL files : not currently implemented

        # Surface

        for row in range(height - 1):
            for col in range(width - 1):
                origin = (row * width) + col
                triangle = [origin,origin+width,origin+1]
                triangles.append(triangle)
                #shades.append((hgt[origin]))
                triangle = [origin+1,origin+width,origin+width+1]
                triangles.append(triangle)
                #shades.append((hgt[origin]))

        #Base

        for row in range(height - 1):
            for col in range(width - 1):
                origin = (row * width) + col + points
                triangle = [origin,origin+1,origin+width]
                triangles.append(triangle)
                #shades.append((hgt[origin]))
                triangle = [origin+1,origin+width+1,origin+width]
                triangles.append(triangle)
                #shades.append((hgt[origin]))

        #Walls

        for col in range(width - 1):
            origin = col
            triangle = [origin,origin+1,origin+points]
            triangles.append(triangle)
            #shades.append((hgt[origin]))
            triangle = [origin+1,origin+points+1,origin+points]
            triangles.append(triangle)
            #shades.append((hgt[origin]))

        for col in range(width - 1):
            origin = ((height - 1) * width) + col
            triangle = [origin,origin+points,origin+1]
            triangles.append(triangle)
            #shades.append((hgt[origin]))
            triangle = [origin+1,origin+points,origin+points+1]
            triangles.append(triangle)
            #shades.append((hgt[origin]))

        for row in range(height - 1):
            origin = row * width
            triangle = [origin,origin+points,origin+width]
            triangles.append(triangle)
            #shades.append((hgt[origin]))
            triangle = [origin+width,origin+points,origin+width+points]
            triangles.append(triangle)
            #shades.append((hgt[origin]))

        for row in range(height - 1):
            origin = (row * width) + (width - 1)
            triangle = [origin,origin+width,origin+points]
            triangles.append(triangle)
            #shades.append((hgt[origin]))
            triangle = [origin+width,origin+width+points,origin+points]
            triangles.append(triangle)
            #shades.append((hgt[origin]))

        print 'Created',len(triangles),'triangles'
        print 'File size should be',80+4+(50*len(triangles)),'bytes'

        print 'Writing STL file...'

        try: # Write data to file
            fp = open(outfile,'wb')
            fp.write(header)
            fp.write(numpy.uint32(len(triangles)))
            for l in range(len(triangles)):
                fp.write(numpy.float32(0.0))
                fp.write(numpy.float32(0.0))
                fp.write(numpy.float32(0.0))
                fp.write(numpy.float32(east[triangles[l][0]]))
                fp.write(numpy.float32(north[triangles[l][0]]))
                fp.write(numpy.float32(hgt[triangles[l][0]]))
                fp.write(numpy.float32(east[triangles[l][1]]))
                fp.write(numpy.float32(north[triangles[l][1]]))
                fp.write(numpy.float32(hgt[triangles[l][1]]))
                fp.write(numpy.float32(east[triangles[l][2]]))
                fp.write(numpy.float32(north[triangles[l][2]]))
                fp.write(numpy.float32(hgt[triangles[l][2]]))
                fp.write(numpy.uint16(0)) # <- data from shades[l] would go here
            fp.close()
        except:
            raise Exception('Could not write data!')

        print 'Complete!'

    except KeyboardInterrupt:
        print 'CTRL+C received...'
     
    finally:
        print 'Bye!'

