# For each cylinder in the scan, find its cartesian coordinates,
# in the scanner's coordinate system.
# Write the result to a file which contains all cylinders, for all scans.
# 03_d_find_cylinders_cartesian
# Claus Brenner, 09 NOV 2012
from lego_robot import *
from math import sin, cos

# Find the derivative in scan data, ignoring invalid measurements.
def compute_derivative(scan, min_dist):
    jumps = [ 0 ]
    for i in xrange(1, len(scan) - 1):
        l = scan[i-1]
        r = scan[i+1]
        if l > min_dist and r > min_dist:
            derivative = (r - l) / 2.0
            jumps.append(derivative)
        else:
            jumps.append(0)
    jumps.append(0)
    return jumps

# For each area between a left falling edge and a right rising edge,
# determine the average ray number and the average depth.
def find_cylinders(scan, scan_derivative, jump, min_dist):
    cylinder_list = []
    on_cylinder = False
    sum_ray, sum_depth, rays = 0.0, 0.0, 0
    average_ray = 0.0
    average_depth = 0.0
    for i in xrange(len(scan_derivative)):
        if abs(scan_derivative[i])<jump and on_cylinder == False:
            pass
        elif abs(scan_derivative[i])<jump and on_cylinder == True:
            rays += 1  
            sum_ray += i
            if scan[i]< min_dist:
                
                sum_depth += scan[i+1]
            else:
                sum_depth += scan[i]


        elif scan_derivative[i] <= jump*(-1):
            if on_cylinder == False:
                on_cylinder = True
                rays += 1 
                sum_ray += i
                if scan[i]< min_dist:                
                    sum_depth += scan[i+1]
                else:
                    sum_depth += scan[i]
            else:                
                rays = 1 
                sum_ray = i
                if scan[i]< min_dist:
                
                    sum_depth = scan[i+1]
                else:
                    sum_depth = scan[i]
        elif scan_derivative[i] >= jump and on_cylinder == False:
            pass
            
        elif scan_derivative[i] >= jump and on_cylinder == True:       
            average_ray = sum_ray/rays
            average_depth = sum_depth/rays
            cylinder_list.append( (average_ray, average_depth) )
            on_cylinder = False
            sum_ray, sum_depth, rays = 0.0, 0.0, 0
            average_ray = 0.0
            average_depth = 0.0
            
            
            
        # --->>> Insert your cylinder code here.
        # Whenever you find a cylinder, add a tuple
        # (average_ray, average_depth) to the cylinder_list.

        # Just for fun, I'll output some cylinders.
        # Replace this by your code.
##        if i % 100 == 0:
##            cylinder_list.append( (i, scan[i]) )
    # --->>> Insert here your previous solution from find_cylinders_question.py.

    return cylinder_list

def compute_cartesian_coordinates(cylinders, cylinder_offset):
    result = []
    angle = 0.0
    for c in cylinders:
        angle = LegoLogfile.beam_index_to_angle(c[0])
        # --->>> Insert here the conversion from polar to Cartesian coordinates.
        # c is a tuple (beam_index, range).
        # For converting the beam index to an angle, use
        # LegoLogfile.beam_index_to_angle(beam_index)
        result.append( (cos(angle)*(c[1]+cylinder_offset),sin(angle)*(c[1]+cylinder_offset)) ) # Replace this by your (x,y)
    return result
        

if __name__ == '__main__':

    minimum_valid_distance = 20.0
    depth_jump = 100.0
    cylinder_offset = 90.0

    # Read the logfile which contains all scans.
    logfile = LegoLogfile()
    logfile.read("robot4_scan.txt")

    # Write a result file containing all cylinder records.
    # Format is: D C x[in mm] y[in mm] ...
    # With zero or more points.
    # Note "D C" is also written for otherwise empty lines (no
    # cylinders in scan)
    out_file = file("cylinders.txt", "w")
    for scan in logfile.scan_data:
        # Find cylinders.
        der = compute_derivative(scan, minimum_valid_distance)
        cylinders = find_cylinders(scan, der, depth_jump,
                                   minimum_valid_distance)
        cartesian_cylinders = compute_cartesian_coordinates(cylinders,
                                                            cylinder_offset)
        # Write to file.
        print >> out_file, "D C",
        for c in cartesian_cylinders:
            print >> out_file, "%.1f %.1f" % c,
        print >> out_file
    out_file.close()
