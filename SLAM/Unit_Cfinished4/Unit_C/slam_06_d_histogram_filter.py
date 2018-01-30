# Histogram implementation of a bayes filter - combines
# convolution and multiplication of distributions, for the
# movement and measurement steps.
# 06_d_histogram_filter
# Claus Brenner, 28 NOV 2012
from pylab import plot, show, ylim
from distribution import *
import copy

def move(distribution, delta):
    """Returns a Distribution that has been moved (x-axis) by the amount of
       delta."""
    return Distribution(distribution.offset + delta, distribution.values)


def convolve(a, b):
    """Convolve distribution a and b and return the resulting new distribution."""
    result = Distribution(0,[])
    a_ = copy.deepcopy(a)
    b_ = copy.deepcopy(b)
    
    All_length = len(a_.values)+len(b_.values)-1
    for i in xrange(len(a_.values)-1):
        b_.values.append((0))
    for i in xrange(len(b_.values)-1):
        a_.values.append((0))
    # insert enough zero to avoid index overflow

    
    # --->>> Put your code here.
    
    # falten
    for i in xrange(All_length):
        summe = 0
        for j in xrange(i+1):
            summe+= a_.values[j]*b_.values[i-j]
        result.values.append((summe))
        
    result = move(result, a_.offset+b_.offset)
    
   
    result.normalize()        
    return result  # Replace this by your own result.

def multiply(a, b):
    """Multiply two distributions and return the resulting distribution."""
    result = Distribution(0,[])
    diff1 = 0.0
    diff2 = 0.0
    
    a_ = copy.deepcopy(a)
    b_ = copy.deepcopy(b)
    
    diff1 = a_.offset-b_.offset

    if diff1>0:
        result.offset = a_.offset
        for i in xrange(diff1):
            a_.values.insert(0,(0))
    elif diff1<0:
        result.offset = b_.offset
        for i in xrange(abs(diff1)):
            b_.values.insert(0,(0))
    else:
        result.offset = a_.offset
        pass

    diff2 = len(a_.values)-len(b_.values)
    if diff2>0:
        for i in xrange(diff2):
            b_.values.append((0))
    elif diff2<0:
        for i in xrange(abs(diff2)):
            a_.values.append((0))
    else:
        pass

    for i in xrange(len(a_.values)):
        if (a_.values[i]*b_.values[i])!=0:
            result.values.append((a_.values[i]*b_.values[i]))
    result.normalize()
    # --->>> Put your code here.
    
    return result  # Modify this to return your result.


# --->>> Copy your convolve(a, b) and multiply(a, b) functions here.



if __name__ == '__main__':
    arena = (0,220)

    # Start position. Exactly known - a unit pulse.
    start_position = 10
    position = Distribution.unit_pulse(start_position)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         linestyle='steps')

    # Movement data.
    controls  =    [ 20 ] * 10

    # Measurement data. Assume (for now) that the measurement data
    # is correct. - This code just builds a cumulative list of the controls,
    # plus the start position.
    p = start_position
    measurements = []
    for c in controls:
        p += c
        measurements.append(p)

    # This is the filter loop.
    for i in xrange(len(controls)):
        # Move, by convolution. Also termed "prediction".
        control = Distribution.triangle(controls[i], 10)
        position = convolve(position, control)#convolve affect the peak more rapidly
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='b', linestyle='steps')

        # Measure, by multiplication. Also termed "correction".
        measurement = Distribution.triangle(measurements[i], 10)
        position = multiply(position, measurement)#multiply affect just the peak
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='r', linestyle='steps')
    #information lose from the movement is balanced
    #from the information gain from measurement
    show()
