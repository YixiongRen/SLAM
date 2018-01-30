# Multiply a distribution by another distribution.
# 06_c_multiply_distribution
# Claus Brenner, 26 NOV 2012
from pylab import plot, show
from distribution import *
import copy
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

if __name__ == '__main__':
    arena = (0,1000)

    # Here is our assumed position. Plotted in blue.
    position_value = 400
    position_error = 100
    position = Distribution.triangle(position_value, position_error)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         color='b', linestyle='steps')

    # Here is our measurement. Plotted in green.
    # That is what we read from the instrument.
    measured_value = 410
    measurement_error = 200
    measurement = Distribution.triangle(measured_value, measurement_error)
    plot(measurement.plotlists(*arena)[0], measurement.plotlists(*arena)[1],
         color='g', linestyle='steps')

    # Now, we integrate our sensor measurement. Result is plotted in red.
    position_after_measurement = multiply(position, measurement)
    plot(position_after_measurement.plotlists(*arena)[0],
         position_after_measurement.plotlists(*arena)[1],
         color='r', linestyle='steps')

    show()
