# Instead of moving a distribution, move (and modify) it using a convolution.
# 06_b_convolve_distribution
# Claus Brenner, 26 NOV 2012
from pylab import plot, show, ylim
from distribution import *
import copy
def move(distribution, delta):
    """Returns a Distribution that has been moved (x-axis) by the amount of
       delta."""
    return Distribution(distribution.offset + delta, distribution.values)

def convolve(a, b):
    """Convolve distribution a and b and return the resulting new distribution."""
    a_ = copy.deepcopy(a)
    b_ = copy.deepcopy(b)
    
    All_length = len(a_.values)+len(b_.values)-1
    for i in xrange(len(a_.values)-1):
        b_.values.append((0))
    for i in xrange(len(b_.values)-1):
        a_.values.append((0))
    # insert enough zero to avoid index overflow

    
    result = Distribution(0,[])
    
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


if __name__ == '__main__':
    arena = (0,100)

    # Move 3 times by 20.
    moves = [20] * 3
    # Start with a known position: probability 1.0 at position 10.
    position = Distribution.unit_pulse(10)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         linestyle='steps')

    # Now move and plot.
    for m in moves:
        move_distribution = Distribution.triangle(m, 2)
        position = convolve(position, move_distribution)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             linestyle='steps')

    ylim(0.0, 1.1)
    show()
