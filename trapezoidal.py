#############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++Built by Miguel Fernandes Guerreiro+++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++22/07/2019++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#############################################################################
#imports
##import numpy
#https://stackoverflow.com/questions/22690637/numerical-integration-over-non-uniform-grids
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.cumtrapz.html

#function
def trapezoid(p1,p2):
    """
requires:p1[tuple()] and p2[tuple()] with (x,y)
of points of trapezoids, in serial order(p1[0]<p2[0])
ensures:area of trapezoid
"""
    return float(p1[1]+p2[1])*float(p2[0]-p1[0])/2.0


def regression_point_calculation(cut_point,points):
    """calculate a point from a regression line
between 2 points
requires:
    cut_point[int()]:
        coordinate to calculate the regression
        value
    points[list(tuple(),tuple())]:
        2 regression points for cut point
ensures:
    tuple(int(),int()) with coordinates for cut point
"""
    m = (points[1][1]-points[0][1])/float(points[1][0]-points[0][0])
    b = points[1][1]-points[1][0]*m
    return (cut_point,m*cut_point+b)


def regression_points_calc(cut_points,points):
    """
requires:
    cut_points    points to be calculated (only x-axis)
    points        tuple points with data to calculate from (x-axis,y-axis)
ensures:
    list of all points calculated, sorted
"""
    pp = []#calculated points
    for p in cut_points:
        if p in [i[0] for i in points]:
            #print "{} present".format(p)
            continue
        if p<min(points,key=lambda x: x[0])[0]:
            #print "{} below".format(p)
            pp.append(regression_point_calculation(p,points[:2]))
        elif p>max(points,key=lambda x: x[0])[0]:
            #print "{} above".format(p)
            pp.append(regression_point_calculation(p,points[-2:]))
        else:
            #print "{} within".format(p)
            pp.append(regression_point_calculation(p,
                                              [min([i for i in points if i[0]>=p],key=lambda x: x[0]),
                                               max([i for i in points if i[0]<=p], key=lambda x: x[0])])
                      )
        #end of cut points calculation loop
    pp.extend(points)#calculated points plus all the rest
    return sorted(pp,key=lambda x:x[0])


def out_of_bounds_non_uniform_integration(points,cut_points=[100,1000,4000]):
    """
requires:
    points        collection of points(x-axis,y-axis)
    cut_points    collection of points(only x-axis) to calculate cut
ensures:
    List()    integrals for diferent sections(within neighboor cut_points)
    Int()     total integral
"""
    pp = regression_points_calc(cut_points,points)
    ppi = [i[0] for i in pp]
    integrals = []#integrals for diferent sections(within neighboor cut_points)
    Integral = 0.0#total integral
    for cp_s in range(0,len(cut_points)-1):#cp_s start index cut_points. 
        cp_e = cp_s + 1#cp_e end index cut_points

        V_s = cut_points[cp_s]#start value cut_points
        V_e = cut_points[cp_e]#end value cut_points

        pi_s = ppi.index(V_s)#start value calculated points
        pi_e = ppi.index(V_e)#end value calculated points
        if pp[pi_e][1]<0:
            break
        integral = 0.0
        for pi in range(pi_s,pi_e):
            integral += trapezoid(pp[pi],pp[pi+1])
        Integral += integral
        integrals.append(integral)
    return integrals,Integral


if __name__=="__main__":
    points = ((0,10),(3,10))#100,1000,4000
##    d,_ = out_of_bounds_non_uniform_integration(points)
    d = regression_points_calc([0,1,2,3,4,5],points)
    print(d)
