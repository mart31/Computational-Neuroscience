# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:55:58 2016

@author: martin
"""


from scipy.integrate import odeint
from math import exp,erf
import matplotlib.pyplot as plt
import time
import sauvegarde_graphique
import random
import numpy
import os

"""
the alpha and beta value are used in the action potential calculation
"""
def alpha_m(v):
    return 0.10*(v + 40)/(1 - exp(-(v +40)/10))

def beta_m(v):
    return 4 * exp(-(v+65)/18)

def alpha_h(v):
    return 0.07 * exp(-(v + 65)/20)

def beta_h(v):
    return 1/(1 + exp(-(v + 35)/10))

def alpha_n(v):
    return 0.01 * (v + 55) / (1 - exp(-(v + 55)/10))

def beta_n(v):
    return 0.125 * exp(-(v + 65)/80)

"""
define the number of pulse used in the simulation
"""    
def pulses_number():
    number_of_pulses = 1
    return number_of_pulses

"""
the 2 next function are used to generate array containing different value of 
time and length. These array are used in the calculation of the action potential
"""    
def letemps(stoptime,numpoints):
    return [stoptime*float(i)/(numpoints-1) for i in range(numpoints)]
    
def le_length_array(length,numpoints):
    return [length*float(i)/(numpoints-1) for i in range(numpoints)]

"""
the pulse function generate a pulse at a given time that last for a finite time
and have a set current.
"""    
def pulse (time):
    nombres_pulses = pulses_number()
    t = 0
    """
    if the current is not high enough to trigger a ap some fonction of this program
    might not fonction properly. I choose 600 because it just high enough to trigger the
    ap.  
    """
    current = 600.0
    start = 0.0
    duration = 1.2
    interval = 25
    
    temp = (time - start) % interval
    finish = (interval * nombres_pulses) + start
    
    if temp < 0:
        temp = 0
        
    if time > finish:
        t = finish - time
    elif temp < duration:
        t = temp
    else:
        t = duration - temp
    
    return current/2*(1.0 + erf(t/0.001))
    
     
"""
The class Axon is the centerpiece of the program. It creat the axon and enable
us to start a simulation and save the data. 
"""
class Axon(object):
    """
    One important note: the value of Ri is 110ohm*cm but since we use
    milli volt and milli amp we need to use the milli ohm*cm
    """
    Ri = 110e3 # milli ohm*cm 
    Cap = 1e-6 # farad/cm^2
    pi = 3.14159265
    nerve_length = 0.0
    
    def __init__(self,diameter,length):
        """
        These attributes are the defining values of the axon. These values
        are needed in order to compute how the action potential vary in time.
        The internodal length is 142.86*diameter and the diameter is 20 time the
        node length.
        """
        self.diameter = diameter
        self.internodal_length = diameter*142.86
        self.node_length = 1e-4
        self.Ga = ((Axon.pi)*(self.diameter)**2)/(4*(self.internodal_length)*(Axon.Ri))
        self.Cm = (Axon.pi)*(Axon.Cap)*(self.diameter)*(self.node_length)
        self.nodes_ranvier = self.nodes(n = 1+(int(length/(self.node_length+self.internodal_length))))
        self.ac = False
        self.vshift = [0 for x in range(0,len(self.nodes_ranvier)/6,1)] # the left shift for the damaged axon
        Axon.nerve_length = length
        self.solution_axon_potential = 0
        self.node_pulse = self.set_node_pulse()*6 # each node contain 6 value hence the *6
        
        
    
        
    def set_node_pulse(self):
        """
        This function locate the node which will receive the pulse. I chosen it to be at
        1/6 the total length. The reason for this is both end of the axon does not 
        connect to anything and the node near the end generate unreliable data
        """
        return int(len(self.vshift)/6)

    
    def save_axon_sol(self,name):
        """
        to use this function on your computer you need to create a folder named
        axon_sol_csv in exactly the same place as this program is save. It will
        save the data generated  in that folder
        """
        path = "%s.csv" %name
        temp = os.path.abspath(path)
        temp2 = os.path.split(temp)[0]
        filename = "axon_sol_csv\%s" % (path)
        filepath = os.path.join(temp2,filename)
        numpy.savetxt(filepath,self.solution_axon_potential,delimiter=",")
        
    
    def set_damage(self,LS):
        """
        2 important point here. If ac is not True the simulation dont take into
        account the left shift caused by the damage. Second each node have 6 values
        which is why we devide self.node_ranvier by 6. There is only one value
        of vshift per node.
        """
        self.ac = True
        if (len(LS) == len(self.nodes_ranvier)/6):
            self.vshift = LS
        else:
            print "the left shift array was not the right length"
            
    def create_distributed_damage(self,vshift,node,deviation,length):
        """
        I created this function to be able to create damage at a specific node and
        have the damage diminish as you get farther from that node.Vshift is the max 
        damage,node is the node receiving the max damage, deviation is how fast the
        damage diminish and length is the length of the nerve
        """
        d = float(deviation)
        damage_array = [0 for x in range(0,length,1)]
        for i in range(0,length,1):
            ratio = exp((-abs(node-i))/d)
            temp = ratio*vshift
            damage_array[i] = temp
        
        self.set_damage(LS = damage_array)
    
    def show_ap_speed(self,t):
        """
        this function calculate the average speed of the ap.
        """
        first_pulse = False
        index = 0
        ap_number = 0
        temps_total  = 0.0
        temps1 = 0.0  
            
        for i in range(0,len(self.solution_axon_potential), 1):
                
            if self.solution_axon_potential[i][0] > 30.0 and index%2 == 0 and first_pulse == False:
                temps1 = t[i]
                index = index+1
                first_pulse = True
                
            elif self.solution_axon_potential[i][len(self.solution_axon_potential[0])-6] > 30.0 and index%2 == 1:
                temps_total = t[i]-temps1 + temps_total
                index = index+1
                ap_number = ap_number + 1
                
            elif self.solution_axon_potential[i][0] < 30.0 and first_pulse == True:
                # wait for the first pulse to be over before taking another measure
                first_pulse = False
            
        if ap_number == 0:
            average_time = 0
        else:
            average_time = temps_total/ap_number
        #print "the total time is %f" %temps_total
        speed_ap = 0.0
        if average_time > 0.01:
            speed_ap = (Axon.nerve_length*10)/average_time
            
                
        return speed_ap    
        
    def set_sol_ap(self,wsol):
        self.solution_axon_potential = wsol
       
    
    def nodes(self,n):
        """
        This fonction create a array with all the initial value needed for the
        axon potential fonction. n is the number of node and each node have 6 value:
        voltage,m coefficient, h coefficient, n coefficient, m damaged and h damaged.
        These values come from the  Hodgkin–Huxley model.
        """
        l = n*6
        tableau = [0 for x in range(l)]
        vi = -65.4945   
        vshift = 0
        
        mi = alpha_m(vi)/(alpha_m(vi)+ beta_m(vi))
        hi = alpha_h(vi)/(alpha_h(vi)+ beta_h(vi))
        ni = alpha_n(vi)/(alpha_n(vi)+beta_n(vi))
        midamaged = alpha_m(vi + vshift)/(alpha_m(vi + vshift)+beta_m(vi + vshift))
        hidamaged = alpha_h(vi + vshift)/(alpha_h(vi + vshift)+beta_h(vi + vshift))
    
        for i in range(l):
            m = i % 6
            if m == 0:
                tableau[i] = vi
            elif m == 1:
                tableau[i] = mi
            elif m == 2:
                tableau[i] = hi
            elif m == 3:
                tableau[i] = ni
            elif m == 4:
                tableau[i] = midamaged
            elif m == 5:
                tableau[i] = hidamaged
            
        return tableau
    
  
    def axon_potential(self,w,t):
        """
        The method axon_potential is where the calcul for the action potential
        occur. The math is pretty long and is more detailled in my report
        """
        m = [0 for x in range(len(w))]
        gacm = self.Ga/self.Cm
        if self.ac == True:
            ac = 1.0
        else:
            ac = 0
    
        for i in range(0,len(w),6):
            vshift = self.vshift[i/6] # get the damage for the specific node
            if i == 0:
                m[i] = -((1-ac)*120*(w[i+1]**3)*w[i+2]*(w[i]-50) \
                            +ac*120*(w[i+4]**3)*w[i+5]*(w[i]-50) \
                            +36*(w[i+3]**4)*(w[i]+77) \
                            +0.25*(w[i]+54.4) \
                            #-pulse(t)+gacm*(w[i]-w[i+6]))
                            #-gacm*(+w[i+6]-65.4945-2*w[i])) # current sink was added
                            +gacm*(w[i]-w[i+6]))
                m[i+1] = alpha_m(w[i]) * (1 - w[i+1]) - beta_m(w[i]) * w[i+1]
                m[i+2] = alpha_h(w[i]) * (1 - w[i+2]) - beta_h(w[i]) * w[i+2]
                m[i+3] = alpha_n(w[i]) * (1 - w[i+3]) - beta_n(w[i]) * w[i+3]
                m[i+4] = alpha_m(w[i]+ vshift) * (1 - w[i+4]) - beta_m(w[i]+ vshift) * w[i+4]
                m[i+5] = alpha_h(w[i]+ vshift) * (1 - w[i+5]) - beta_h(w[i]+ vshift) * w[i+5]
            elif i == (len(w)-6):
                m[i] = -((1-ac)*120*(w[i+1]**3)*w[i+2]*(w[i]-50) \
                            +ac*120*(w[i+4]**3)*w[i+5]*(w[i]-50) \
                            +36*(w[i+3]**4)*(w[i]+77) \
                            +0.25*(w[i]+54.4)\
                            #-(gacm*(w[i-6]-65.4945-2*w[i]))) # current sink was added
                            -(gacm*(w[i-6]-w[i])))
                m[i+1] = alpha_m(w[i]) * (1 - w[i+1]) - beta_m(w[i]) * w[i+1]
                m[i+2] = alpha_h(w[i]) * (1 - w[i+2]) - beta_h(w[i]) * w[i+2]
                m[i+3] = alpha_n(w[i]) * (1 - w[i+3]) - beta_n(w[i]) * w[i+3]
                m[i+4] = alpha_m(w[i] + vshift) * (1 - w[i+4]) - beta_m(w[i]+ vshift) * w[i+4]
                m[i+5] = alpha_h(w[i] + vshift) * (1 - w[i+5]) - beta_h(w[i]+ vshift) * w[i+5]
            elif i == self.node_pulse:
                m[i] = -((1-ac)*120*(w[i+1]**3)*w[i+2]*(w[i]-50) \
                            +ac*120*(w[i+4]**3)*w[i+5]*(w[i]-50) \
                            +36*(w[i+3]**4)*(w[i]+77) \
                            +0.25*(w[i]+54.4) \
                            -pulse(t)-(gacm*(w[i-6]+w[i+6]-2*w[i])))
                m[i+1] = alpha_m(w[i]) * (1 - w[i+1]) - beta_m(w[i]) * w[i+1]
                m[i+2] = alpha_h(w[i]) * (1 - w[i+2]) - beta_h(w[i]) * w[i+2]
                m[i+3] = alpha_n(w[i]) * (1 - w[i+3]) - beta_n(w[i]) * w[i+3]
                m[i+4] = alpha_m(w[i] + vshift) * (1 - w[i+4]) - beta_m(w[i]+ vshift) * w[i+4]
                m[i+5] = alpha_h(w[i] + vshift) * (1 - w[i+5]) - beta_h(w[i]+ vshift) * w[i+5]
            else:
                m[i] = -((1-ac)*120*(w[i+1]**3)*w[i+2]*(w[i]-50) \
                            +ac*120*(w[i+4]**3)*w[i+5]*(w[i]-50) \
                            +36*(w[i+3]**4)*(w[i]+77) \
                            +0.25*(w[i]+54.4) \
                            -(gacm*(w[i-6]+w[i+6]-2*w[i])))
                m[i+1] = alpha_m(w[i]) * (1 - w[i+1]) - beta_m(w[i]) * w[i+1]
                m[i+2] = alpha_h(w[i]) * (1 - w[i+2]) - beta_h(w[i]) * w[i+2]
                m[i+3] = alpha_n(w[i]) * (1 - w[i+3]) - beta_n(w[i]) * w[i+3]
                m[i+4] = alpha_m(w[i] + vshift) * (1 - w[i+4]) - beta_m(w[i]+ vshift) * w[i+4]
                m[i+5] = alpha_h(w[i] + vshift) * (1 - w[i+5]) - beta_h(w[i]+ vshift) * w[i+5]
        
        return m




"""
The nerve is used to create a bundle of nerve and start the simulation. There
are some function in this class that are used to visualize the result of the simulation
"""
class Nerve(object):
   
    """
    There are 3 value that are essential in this class. First the length of the nerve.
    Second the diameter_array of all the axon. Third the time array which containt a
    set amount of value all separated exactly by delta t. The time array is with
    letemps function at the beginning of this module.
    """
    def __init__(self,length,diameter_array):
        self.length = length
        self.number_axon = len(diameter_array)
        self.axon_array = self.create_axons(diameters = diameter_array,length = length)
        
        
    def create_axons(self,diameters,length):
        axons = [Axon(diameter = diameters[x],length = length) for x in range(len(diameters))]
        return axons
    
    def start_simulation(self,t):
        debut = time.time()
        temps = t

        for x in range(0,self.number_axon,1):
            sol = odeint(self.axon_array[x].axon_potential,self.axon_array[x].nodes_ranvier,temps,atol = 1.0e-12,rtol = 1.0e-12)
            self.axon_array[x].set_sol_ap(wsol = sol)           
        
        fin = time.time()
        print "the time it took to execute the simulation: %f second" %(fin-debut)
    
    """
    the plot_node_ap is used to graph the different ap at the different node. This
    function was used mainly to see if the action potential was travelling correctly
    through the axon. t is again the time array, start and end are the beginning and
    end of the plot.
    """
    def plot_node_ap(self,t,start,end,axon,name = "none"):
        temps = t
        a = axon
      
        for i in range(0,len(self.axon_array[a].solution_axon_potential[0]),6):
            plt.plot(temps,self.axon_array[a].solution_axon_potential[:,i])
            
        plt.legend(loc = 'upper left')
        plt.title('V par rapport au temps')
        plt.ylabel('V (mV)')
        plt.xlabel('le temps en mS')
        plt.xlim([start,end])
        
        if name != "none":
            sauvegarde_graphique.save(name, ext="png", close=False, verbose=True)        
        
        plt.show()
       
    """
    m,h and n are the coefficient needed to produce the action potential of the 
    Hodgkin–Huxley model. This function enable us to see these value fluctuate in
    time. Like the plot node ap you have to choose the specific node and also
    choose the axon you want to look at.
    """
    def show_graph_mhn(self,start,end,temps,n,axon,name = "none"):
        # node must have value 0,6,12,18,24 etc...
        wsol = self.axon_array[axon].solution_axon_potential
        if n == "last":
            node = len(wsol[0])-6
        else:
            node = n*6        
        
        
        plt.figure(self.number_axon+1) 
        plt.plot(temps, wsol[:,node+1]**3*wsol[:,node+2], label = 'm3h')
        plt.plot(temps, wsol[:,node+1], label = 'm')
        plt.plot(temps, wsol[:,node+2], label = 'h')
        plt.plot(temps, wsol[:,node+3], label = 'n')
        plt.plot(temps, wsol[:,node+3]**4, label = 'n4')
        plt.legend(loc = 'upper right')
        plt.xlabel('Time, t (ms)')
        plt.ylabel('Nav open probability')
        plt.xlim([start,end])
        if name != "none":
            sauvegarde_graphique.save(name, ext="png", close=False, verbose=True) 
        
        plt.show()      
    
    
    def add_damage(self,LS,axon):
        """
        LS is a array containing a left shift value for each node of the 
        specified axon
        """
        # the damage is binary
        self.axon_array[axon].set_damage(LS = LS)
        
    def create_random_damage(self,intensity = None):
        """
        this function use the create distributed damage function to damage
        randomly the nerve
        """
        if intensity == None:
            dmg_amplifier = 1
        else:
            dmg_amplifier = intensity
            print "damage intensity is: %f" %intensity
            
        vmax = int(11*dmg_amplifier)
       
        
        for x in range(0,self.number_axon,1):
            length = len(self.axon_array[x].nodes_ranvier)/6
            dmax = int(length*dmg_amplifier/2)            
            
            v = random.randrange(0, vmax, 1)/10.0
            n = random.randrange(0, length, 1)
            if dmax > 2:
                d = random.randrange(1, dmax, 1)
            else:
                d = random.randrange(1,2,1)
            
            
            self.axon_array[x].create_distributed_damage(vshift = v,node = n,deviation = d,length = length)
            
    """
    The function graph_nerve_damage help visualize the amounth of damage that is
    present at every point in the nerve
    """       
    def graph_nerve_damage(self):
        numpoints = 1000
        l = le_length_array(self.length,numpoints)
        dmg = [0 for i in range(0,numpoints,1)]
        for y in range(0,numpoints,1):
            temp = 0.0
            for x in range(0,self.number_axon,1):
                node = len(self.axon_array[x].nodes_ranvier)/6
                divider = numpoints/node
                index = int(y/divider)
                if index >= node:
                    index = node-1
                
                temp = self.axon_array[x].vshift[index] + temp
            dmg[y] = temp/self.number_axon
        
        plt.plot(l,dmg)
        plt.show()                            
      
    """
    The next 2 function are the least user friendly of the module. If you want to use
    these 2 function without modifying them you need to create a axon_sol_csv folder in
    the same location that this program is located
    """
    def save_axon_damage(self):
        length_max = len(self.axon_array[0].vshift)
        dmg = [[0 for y in range(0,length_max,1)] for x in range(0,self.number_axon,1)]
        
        for x in range(0,self.number_axon,1):
            for y in range(0,len(self.axon_array[x].vshift),1):
                dmg[x][y] = self.axon_array[x].vshift[y]
        
        """
        numpy.savetxt can save a 2d array only if the length of every row
        is the same. To solve this i created artificial node of value 0
        """
        
        path = "damage_array_%dcm.csv" %int(self.length) 
        temp = os.path.abspath(path)
        temp2 = os.path.split(temp)[0]
        filename = "axon_sol_csv\%s" % (path)
        filepath = os.path.join(temp2,filename)
        numpy.savetxt(filepath,dmg,delimiter=",")
    
    def save_simulation_result(self):
        temp = int(self.length)
        for x in range(0,len(test_nerve.axon_array),1):
            filename = "axon%d_%f_sol_%dcm" %(x,test_nerve.axon_array[x].diameter,temp)
            test_nerve.axon_array[x].save_axon_sol(name = filename)

    
        
"""
This section initialize a nerve and the values needed for the simulation. There 3
critical values needed for the simulation. First the time array, second the diameters array
and third and the length of the nerve.
"""     
temps = letemps(10,2000)  # the time array 8 millisecond 20000 interval
#d1 = [15e-4,15.5e-4,16e-4,16.5e-4,17e-4,17.5e-4,18e-4,18.5e-4,19e-4,19.5e-4] 
#d2 = [20e-4,20.5e-4,21e-4,21.5e-4,22e-4,22.5e-4,23e-4,23.5e-4,24e-4,24.5e-4]

d1 = [8e-4,8e-4,8e-4,8e-4,8e-4,11e-4,11e-4,11e-4,11e-4,11e-4] 
d2 = [20e-4,20e-4,20e-4,20e-4,20e-4,25e-4,25e-4,25e-4,25e-4,25e-4]
diameters = d1 + d2 # the diameters of all the axon in the nerve
#diameters = [8e-4]
l = 30.0 # the length of the nerve

test_nerve = Nerve(length = l,diameter_array = diameters) # this is our nerve

"""
the next part add some damage to the nerve. The for loop create constant damage 
at a precise location on the nerve. The create random damage does exactly as it said
put both into comment if you want a intact nerve
"""
#test_nerve.create_random_damage(intensity = 1.5)

for y in range(0,test_nerve.number_axon,1):
    dmg = [0.0 for x in range(0,len(test_nerve.axon_array[y].vshift))]
    third = int(len(test_nerve.axon_array[y].vshift)/3)
    for z in range(third,len(test_nerve.axon_array[y].vshift)-third,1):
        dmg[z] = 1.0
    test_nerve.axon_array[y].set_damage(LS = dmg)
        


"""
this next part start the simulation and save the result
"""
test_nerve.start_simulation(t = temps) 
test_nerve.save_simulation_result()
test_nerve.save_axon_damage()


"""
This part plot the individual action potential of every node for every action.
Note that in the data analysis module we cut the first 1/6 and last 1/6 of the nerve
so the data might look different. We cut both end in order to remove the error generated
by having a dead end in the axon.
"""
for x in range(0,test_nerve.number_axon,1):
    print "\n node ap of axon: %d" %x
    test_nerve.plot_node_ap(t = temps,start = 0.0,end = 8.0,axon = x)
   

for x in range(0,test_nerve.number_axon,1):
    speed_ap = test_nerve.axon_array[x].show_ap_speed(t = temps)
    print "axon %d ap speed: %f m/s" %(x,speed_ap)







