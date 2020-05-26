#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on Thu May 14 11:25:05 2020

#@author: pengwah

from vpython import *
import math
import random
import time

PI = math.pi
allAtoms = []
allWalls = []
allOutlines = []
scene.autoscale = False
temperature = 100
numCollisions = 0
totalForce = 0
memory = []
timestamps = []

def updatePressure(a):
    global memory
    global timestamps
    if (len(memory) < 100):
        memory.append(a)
        timestamps.append(time.time())
    else:
        memory[(numCollisions-1) % 100] = a
        timestamps[(numCollisions-1) % 100] = time.time()



scene.background = color.white


class Atom:
    velocity = None
    nucleus = None
    
    xy = 0
    z = 0
    velocity = 0
    def __init__(self, m, x,y,z):
        a = random.randrange(1,100)
        self.nucleus = sphere(pos = vector(x,y,z), radius = m**(1/3), color = color.black)
        self.xy = random.randrange(0,360) * PI/180
        self.z = random.randrange(-90,90) * PI/180
        standard_v = getStandardV(m**(1/3))
        v = getActualV(standard_v)/40
        self.velocity = vector(v * math.cos(self.xy)*math.cos(self.z), v * math.sin(self.xy)*math.cos(self.z), v * math.sin(self.z))
    def move(self):
        
        self.nucleus.pos = self.nucleus.pos + self.velocity
    
    def remove(self):
        self.nucleus.visible = False
        del self.nucleus


x_d = 120
y_d = 120
z_d = 120
thickness = min(x_d, y_d, z_d)/240
 


colorDict = {"hydrogen": color.black, "helium": color.red, "nitrogen": color.orange, "oxygen": color.yellow, "water": color.cyan, "carbon dioxide": color.black, "methane": color.green}
        


class Hydrogen(Atom):

    def __init__(self):
        super(Hydrogen, self).__init__(1, random.randrange(100, x_d*1000 - 100)/1000, random.randrange(100, y_d*1000 - 100)/1000, random.randrange(100, z_d*1000 - 100)/1000)
        self.nucleus.color = colorDict["hydrogen"]
        
class Helium(Atom):

    def __init__(self):
        super(Helium, self).__init__(4, random.randrange(100, x_d*1000 - 100)/1000, random.randrange(100, y_d*1000 - 100)/1000, random.randrange(100, z_d*1000 - 100)/1000)
        self.nucleus.color = colorDict["helium"]
        

class Nitrogen(Atom):

    def __init__(self):
        super(Nitrogen, self).__init__(14, random.randrange(100, x_d*1000 - 100)/1000, random.randrange(100, y_d*1000 - 100)/1000, random.randrange(100, z_d*1000 - 100)/1000)
        self.nucleus.color = colorDict["oxygen"]
        
class Oxygen(Atom):

    def __init__(self):
        super(Oxygen, self).__init__(16, random.randrange(100, x_d*1000 - 100)/1000, random.randrange(100, y_d*1000 - 100)/1000, random.randrange(100, z_d*1000 - 100)/1000)
        self.nucleus.color = colorDict["oxygen"]
        


class Water(Atom):
    
    def __init__(self):
        super(Water, self).__init__(18, random.randrange(100, x_d*1000 - 100)/1000, random.randrange(100, y_d*1000 - 100)/1000, random.randrange(100, z_d*1000 - 100)/1000)
        self.nucleus.color = colorDict["water"]


class CO2(Atom):
    
    def __init__(self):
        super(CO2, self).__init__(44, random.randrange(100, x_d*1000 - 100)/1000, random.randrange(100, y_d*1000 - 100)/1000, random.randrange(100, z_d*1000 - 100)/1000)
        self.nucleus.color = colorDict["carbon dioxide"]
        
class Methane(Atom):
    
    def __init__(self):
        super(Methane, self).__init__(16, random.randrange(100, x_d*1000 - 100)/1000, random.randrange(100, y_d*1000 - 100)/1000, random.randrange(100, z_d*1000 - 100)/1000)
        self.nucleus.color = colorDict["methane"]


class Outline:
    line = None
    def __init__(self, pos, size, color):
        self.line = box(pos = pos, size = size, color = color)

    def remove(self):
        self.line.visible = False
        del self.line
    
def make_outlines():
    bot0 = Outline(pos = vector(x_d/2, 0, 0), size = vector(x_d, thickness, thickness), color = color.black)
    allOutlines.append(bot0)
    bot1 = Outline(pos = vector(x_d, y_d/2, 0), size = vector(thickness, y_d, thickness), color = color.black)
    allOutlines.append(bot1)
    bot2 = Outline(pos = vector(x_d/2, y_d, 0), size = vector(x_d, thickness, thickness), color = color.black)
    allOutlines.append(bot2)
    bot3 = Outline(pos = vector(0, y_d/2, 0), size = vector(thickness, y_d, thickness), color = color.black)
    allOutlines.append(bot3)

    mid0 = Outline(pos = vector(0, 0, z_d/2), size = vector(thickness, thickness, z_d), color = color.black)
    allOutlines.append(mid0)
    mid1 = Outline(pos = vector(x_d, 0, z_d/2), size = vector(thickness, thickness, z_d), color = color.black)
    allOutlines.append(mid1)
    mid2 = Outline(pos = vector(x_d, y_d, z_d/2), size = vector(thickness, thickness, z_d), color = color.black)
    allOutlines.append(mid2)
    mid3 = Outline(pos = vector(0, y_d, z_d/2), size = vector(thickness, thickness, z_d), color = color.black)
    allOutlines.append(mid3)
    
    top0 = Outline(pos = vector(x_d/2, 0, z_d), size = vector(x_d, thickness, thickness), color = color.black)
    allOutlines.append(top0)
    top1 = Outline(pos = vector(x_d, y_d/2, z_d), size = vector(thickness, y_d, thickness), color = color.black)
    allOutlines.append(top1)
    top2 = Outline(pos = vector(x_d/2, y_d, z_d), size = vector(x_d, thickness, thickness), color = color.black)
    allOutlines.append(top2)
    top3 = Outline(pos = vector(0, y_d/2, z_d), size = vector(thickness, y_d, thickness), color = color.black)
    allOutlines.append(top3)

    




class Wall:
    orientation = None
    positive = None
    wall = None
    thickness = 0.01
    def __init__(self, o, p):
        self.orientation = o
        self.positive = p
        if (o == "xy" and  p == True):
            self.wall = box(pos = vector(x_d/2, y_d/2, z_d), size = vector(x_d, y_d, self.thickness), color = color.blue)
        elif (o == "xy" and  p == False):
            self.wall = box(pos = vector(x_d/2, y_d/2, 0), size = vector(x_d, y_d, self.thickness), color = color.blue)
        elif (o == "yz" and  p == True):
            self.wall = box(pos = vector(x_d, y_d/2, z_d/2), size = vector(self.thickness, y_d, z_d), color = color.blue)
        elif (o == "yz" and  p == False):
            self.wall = box(pos = vector(0, y_d/2, z_d/2), size = vector(self.thickness, y_d, z_d), color = color.blue)
        elif (o == "xz" and  p == True):
            self.wall = box(pos = vector(x_d/2, y_d, z_d/2), size = vector(x_d, self.thickness, z_d), color = color.blue)
        elif (o == "xz" and  p == False):
            self.wall = box(pos = vector(x_d/2, 0, z_d/2), size = vector(x_d, self.thickness, z_d), color = color.blue)
        self.wall.opacity = 0.2
  
    def remove(self):
        self.wall.visible = False
        del self.wall



def start_pause(b):
    if (b.text == "Start simulation"):
        b.text = "Pause simulation"
    else:
        b.text = "Start simulation"
start_button = button(text="Start simulation", bind=start_pause)
scene.append_to_caption("\n\n")



def reset(b):
    numAtoms = len(allAtoms)
    for x in range(numAtoms):
        target = allAtoms.pop()
        target.remove()
    for x in range(12):
        target = allOutlines.pop()
        target.remove()
    for x in range(6):
        target = allWalls.pop()
        target.remove()
    print("Collisions: " + str(numCollisions) + ", total force: " + str(totalForce))
reset_button = button(text="Clear simulation", bind=reset, disabled=True)


#Making the menu
def M(m):
    pass
#    if (m.selected is None or m.selected != "Select a molecule or create your own using the slider:"):
 #   add_atom.disabled = False
molecule_menu = menu( choices=['Select a molecule or create your own using the slider:', 'Custom', 'Hydrogen', 'Helium', 'Nitrogen', 'Oxygen', 'Water', 'Carbon Dioxide', 'Methane'], bind=M, selected='Select a molecule or create your own using the slider:')
#scene.append_to_caption('\n\n')




def setsize(s):
    wt.text = "Custom particle mass: " + str(s.value) + " amu"
    
mass_slider = slider(min=1, max=200, value=1, length=220, step = 1, bind=setsize, right=15)

wt = wtext(text="Custom particle mass: " + str(mass_slider.value) + " amu")



scene.append_to_caption('\n\n')
              
#I know a few of these aren't atoms but i couldn't care less i only decided to add molecules after everything was named "atom" and i ain't taking the time to fix that
def add_atom(b):
    atomq = None
    if (molecule_menu.selected == "Hydrogen"):
        atomq = Hydrogen()
    elif (molecule_menu.selected == "Helium"):
        atomq = Helium()
    elif (molecule_menu.selected == "Nitrogen"):
        atomq = Nitrogen()
    elif (molecule_menu.selected == "Oxygen"):
        atomq = Oxygen()
    elif (molecule_menu.selected == "Water"):
        atomq = Water()
    elif (molecule_menu.selected == "Carbon Dioxide"):
        atomq = CO2()
    elif (molecule_menu.selected == "Methane"):
        atomq = Methane()
    elif (molecule_menu.selected == "Custom"):
        atomq = Atom(mass_slider.value, random.randrange(100, x_d*1000 - 100)/1000, random.randrange(100, y_d*1000 - 100)/1000, random.randrange(100, z_d*1000 - 100)/1000)
    
    if(atomq is not None):
        #print(atomq.nucleus.pos)
    #atom = Atom(0.2,x_d/2,y_d/2,z_d/2)
        allAtoms.append(atomq)
    
     
def remove_atom(b):
    try:
        target = allAtoms.pop(len(allAtoms) - 1)
        target.remove()
    except:
        pass
    

add_atom = button( bind=add_atom, text="Add atom", disabled=True)
remove_atom = button(bind = remove_atom, text="Remove atom", disabled=True)




numAtoms = wtext(text="Number of atoms: " + str(len(allAtoms)))



scene.append_to_caption('\n\n')






def x_toggle(s):
    print(len(allWalls))
    xlabel.text = "Container length: " + str(s.value)
    global x_d
    x_d = s.value
    for atom in allAtoms:
        if (atom.nucleus.pos.x > x_d):
            atom.nucleus.pos.x = .9 * x_d
    
    
x_slider = slider(min=50, max=100, value=120, length=220, step = 1, bind=x_toggle, right=15)
xlabel = wtext(text="Container length: " + str(x_slider.value)+ "           ")



def temp_toggle(s):
    global temperature
    templabel.text = "Temperature: " + str(temp_slider.value)
    differential = s.value/temperature
    temperature = s.value
    for atom in allAtoms:
        atom.velocity *= (differential**.5)
        print(atom.velocity)
temp_slider = slider(min=10, max=200, value=100, length=220, step=1, bind=temp_toggle, right=15)
templabel = wtext(text="Temperature: " + str(temp_slider.value), right = 15)



scene.append_to_caption("\n\n")



def y_toggle(s):
    print(len(allWalls))
    ylabel.text = "Container height: " + str(s.value)
    global y_d
    y_d = s.value
    for atom in allAtoms:
        if (atom.nucleus.pos.y > y_d):
            atom.nucleus.pos.y = .9 * y_d
    
    
    
y_slider = slider(min=50, max=100, value=120, length=220, step = 1, bind=y_toggle, right=15)
ylabel = wtext(text="Container height: " + str(y_slider.value))
volume = wtext(text = "                                                     Volume: " + str(x_d * y_d * z_d))
scene.append_to_caption("\n\n")





def z_toggle(s):
    #print(len(allWalls))
    zlabel.text = "Container width: " + str(s.value)
    global z_d
    z_d = s.value
    for atom in allAtoms:
        if (atom.nucleus.pos.z > z_d):
            atom.nucleus.pos.z = .9 * z_d
    
    
z_slider = slider(min=50, max=100, value=120, length=220, step = 1, bind=z_toggle, right=15)
zlabel = wtext(text="Container width: " + str(z_slider.value) + "           ", right = 15)



pressure_ratio = wtext(text="                                           Pressure ratio: 0")


scene.append_to_caption("\n\n")








def calculateSA():
    return ((x_d * y_d) + (x_d * z_d) + (y_d * z_d))*2


#gets the exact average speed for a molecule with certain radius
def getStandardV(radius):
    boltz = 1.38 * (10**-23)
    kinetic_energy = 1.5 * boltz * temperature
    kinetic_energy *= 2
    kinetic_energy /= (radius ** 3)
    speed = kinetic_energy ** .5
    print("Speed:", speed * (10 ** 12))
    return (speed * (10 ** 12))


#varies molecule speed because of normal curve
def getActualV(v):
    new_v = int(v)
    standard_dev = int(.01*v)
    #print("Deviation:", standard_dev)
    rng = random.randint(0, 999)
    print(rng)
    if (rng < 23):
        multiplier = -3
    elif (rng < 160):
        multiplier = -2
    elif (rng < 500):
        multiplier = -1
    elif (rng >= 977):
        multiplier = 2
    elif (rng >= 840):
        multiplier = 1
    elif (rng >= 500):
        multiplier = 0
    print("Dev:", multiplier * standard_dev)
    finalRng = random.randint(new_v + (multiplier * standard_dev), new_v + ((multiplier + 1) * standard_dev))
    print("Final speed:", finalRng)
    return finalRng





        




        

print("making wall")
#wall = box(pos= vector(6,0,0), size = vector(0.2, 12, 12), color = color.blue)
"""
ball = sphere(pos=vector(0,0,0), radius = 1, color=color.green)
ball2 = sphere(pos=vector(6,0,0), radius = 1, color=color.blue)
ball3 = sphere(pos=vector(0,6,0), radius = 1, color=color.white)
ball3 = sphere(pos=vector(0,0,6), radius = 1, color=color.yellow)

"""




def makeWalls():
    for x in (["xy", "yz", "xz"]):
        for y in ([True, False]):
            yeet = Wall(x,y)
            allWalls.append(yeet)

def wallCollision(atom, wallObj):
    global numCollisions
    global totalForce
    
    
    if (wallObj.orientation == "xy"):
        atom.velocity.z *= -1
        numCollisions+= 1
        
        
        momentum = (atom.nucleus.radius ** 3) * abs(atom.velocity.z) * 2
        updatePressure(momentum)
        
        
        #totalForce += momentum 
    elif (wallObj.orientation == "yz"):
        atom.velocity.x *= -1
        numCollisions+= 1
        momentum = (atom.nucleus.radius ** 3) * abs(atom.velocity.x) * 2
        updatePressure(momentum)
        #totalForce += momentum
    elif (wallObj.orientation == "xz"):
        atom.velocity.y *= -1
        numCollisions+= 1
        momentum = (atom.nucleus.radius ** 3) * abs(atom.velocity.y) * 2
        updatePressure(momentum)
        #totalForce += momentum
    
    
    
    
    
    """
def checkWallCollides():
    for sph in allBalls:
        for wal in allWalls:
            if (mag(sph.ball.pos - wal.wall.pos) <= sph.ball.radius):
                wallCollision(sph.ball, wal)
"""

def checkWallCollides():    
    for atom in allAtoms:
        currentNucleus = atom.nucleus
        for wallObj in allWalls:
            currentWall = wallObj.wall
            if (wallObj.orientation == "xy"):
                if (currentWall.pos.z == 0 and atom.velocity.z < 0 and atom.nucleus.pos.z < atom.nucleus.radius):
                    wallCollision(atom, wallObj)
                elif (currentWall.pos.z == z_d and atom.velocity.z > 0 and atom.nucleus.pos.z > (z_d - atom.nucleus.radius)):
                    wallCollision(atom, wallObj)
            elif (wallObj.orientation == "yz"):
                if (currentWall.pos.x == 0 and atom.velocity.x < 0 and atom.nucleus.pos.x < atom.nucleus.radius):
                    wallCollision(atom, wallObj)
                elif (currentWall.pos.x == x_d and atom.velocity.x > 0 and atom.nucleus.pos.x > (x_d - atom.nucleus.radius)):
                    wallCollision(atom, wallObj)
            elif (wallObj.orientation == "xz"):
                if (currentWall.pos.y == 0 and atom.velocity.y < 0 and atom.nucleus.pos.y < atom.nucleus.radius):
                    wallCollision(atom, wallObj)
                elif (currentWall.pos.y == y_d and atom.velocity.y > 0 and atom.nucleus.pos.y > (y_d - atom.nucleus.radius)):
                    wallCollision(atom, wallObj)


def elastic(m1, v1, m2, v2):
    final1 = (((m1-m2)/(m1 + m2))*v1)  +  ((2*m2/(m1 + m2)) * v2)
    final2 = ((2*m1/(m1 + m2))*v1) + (((m2-m1)/(m1 + m2)) * v2)
    return[final1, final2]

def atomCollision(atom1, atom2):
    if (mag(atom1.nucleus.pos - atom2.nucleus.pos) <= atom1.nucleus.radius + atom2.nucleus.radius):
        #print("Collided!")
        m1 = atom1.nucleus.radius ** 3
        m2 = atom2.nucleus.radius ** 3
        
        #x
        elasticx = elastic(m1, atom1.velocity.x, m2, atom2.velocity.x)
        atom1.velocity.x = elasticx[0]
        atom2.velocity.x = elasticx[1]
        #y
        elasticy = elastic(m1, atom1.velocity.y, m2, atom2.velocity.y)
        atom1.velocity.y = elasticy[0]
        atom2.velocity.y = elasticy[1]
        #z
        elasticz = elastic(m1, atom1.velocity.z, m2, atom2.velocity.z)
        atom1.velocity.z = elasticz[0]
        atom2.velocity.z = elasticz[1]
        


def checkAtomCollides():
    if (len(allAtoms) > 0):
        for atom1 in range(len(allAtoms)):
            for atom2 in range(atom1 + 1, len(allAtoms)):
                atomCollision(allAtoms[atom1], allAtoms[atom2])
            





makeWalls()
make_outlines()


while(True):
    rate(30)
    
    
    if (len(memory) == 0 or len(timestamps) == 1):
        pressure_ratio.text = "                                           Pressure ratio: 0"
    else:
        pressure_ratio.text = "                                           Pressure ratio: " + str(1000000*(sum(memory)/len(memory))/(max(timestamps) - min(timestamps))/calculateSA())
    numAtoms.text = "Moles of gas present: " + str(len(allAtoms))
    volume.text = "                                                     Volume: " + str(x_d * y_d * z_d)

    for x in range(6):
        try:
            target = allWalls.pop()
            target.remove()
        except:
            print("yeet")
    for x in range(12):
        try:
            target = allOutlines.pop()
            target.remove()
        except:
            print("no")
    makeWalls()
    make_outlines()
    if(start_button.text == "Pause simulation"):
        reset_button.disabled = True
        #print(molecule_menu.selected) 
        #print(type(molecule_menu.selected))
        if (molecule_menu.selected is None or molecule_menu.selected != "Select a molecule or create your own using the slider:"):
            add_atom.disabled = False
        if (len(allAtoms) > 0):
            #print(allAtoms[0].nucleus.pos, allAtoms[0].velocity)
            remove_atom.disabled = False
        for atom in allAtoms:
            atom.move()
        checkWallCollides()
        checkAtomCollides()
        #print(mass_slider.value)
        #rgb = wtext(pos=scene.title_anchor, text="  Mass of custom molecule: " + str(mass_slider.value) + " ")
    else:
        reset_button.disabled = False
