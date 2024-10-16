import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import pandas as pd
import math
import itertools

class Planet:
    def __init__(self, name, position, velocity, mass):
        self.name = name
        self.position = position
        self.velocity = velocity
        self.mass = mass

b1 = Planet('b1',[2,0],[0,-1],3)
b2 = Planet('b2',[-5,0],[0,0],10)
b3 = Planet('b3',[5,-2],[0,1],10)
b4 = Planet('b4',[0,0],[0,-1],1)
b5 = Planet('b5',[2,0],[0,1],1)
b6 = Planet('b6',[-2,0],[0,0],10)

planets_list = [b1,b2,b3]


# The purpose of the main function it to apply the euler's method in order to model the trajectories or planets. These trajectories are stored in a dataframe.
def main(planets_list, delta, iterations):
    
    grav_constant = -1
    
    df = pd.DataFrame()

    # Adding initial conditions to dataframe
    for planet in planets_list:
        df[str(planet.name) + 'dx'] = [planet.position[0]]
        df[str(planet.name) + 'dy'] = [planet.position[1]]

    for planet in planets_list:
        df[str(planet.name) + 'vx'] = [planet.velocity[0]]
        df[str(planet.name) + 'vy'] = [planet.velocity[1]]

    # Creating seconds column
    df['seconds'] = [0]
    
    # Creating pairs of planets (important when calculating the force between two planets)
    planets_names = [planet.name for planet in planets_list]
    pairs = list(itertools.combinations(planets_names, 2))
    planets_dict = {planet.name: planet for planet in planets_list}

    # Creating empty lines in the dataframe. Number of lines equivalent to number of iterations
    for _ in range(iterations):
        new_row = [None] * len(df.columns)
        df.loc[len(df)] = new_row

    # Creating force columns with empty values

    for pair in pairs:
        p1, p2 = pair
        df[str(p1) + str(p2) + 'fx'] = [0.0] * len(df)
        df[str(p1) + str(p2) + 'fy'] = [0.0] * len(df)
        df[str(p2) + str(p1) + 'fx'] = [0.0] * len(df)
        df[str(p2) + str(p1) + 'fy'] = [0.0] * len(df)

    # This part of the code applies the Euler's method for differential equations in the Newton's force equations F = gm1m2/d*d
    for i in range(iterations):

        # tn+1 = tn + h
        df.loc[i + 1, 'seconds'] = df.loc[i, 'seconds'] + delta
        
        # Calculating distance between pairs of planets and saving them in the dataframe
        for pair in pairs:
            p1, p2 = pair
            df[str(p1) + str(p2) + 'dx'] = df[str(p1) + 'dx'] - df[str(p2) + 'dx']
            df[str(p1) + str(p2) + 'dy'] = df[str(p1) + 'dy'] - df[str(p2) + 'dy']
            distance = math.sqrt((df[str(p1) + str(p2) + 'dx'].iloc[0])**2 + (df[str(p1) + str(p2) + 'dy'].iloc[0])**2)
            
            # Getting unit vector that points from one body to another one
            ux = df[str(p1) + str(p2) + 'dx'].iloc[i] / distance
            uy = df[str(p1) + str(p2) + 'dy'].iloc[i] / distance
    
            # Getting mass for each planet and calculating force between them
            mass1 = planets_dict[p1].mass
            mass2 = planets_dict[p2].mass
    
            force = grav_constant * mass1 * mass2 / (distance**2)
    
            # Saving force value and directions for planet 1 and 2.
            # If the force from A to B is F, the force from B to A is -F
        
            df.loc[i, str(p1) + str(p2) + 'fx'] = force * ux
            df.loc[i, str(p1) + str(p2) + 'fy'] = force * uy
            df.loc[i, str(p2) + str(p1) + 'fx'] = -force * ux
            df.loc[i, str(p2) + str(p1) + 'fy'] = -force * uy

        # Applying Newton method to planet position and velocity
        for planet in planets_list:
            list_x = []
            list_y = []
            other_planets = planets_list.copy()
            other_planets.remove(planet)
            
            
            df.loc[i + 1, str(planet.name) + 'dx'] = df.loc[i, str(planet.name) + 'vx'] * delta + df.loc[i, str(planet.name) + 'dx']
            df.loc[i + 1, str(planet.name) + 'dy'] = df.loc[i, str(planet.name) + 'vy'] * delta + df.loc[i, str(planet.name) + 'dy']

            for other_planet in other_planets:
                list_x.append(df[str(planet.name) + str(other_planet.name) + 'fx'].iloc[i])
                list_y.append(df[str(planet.name) + str(other_planet.name) + 'fy'].iloc[i])

            sum_fx = sum(list_x)
            sum_fy = sum(list_y)

            ax = sum_fx / planet.mass
            ay = sum_fy / planet.mass

            df.loc[i + 1, str(planet.name) + 'vx'] = ax * delta + df.loc[i, str(planet.name) + 'vx']
            df.loc[i + 1, str(planet.name) + 'vy'] = ay * delta + df.loc[i, str(planet.name) + 'vy']
   
    return df

# This part of the code is entirely focused on animating the trajectories of the planets.

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def plot_planet_trajectories(df, planets_list):
    # Figure and axis setup
    fig, axis = plt.subplots()

    # Get the min and max for all planet positions (for x and y limits)
    x_min = min(df[[f'{planet.name}dx' for planet in planets_list]].min())
    x_max = max(df[[f'{planet.name}dx' for planet in planets_list]].max())
    y_min = min(df[[f'{planet.name}dy' for planet in planets_list]].min())
    y_max = max(df[[f'{planet.name}dy' for planet in planets_list]].max())

    axis.set_xlim([x_min, x_max])
    axis.set_ylim([y_min, y_max])

    # Create a Line2D object for each planet
    lines = []
    for planet in planets_list:
        line, = axis.plot([], [], label=planet.name)  # Label for each planet
        lines.append(line)

    # Function to update the plot
    def update_data(frame):
        for i, planet in enumerate(planets_list):
            lines[i].set_data(df[f'{planet.name}dx'][:frame], df[f'{planet.name}dy'][:frame])
        return lines

    # Create animation
    # User can ajust velocity of animation by changing number below
    velocity = 100
    animation = FuncAnimation(fig=fig, func=update_data, frames=len(df), interval=(velocity*df['seconds'].iloc[-1]/len(df)))

    # Add legend to differentiate the planets
    axis.legend()

    # Show animation
    plt.show()

#Running the code

plot_planet_trajectories(main(planets_list, 0.1, 1000), planets_list)
