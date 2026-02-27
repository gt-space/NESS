import numpy as np
import matplotlib.pyplot as plt

class conicalContour:
    def __init__(self, chmbR, chmbL, contAngle, throatR, throatL, expAngle, exitR, numpts = 50):
        self.chmbR = chmbR
        self.chmbL = chmbL
        self.contAngle = contAngle
        self.throatR = throatR
        self.throatL = throatL
        self.expAngle = expAngle
        self.exitR = exitR
        self.numpts = numpts
        self.z_coords = []
        self.r_coords = []

    def makeContour(self):
        contAngle_rad = np.radians(self.contAngle)
        expAngle_rad = np.radians(self.expAngle)

        x_throat_start = -self.throatL / 2
        x_throat_end = self.throatL / 2

        x_length_conv = (self.chmbR - self.throatR) / np.tan(contAngle_rad)
        x_chamber_end = x_throat_start - x_length_conv

        x_chamber_start = x_chamber_end - self.chmbL

        x_length_div = (self.exitR - self.throatR) / np.tan(expAngle_rad)
        x_exit = x_throat_end + x_length_div

        self.z_coords = np.linspace(x_chamber_start, x_exit, self.numpts)
        self.r_coords = np.zeros(self.numpts)

        # Calculate y-coordinates for each x-coordinate
        for x in self.z_coords:
            if x <= x_chamber_end: # within cylindrical portion of the chamber
                y = self.chmbR
                
            elif x <= x_throat_start: # in the converging section
                m = -np.tan(contAngle_rad)
                y = self.chmbR + m * (x - x_chamber_end)
                
            elif x <= x_throat_end: # within the throat
                y = self.throatR
         
            else:
                # in diverging section
                m = np.tan(expAngle_rad)
                y = self.throatR + m * (x - x_throat_end)
                
            self.r_coords.append(y)
        return self.z_coords, self.r_coords

    def plotContour(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.z_coords, self.r_coords, label='Upper Contour')
        plt.plot(self.z_coords, -np.array(self.r_coords), label='Lower Contour')
        plt.title('Conical Nozzle Contour')
        plt.xlabel('Axial Distance (in)')
        plt.ylabel('Radial Distance (in)')
        plt.grid()
        plt.axis('equal')
        plt.legend()
        plt.show()
    def saveContour(self, filename="conical_contour.csv"):
        data = np.column_stack((self.z_coords, self.r_coords))
        np.savetxt(filename, data, delimiter=",", header="X Position (in), Y Position (in)", comments="")