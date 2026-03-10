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
        temp_r_coords = []

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
                
            temp_r_coords.append(y)
        self.r_coords = temp_r_coords
        return self.z_coords, self.r_coords

    def plotContour(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.z_coords, self.r_coords, color = 'blue')
        plt.plot(self.z_coords, -np.array(self.r_coords), color = 'blue')
        plt.title('Conical Nozzle Contour')
        plt.xlabel('Axial Distance (in)')
        plt.ylabel('Radial Distance (in)')
        plt.grid()
        plt.axis('equal')
        plt.legend()
        plt.show()

    def saveContour(self, filename="conical_contourv2.txt"):
        contAngle_rad = np.radians(self.contAngle)
        expAngle_rad = np.radians(self.expAngle)
        x_throat_start = -self.throatL / 2
        x_throat_end = self.throatL / 2
        x_length_conv = (self.chmbR - self.throatR) / np.tan(contAngle_rad)
        x_chamber_end = x_throat_start - x_length_conv
        x_chamber_start = x_chamber_end - self.chmbL
        x_length_div = (self.exitR - self.throatR) / np.tan(expAngle_rad)
        x_exit = x_throat_end + x_length_div

        x_crit = [
            x_chamber_start,  # Chamber start (injector)
            x_chamber_end,    # Chamber end / start of convergence
            x_throat_start,   # Throat start
            x_throat_end,     # Throat end
            x_exit            # Nozzle exit
        ]

        y_crit = [
            self.chmbR,       # Chamber radius
            self.chmbR,       # Chamber radius
            self.throatR,     # Throat radius
            self.throatR,     # Throat radius
            self.exitR        # Exit radius
        ]
        z_crit = np.zeros(len(x_crit))
        data = np.column_stack((x_crit, y_crit, z_crit))
        np.savetxt(
            filename,
            data,
            delimiter=",",
        )
