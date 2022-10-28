#
# dynamic GPS signal simulation
# Tang Nguyen Developer
#

import subprocess
import numpy as np
from math import cos, sin
import datetime
import csv
from pymap3d import geodetic2ecef
import matplotlib.pyplot as plt


class Location:
    def __init__(self, latitude, longitude, height, r=1):
        self.__latitude = latitude
        self.__longitude = longitude
        self.__height = height
        self.r = r
        self.__list_point_x = []
        self.__list_point_y = []
        self.__list_point_z = []


    def __del__(self):
        self.__list_point_x.clear()
        self.__list_point_y.clear()
        self.__list_point_z.clear()


    @property
    def list_point(self):
        return np.concatenate((
            np.array(self.__list_point_x),
            np.array(self.__list_point_y),
            np.array(self.__list_point_z)
        ), axis=0)


    def __RotationAngle(self, rad, xyz):
        
        rotation = None
        p = None

        if xyz == 0:
            p = np.array([0, self.r, 0, 1]).reshape(4, 1)
            rotation = np.array([
                [1, 0, 0, 0],
                [0, cos(rad), -sin(rad), 0],
                [0, sin(rad), cos(rad), 0],
                [0, 0, 0, 1]
            ])

        elif xyz == 1:
            p = np.array([0, 0, self.r, 1]).reshape(4, 1)
            rotation = np.array([
                [cos(rad), 0, sin(rad), 0],
                [0, 1, 0, 0],
                [-sin(rad), 0, cos(rad), 0],
                [0, 0, 0, 1]
            ])

        elif xyz == 2:
            p = np.array([self.r, 0, 0, 1]).reshape(4, 1)
            rotation = np.array([
                [cos(rad), -sin(rad), 0, 0],
                [sin(rad), cos(rad), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ])

        else:
            return None
            
        p_new = rotation @ p
        p_new[2] *= 1e5 
        return (np.array([self.__latitude, self.__longitude, self.__height, 0]).reshape(4, 1) * 1e5 + p_new) / 1e5


    def GeneratePoint(self, rotation_angle=45):

        alphas = np.arange(1,360/rotation_angle+1,1)*rotation_angle
        rads = np.deg2rad(alphas)

        for rad in rads:
            p_x = self.__RotationAngle(rad=rad, xyz=0)
            p_y = self.__RotationAngle(rad=rad, xyz=1)
            p_z = self.__RotationAngle(rad=rad, xyz=2)

            self.__list_point_x.append(p_x[:3].reshape(3,))
            self.__list_point_y.append(p_y[:3].reshape(3,))
            self.__list_point_z.append(p_z[:3].reshape(3,))


def GetDateTime():
    utcTimeDelta = datetime.timedelta(hours=7)
    utcTZObject = datetime.timezone(utcTimeDelta, name="utc")
    d1 = datetime.datetime.now(utcTZObject)
    
    day, month, years = d1.day, d1.month, d1.year
    hour, minute, second = d1.hour, d1.minute, d1.second

    dt = [str(day), str(month), str(hour), str(minute), str(second)]
    for i in range(5):
        if len(dt[i]) < 2:
            dt[i] = f"0{dt[i]}"
    
    return f"{years}/{dt[1]}/{dt[0]},{dt[2]}:{dt[3]}:{dt[4]}"


def ShowGlobular():

    lo = Location(0, 0, 0)
    lo.GeneratePoint(1)

    ax = plt.axes(projection='3d')

    xyz = np.array(lo.list_point) * 1e5
    n_point = int(xyz.shape[0]/3)
    xyz = xyz.T

    color = ["blue", "green", "orange"]

    for i in range(3):
        ax.plot3D(xyz[0][i*n_point:(i+1)*n_point], xyz[1][i*n_point:(i+1)*n_point], xyz[2][i*n_point:(i+1)*n_point]/1e5, color[i])
    '''
    ax.plot3D([0, 5,], [0, 5], [0, 0], "black")
    for t in range(6):
        for i in range(3):
            ax.plot3D((xyz[0][i*n_point:(i+1)*n_point])+t, (xyz[1][i*n_point:(i+1)*n_point])+t, xyz[2][i*n_point:(i+1)*n_point]/1e5, color[i])
    '''

    ax.legend(["Rotated by latitude", "Rotated by longitude", "Rotated by height"])
    ax.set_xlabel("Latitude")
    ax.set_ylabel("Longitude")
    ax.set_zlabel("Height")
    ax.set_title("R = 1.112 meters")
    plt.show()


def GenerateDynamicGPSsignal(start_lo, leap):

    alpha = 0.5
    num_Leap = int(1e2)
    duration = int((360/alpha)*num_Leap)*3
    time_slot = 0

    f = open("dynamic_point.csv", "w")
    writer = csv.writer(f)

    for i in range(0, num_Leap):
        current = (start_lo + (leap*i)) / 1e6
        current_location = Location(latitude=current[0], longitude=current[1], height=current[2], r=1e-2)
        current_location.GeneratePoint(alpha)

        list_point = current_location.list_point
        n_point_of_axis = int(list_point.shape[0]/3)
        
        for point in range(n_point_of_axis):
            for axis in range(3):
                p = list_point[point+(axis*n_point_of_axis)]
                x, y, z = geodetic2ecef(p[0], p[1], p[2])
                p = np.array([time_slot / 10, x, y, z])
                writer.writerow(p)
                time_slot += 1

                print(f"\033[A\033[AGenerating {duration} time slot \n Location Point ...", p)

    f.close()


if __name__=="__main__":

    #ShowGlobular()

    
    print("Please wait...\n")

    leap = 1e3 #1e-3 * 1e6
    start_lo = np.array([49.035618, 31.322188, 100.0])*1e6
    GenerateDynamicGPSsignal(start_lo, leap)
    
    
    subprocess.run(["../LimeGPS",
        "-e", "../../brdc3000.22n",
        "-u", "dynamic_point.csv",
        "-T", GetDateTime(),
        "-a", "1.0"
    ])
    

