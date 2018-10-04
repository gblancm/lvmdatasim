from Telescope import Telescope
try:
    myTel = Telescope("LVM160-SCI-S")
    print("Telescope Import Successful")
except:
    print("Telescope Import Failed")
