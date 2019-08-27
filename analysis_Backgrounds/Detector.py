#!/Users/Lenardo/anaconda/envs/py36/bin/python

class XenonDetector:
  # Natural abundances
  isotopes = {
     '124':0.000952,
     '126':0.000890,
     '128':0.019102,
     '129':0.264006,
     '130':0.040710,
     '131':0.21232,
     '132':0.269086,
     '134':0.104257,
     '136':0.088573,
  }

  def __init__( self ):
    self.distance = 100 #m to core of reactor
    self.mass     = 1000 #kg of xenon
    self.height   = 1 #m
    self.width    = 1 #m
    self.ComputeSolidAngleFraction()

  def ComputeSolidAngleFraction(self):
    # We'll approximate that the detector is far 
    # enough away for small-angle approximations
    surfaceArea = (self.distance * 100)**2 * 4 * 3.1415927 # cm^2
    self.fluxFactor = 1./surfaceArea
    self.solidAngleFraction = (self.height * 100 * self.width * 100) / surfaceArea





class ArgonDetector:
  isotopes = {
    '36':0.003373,
    '40':0.996627,
  }

  def __init__( self ):
    self.distance = 100
    self.mass = 1000
    self.height = 1
    self.width = 1
    self.ComputeSolidAngleFraction()

  def ComputeSolidAngleFraction(self):
    # We'll approximate that the detector is far 
    # enough away for small-angle approximations
    surfaceArea = (self.distance*100)**2 * 4 * 3.1415927 #cm^2
    self.fluxFactor = 1./surfaceArea # cm^-2
    self.solidAngleFraction = (self.height * 100 * self.width * 100) / surfaceArea


class GeDetector:
  isotopes = {
    '70':0.2038,
    '72':0.2731,
    '73':0.0776,
    '74':0.3672,
    '76':0.0786,    
  }

  def __init__( self ):
    self.distance = 100
    self.mass = 1000
    self.height = 1
    self.width = 1
    self.ComputeSolidAngleFraction()

  def ComputeSolidAngleFraction( self ):
    surfaceArea = (self.distance*100)**2 * 4 * 3.1415927 #cm^2
    self.fluxFactor = 1./surfaceArea # cm^-2
    self.solidAngleFraction = (self.height * 100 * self.width * 100) / surfaceArea
    
