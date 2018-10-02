class Exposure(object):
	"""
	Exposure Class:

	Parameters:
	-----------

	exptime: float or int
	Exposure time in seconds
	"""

	def __init__(self):

		self.exptime = 900
		self.ra = 0.0
		self.dec = 0.0
		self.theta = 0.0
		self.airmass = 1.0

