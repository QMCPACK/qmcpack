
import math
import numpy as np

class averager:
  """Compute average, variance, and standard error for a set of data"""
  def __init__(self):
    self.sum = 0.0
    self.sum_sq = 0.0
    self.norm = 0

  def add_value(self,v):
    self.sum += v
    self.sum_sq += v*v
    self.norm += 1

  def average(self):
    if (self.norm == 0):
      return 0.0
    else:
      return self.sum/self.norm

  def variance(self):
    var = 0.0
    if (self.norm != 0):
      var = (self.sum_sq - self.sum*self.sum/self.norm)/self.norm
    return var

  def std_dev(self):
    return math.sqrt(self.variance())

  def error(self):
    err = 0.0
    if (self.norm > 1):
      var = self.variance()
      err = np.sqrt(var/(self.norm-1))
    return err

  def merge(self, other):
      self.sum += other.sum
      self.sum_sq += other.sum_sq
      self.norm += other.norm

