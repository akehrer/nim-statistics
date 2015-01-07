import algorithm
import math

const
  NAN = 0.0/0.0 # floating point not a number (NaN)

## Some additional functions from math.h are needed that aren't included in the math module

proc cIsNaN(x: float): int {.importc: "isnan", header: "<math.h>".}
  ## returns non-zero if `x` is not a number

proc isNaN*(x: float): bool =
  ## converts the integer result from `cIsNaN` to a boolean
  if cIsNaN(x) == 1: true
  else: false

proc erf*(x: float): float {.importc: "erf", header: "<math.h>".}
  ## computes the error function (also called the Gauss error function)

proc erfc*(x: float): float {.importc: "erfc", header: "<math.h>".}
  ## computes the complementary error function (also called the Gauss error function)

## These are some additional descriptive statistics not found in the math module

proc median*(x: openArray[float]): float = 
  ## computes the median of the elements in `x`. 
  ## If `x` is empty, NaN is returned.
  
  var sx = @x # convert to a sequence since sort() won't take an openArray
  sx.sort(system.cmp[float])
  
  try:
    if sx.len mod 2 == 0:
      var n1 = sx[(sx.len - 1) div 2]
      var n2 = sx[sx.len div 2]
      result = (n1 + n2) / 2.0
    else:
      result = sx[(sx.len - 1) div 2]
  except IndexError:
    result = NAN

proc quantile*(x: openArray[float], frac: float): float = 
  ## computes the quantile value of `x` determined by the fraction `frac`
  ## `frac` must be between 0 and 1 so for the 25th quantile 
  ## the value should be 0.25, any other value returns `NaN`
  if x.len == 0:
    result = NAN  
  elif frac < 0.0 or frac > 1.0:
    result = NAN
  elif frac == 0.0:
    result = x.min
  else:
    var sx = @x # convert to a sequence since sort() won't take an openArray
    sx.sort(system.cmp[float])

    var n = sx.len - 1  # max index
    var i = int(math.floor(float(n) * frac))  # quantile index

    if i == n:
      result = sx[n]
    elif sx.len mod 2 == 0:
      # even length
      var n1 = sx[i]
      var n2 = sx[i+1]
      result = (n1 + n2) / 2.0
    else:
      # odd length
      result = sx[i]


## Gaussian (Normal) Distribution
## The functions are from http://en.wikipedia.org/wiki/Normal_distribution
type
  GaussDist* = object
    mu, sigma: float

proc NormDist*(): GaussDist = 
  ## A Normal Distribution is a special form of the Gaussian Distribution with
  ## mean 0.0 and standard deviation 1.0
  result.mu = 0.0
  result.sigma = 1.0

proc mean*(g: GaussDist): float =
  result = g.mu

proc standardDeviation*(g: GaussDist): float = 
  result = g.sigma

proc variance*(g: GaussDist): float = 
  result = math.pow(g.sigma, 2)

proc pdf*(g: GaussDist, x: float): float = 
  var numer, denom: float

  numer = math.exp(-(math.pow((x - g.mu), 2)/(2 * math.pow(g.sigma, 2))))
  denom = g.sigma * math.sqrt(2 * math.PI)
  result = numer / denom

proc cdf*(g: GaussDist, x: float): float = 
  var z: float

  z = (x - g.mu) / (g.sigma * math.sqrt(2))
  result = 0.5 * (1 + erf(z))

when isMainModule:
  # Test median()
  var med1: array[0..6, float]
  var med2: array[0..7, float]
  var med3 = newSeq[float]()
  var med4: array[1, float]
  var med5: array[2, float]
  med1 = [1.4, 3.6, 6.5, 9.3, 10.2, 15.1, 2.2]
  med2 = [1.4, 3.6, 6.5, 9.3, 10.2, 15.1, 2.2, 0.5]
  med4 = [2.3]
  med5 = [2.2, 2.5]

  assert(median(med1) == 6.5)
  assert(median(med2) == 5.05)
  assert(isNaN(median(med3)))  # test an empty sequence
  assert(median(med4) == 2.3)
  assert(median(med5) == 2.35)  

  # Test quantile()

  assert(quantile(med1, 0.5) == median(med1))
  assert(quantile(med2, 0.5) == median(med2))
  assert(quantile(med2, 0.75) == 9.75)
  assert(abs(quantile(med2, 0.9) - 12.65) < 1e-8)
  assert(abs(quantile(med2, 0.1) - 0.95) < 1e-8)
  assert(quantile(med2, 1.0) == med2.max)
  assert(quantile(med2, 0.0) == med2.min)
  
  # Test GaussDist
  var n = NormDist()
  var gnorm = GaussDist(mu: 0.0, sigma: 1.0)

  assert(n.mean == gnorm.mean)
  assert(n.standardDeviation == gnorm.standardDeviation)
  assert(n.variance == gnorm.variance)
  assert(abs(n.pdf(0.5) - 0.3520653267642995) < 1e-8)
  assert(abs(n.cdf(0.5) - 0.6914624612740131) < 1e-8)

  echo "SUCCESS: Tests passed!"
