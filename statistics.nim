import algorithm
import math

const
  NAN = 0.0/0.0 # floating point not a number (NaN)

## Some additional functions from math.h are needed that aren't included in the math module

proc cIsNaN(x: float): int {.importc: "isnan", header: "<math.h>".}
  ## Returns non-zero if `x` is not a number


proc isNaN*(x: float): bool =
  ## Converts the integer result from `cIsNaN` to a boolean
  if cIsNaN(x) != 0: true
  else: false


proc erf*(x: float): float {.importc: "erf", header: "<math.h>".}
  ## Computes the error function (also called the Gauss error function)


proc erfc*(x: float): float {.importc: "erfc", header: "<math.h>".}
  ## Computes the complementary error function (also called the Gauss error function)


## These are some additional descriptive statistics not found in the math module

proc standardDeviation*(x: openArray[float]): float =
  ## Computes the standard deviation of `x`
  result = math.sqrt(math.variance(x))


proc unbiasedVariance*(x: openArray[float]): float =
  ## Computes the unbiased estimate sample variance of `x`
  ## If the length of `x` is lest than 2, NaN is returned.
  result = 0.0
  var n = x.len
  var xbar = math.mean(x)
  var s2: float
  for i in x:
    s2 += math.pow((i - xbar), 2)
  result = 1/(n-1) * s2
  

proc median*(x: openArray[float]): float = 
  ## Computes the median of the elements in `x`. 
  ## If `x` is empty, NaN is returned.
  if x.len == 0:
    return NAN
  
  var sx = @x # convert to a sequence since sort() won't take an openArray
  sx.sort(system.cmp[float])
  
  if sx.len mod 2 == 0:
    var n1 = sx[(sx.len - 1) div 2]
    var n2 = sx[sx.len div 2]
    result = (n1 + n2) / 2.0
  else:
    result = sx[(sx.len - 1) div 2]


proc quantile*(x: openArray[float], frac: float): float = 
  ## Computes the quantile value of `x` determined by the fraction `frac`
  ## If `x` is empty, NaN is returned.
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


proc skewness*(x: openArray[float]): float = 
  ## Computes the skewness of `x` as the adjusted Fisher-Pearson 
  ## standardized moment coefficient.
  ## If the length of `x` is lest than 3, NaN is returned.
  if x.len < 3:
    return NAN

  var xbar = math.mean(x)
  var n = float(x.len)

  var lhs = math.pow(n, 2) / ((n - 1.0) * (n - 2.0))
  var m3: float
  var s3: float

  for i in x:
    m3 += math.pow((i - xbar), 3)
    s3 += math.pow((i - xbar), 2)

  m3 *= 1/n
  s3 *= 1/(n-1)
  s3 = math.pow(s3, 3/2)

  result = lhs * m3 / s3


proc kurtosis*(x: openArray[float]): float = 
  ## Computes the population excess kurtosis using sample `x`.
  ## If the length of `x` is lest than 4, NaN is returned.
  if x.len < 4:
    return NAN

  var xbar = math.mean(x)
  var s = unbiasedVariance(x)
  var n = float(x.len)

  var lhs = ((n+1) * n) / ((n-1) * (n-2) * (n-3))
  var rhs = 3 * (math.pow((n-1), 2) / ((n-2) * (n-3)))

  var cen: float

  for i in x:
    cen += math.pow((i - xbar), 4)

  cen *= 1/math.pow(s, 2)

  result =  lhs * cen - rhs


## Gaussian (Normal) Distribution
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


proc median*(g: GaussDist): float = 
  result = g.mu


proc standardDeviation*(g: GaussDist): float = 
  result = g.sigma


proc variance*(g: GaussDist): float = 
  result = math.pow(g.sigma, 2)


proc skewness*(g: GaussDist): float = 
  result = 0.0


proc kurtosis*(g: GaussDist): float = 
  result = 0.0


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
  # Setup some test data
  var data1: array[0..6, float]
  var data2: array[0..7, float]
  var data3 = newSeq[float]()
  var data4: array[1, float]
  var data5: array[2, float]
  data1 = [1.4, 3.6, 6.5, 9.3, 10.2, 15.1, 2.2]
  data2 = [1.4, 3.6, 6.5, 9.3, 10.2, 15.1, 2.2, 0.5]
  data4 = [2.3]
  data5 = [2.2, 2.5]
  
  # Test median()
  assert(median(data1) == 6.5)
  assert(median(data2) == 5.05)
  assert(isNaN(median(data3)))  # test an empty sequence
  assert(median(data4) == 2.3)
  assert(median(data5) == 2.35)  

  # Test quantile()
  assert(quantile(data1, 0.5) == median(data1))
  assert(quantile(data2, 0.5) == median(data2))
  assert(quantile(data2, 0.75) == 9.75)
  assert(abs(quantile(data2, 0.9) - 12.65) < 1e-8)
  assert(abs(quantile(data2, 0.1) - 0.95) < 1e-8)
  assert(quantile(data2, 1.0) == data2.max)
  assert(quantile(data2, 0.0) == data2.min)

  # Test skewness()
  assert(abs(skewness(data1) - 0.5658573331608636) < 1e-8)
  assert(abs(skewness(data2) - 0.678154254246652) < 1e-8)
  assert(isNaN(skewness(data3)))
  assert(isNaN(skewness(data4)))
  assert(isNaN(skewness(data5)))

  # Test kurtosis()
  assert(abs(kurtosis(data1)- -0.6022513181999382) < 1e-8)
  assert(abs(kurtosis(data2)- -0.5461679512306573) < 1e-8)
  assert(isNaN(kurtosis(data3)))
  assert(isNaN(kurtosis(data4)))
  assert(isNaN(kurtosis(data5)))
  
  # Test GaussDist
  var n = NormDist()
  var gnorm = GaussDist(mu: 0.0, sigma: 1.0)

  assert(n.mean == gnorm.mean)
  assert(n.median == gnorm.median)
  assert(n.standardDeviation == gnorm.standardDeviation)
  assert(n.variance == gnorm.variance)
  assert(abs(n.pdf(0.5) - 0.3520653267642995) < 1e-8)
  assert(abs(n.cdf(0.5) - 0.6914624612740131) < 1e-8)

  echo "SUCCESS: Tests passed!"
