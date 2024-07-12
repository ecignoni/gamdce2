# GAMDCE2

A simple and fast implementation of 2nd order Cumulant Expansion reweighting for GaMD simulations.

### Toy example in 1D

```python
n_samples = 10000000
# collective variable
cv = np.random.normal(size=n_samples)
# bias (you should read it from the gamd log)
bias = np.sin(cv) + np.random.normal(size=n_samples) * 0.1
# number of grid bins
n_bins = 500
# minimum and maximum value of the fes grid
fes_min = np.min(cv)
fes_max = np.max(cv)
delta = ( fes_max - fes_min ) / 5
fes_min = fes_min - delta
fes_max = fes_max + delta
# thermodynamic beta (should have same units as bias)
beta = 1.0

fes = gamdce2.reweight_1d(n_samples, cv, bias, n_bins, fes_min, fes_max, beta)
```

### Toy example in 2D

```python
n_samples = 10000000
# collective variables
cv1 = np.random.normal(size=n_samples)
cv2 = np.random.normal(size=n_samples) * 0.5
# bias (you should read it from the gamd log)
bias = np.sin(cv1) + np.random.normal(size=n_samples) * 0.1
# number of grid bins for both cvs
n_bins1, n_bins2 = 100, 50

def get_min_max_fes(cv, delta=5):
    fes_min = np.min(cv)
    fes_max = np.max(cv)
    delta = ( fes_max - fes_min ) / delta
    fes_min = fes_min - delta
    fes_max = fes_max + delta
    return fes_min, fes_max

# minimum and maximum value of the fes grid for both cvs
fes_min1, fes_max1 = get_min_max_fes(cv1)
fes_min2, fes_max2 = get_min_max_fes(cv2)

# thermodynamic beta (should have same units as bias)
beta = 1.0

fes = gamdce2.reweight_2d(n_samples, cv1, cv2, bias, n_bins1, n_bins2, fes_min1, fes_max1, fes_min2, fes_max2, beta)
```

### Helper to read AMBER GaMD log

Returns a pandas DataFrame with the log values:

```python
df = gamdce2.read_gamd_log("PATH_TO_LOG")
```
