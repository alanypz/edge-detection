# edge-detection

Implementation of canny and sobel edge detection algorithms.

## canny

From command line:
```
gcc canny.c
./a.out
```

Produces three .pgm image files for magnitude, peaks, and thresholds.

## sobel

From command line:
```
gcc sobel.c
./a.out
```

Command line output:
```
Setting sigma value to default: 1.000000
Setting threshold value to: 0.020000

High threshold: 174
Low threshold: 60
```

Produces three .pgm image files for magnitude, peaks, and thresholds.