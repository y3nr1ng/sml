# Single molecule localization
Project description.


## Getting Started
The following instructions will get a copy of this project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy it on a live system.

### Prerequisites
You will need a functional Python environment, [conda](https://conda.io/docs/) is the preferred package management system for us.

### Installing
We first create a functional environment based on the [environments.yml](environments.yml) file.
```
conda env create -f environments.yml
```
To keep the environment up-to-date for any upstream bugs, package versions are not strictly fixed in the environment definition.

The newly created environment is now called `sml-dev`, activate it.
```
conda activate sml-dev
```

Next, we need to install a local copy of it, we are going to use `conda-develop`, for `pip` users, ignore this step.
```
conda install conda-build
conda develop .
```
This will create a `.pth` file in current environment, linking this repository to make it appears as installed.


## Running the tests
Explain how to run the automated tests.

To test whether everything goes as expected, go to [tests](tests) folder, and
```
python -m sml.cli.pix input.yml
```
This requests Python to use `sml.cli.pix` module, and exexcute it with command line arguments `input.yml`, which points to the file under [tests](tests). For now (2018/08/01), it should print parsed results.
```
{'FIT_MAX_DW': 100,
 'FIT_MAX_DX': 100,
 'FIT_MAX_DY': 100,
 'I_THRES_MAX': 450,
 'I_THRES_MIN': 50,
 'MAX_DI_I': 0.5,
 'MIN_SN_RATIO': 10,
 'NFSEP_IDP': 3,
 'NM_PIXEL_X': 140.8451,
 'NM_PIXEL_Y': 134.2282,
 'SOLVER_DZ': 1,
 'SOLVER_Z1': -1000,
 'SOLVER_Z2': 1000,
 'X_FIND_PIXELS': 5,
 'Y_FIND_PIXELS': 5,
 'frames': {'end': 10, 'start': 1},
 'i_photon': 0.1,
 'image_format': 'TIFF',
 'input': {'calibration_file': 'data/image.cab',
           'image_file': 'data/image.tif'},
 'mode': 'SML',
 'n_dim': '2D',
 'output': {'candidate': 'candidate',
            'fsum': 'Fsum',
            'spot': 'spot',
            'spoth': 'spot_high',
            'statistics': 'Fstats'},
 'roi': {'x1': 50, 'x2': 229, 'y1': 50, 'y2': 229},
 'scan_algorithm': 'regional_max',
 'verbose': True}
```

## Deployment
Additional notes about how to deploy this project on a live system.


## Built with


## Contributing


## Versioning


## Authors


## License
This project is licensed under the Apache 2.0 License - see the [LICENSE](LICENSE) file for details.


## Acknowledgments
Hat tip to anyone whose code was used/inspired.

---

## pix
### Description
This is the main code for single molecule localization. It reads a set
of frames of input images (in TIFF format), using differential frame
algorithm to find all the spots, and generate the output files. This
code is parallelized with OpenMP and MPI. With the fitting parameters
of the calibration curve, this code can also do 3D single molecule
localization.

### Input file
```
  1                    ! image file format: 1:tiff, 2:raw
  data/image.tif       ! image filename.
  data/image.cab       ! calibration parameter filename.
  spot                 ! output data filename for normal spots.
  spotH                ! output data filename for high intensity spots.
  Fsts                 ! output data filename for frame event statistics.
  Fsum                 ! output data filename for frame sums.
  candidate            ! output data filename for candidate events (pixels).
  0                    ! mode: 0:2D, 1:3D
  1 10                 ! frameID1, frameID2.
  50 229 50 229        ! frame_x1, frame_x2, frame_y1, frame_y2.
  5 5                  ! x_find_pixels, y_find_pixels.
  50 450               ! threshold1, threshold2 (in # of photons)
  3                    ! max. frame separation for identical particles
  140.8451 134.2282    ! nm per pixels in (x,y) directions
  0.1                  ! factor to convert pixel intensity to # of photons
  0                    ! rmode: 0:spot fitting, 1:fit & output spot
  1                    ! alg: 0:regional max, 1:max intensity
  100 100 100          ! max_dx, max_dy, max_dw (nm)
  10                   ! min_SN (Signal/Noise ratio)
  0.5                  ! max ratio of d(Intensity)/Intensity
  -1000 1000 1         ! Solve z-coord.: [z1,z2],dz
  0                    ! verbose message
```


## calbfit:
### Description:
Use the calibration data file (i.e., the length of X-axis or Y-axis
of ellipse shaped spot versus its z position) to fit to the function
form of the calibration curve. The fitting function is:

f(z) = w0*sqrt(1+((z-c)/d)**2*(1+A*((z-c)/d)+B*((z-c)/d)**2))

where "z" is the z position of the spot, and w0, A, B, c, d are the
fitting parameters.

### Usage:
```
./calbfit [-v] <calb_rawX_file> <calb_rawY_file> <outfn>
```
where <calb_rawX_file> and <calb_rawY_file> are the list of ellipical
axis length (for X and Y, respectively) and spot z position data, which
are measured in experiments. The data format is:
```
# z position  axis length
0.00000000    1658.56967455
0.99945975    1656.37265721
1.99891950    1654.17752838
2.99837925    1651.98428806
....          ....
```

### Output
The output file contains the fitting results of the calibration curve,
for both X-axis and Y-axis. For example:

```
w0x =  3.7951032353E+02 +- 1.9727E-10
WxA = -2.4912708613E-12 +- 2.4157E-12
WxB =  2.5000000000E-01 +- 1.7360E-12
Wxc =  1.1632308250E+03 +- 4.7015E-10
Wxd =  4.4804058223E+02 +- 5.9956E-10
w0y =  4.1035075257E+02 +- 1.6612E-10
WyA =  2.9161639768E-12 +- 3.4227E-12
WyB =  2.5000000000E-01 +- 2.1385E-12
Wyc =  6.3270737668E+02 +- 5.9411E-10
Wyd =  5.6715914417E+02 +- 7.2582E-10
```

## aimg
### Description
Use random number generator to generate an artificial image with spots,
which could be used for test. The image size, number of spots, and spot
size (in pixels) are fixed. The spot positions, spot intensity, spot
width, and the noise intensity are randomly generated.

### Input file
```
  testimg      ! output image filename (for both JPG and RAW)
  32339        ! random number seed
  200 200      ! iH, iW: image size (pixel)
  0.5 0.9      ! Signal intensity range: [0.0-1.0]
  0.0 0.1      ! Noise intensity range:  [0.0-1.0]
  0.5 2.5      ! Spot width range (pixel)
  20           ! number of spots
  9            ! spot size (pixels)
  10           ! mesh size of each pixel to generate spots
```

### Output files
- testimg.tab:
  List of randomly generated positions, width, intensity, and S/N ratio
  of each spot.

- testimg.jpg:
  The artificially generated image with spots.

- testimg.txt:
  The list of (x,y) position and intensity of each pixel of the JPEG image.


## pclst
### Description
This code finds the cluster of spots from the data file of spot list,
compute the correlation of spot positions of this cluster, and output
the cluster information to "xcor.dat".

### Input file
```
  0000_spot.txt           ! input spot data filename
  0000_Fsts.txt           ! input frame statistics filename
  0     0                 ! left-bottom corner of the image (nm)
  20000 20000             ! right-top   corner of the image (nm)
  0.3                     ! max. event density for starting frame
  2000                    ! max. distances of spots to form a cluster (nm)
  100                     ! min # of spots within a cluster to draw
  20                      ! # of rectangles to draw
  150                     ! # of clusters to keep in the first stage
  100                     ! n_intvls for histogram for xcor
  0                       ! verbose message output
```

### Output file
    xcor.dat:
    - N_REC:   number of clusters found.
    - REC_ID:  cluster ID
    - REC_C1:  coordinate of left-bottom corner of the cluster.
    - REC_C2:  coordinate of right-top corner of the cluster.
    - REC_NP:  number of spots in this cluster.
    - REC_NH:  number of histogram intervals of spot position correlation.
    - REC_F0:  the starting frame of input images for cluster analysis.
    - REC_DR:  interval size of the spot position correlation histogram.
    - REC_SP:  list of spot IDs of this cluster
    - list the histogram of spot position correlation.


## pspot
### Description
This code graphically displays the found spots on the screen, shows the
the clusters of spots from data file "xcor.dat", and generate the output
image and cluster information. This code is developed by OpenGL. It
generates 3 kinds of outputs:

- 0:spot mode:
  Display the found spots from data files and show the clusters of spots
  from data file "xcor.dat". In this mode user can interactively circle
  new clusters or delete the found clusters by mouse, in case that some
  of the clusters may not be determined correctly.

- 1:spotmesh mode:
  Defining a grid on the image (each mesh has size N(nm) x N(nm)), count
  the number of spots (events) within each mesh. Then use the spots count
  of each mesh as intensity to generate the spot image.

- 2:pixel mode:
  Generate the image of sum of all source (input) frames of images.


### Usage
```
./pspot [-v] <input>
 -v: display the image on the screen.
```
Without the `-v` flag, the code just process data, generate the output
images and data, and then stop. With the `-v` flag, it shows the image
in the screen. For 0:spot mode, user can use mouse to manually select
or delete the spot clusters. The final output images will be generate
when quit the code.

### Input file:
```
  0                   ! data format: 0:spot, 1:spotmesh, 2:pixel
  spot.txt            ! input data (spot list, or pixel sum (for mode 2))
  Fsts.txt            ! frame statistics filename
  spotmesh.txt        ! output spot meshed filename
  0.3                 ! max. event density for starting frame
  0                   ! image size: 0:full, 1:scaled
  512 512             ! image size  (pixel)
  20000 20000         ! image range (nm)
  100                 ! n_intvls for histogram for xcor
  20                  ! mesh size of the grid (nm)
  0.5                 ! gamma of the gray scale image
  /path/to/font72.glf ! the path of the font for display
```

### Output file:
- For mode 0: spot:
      Generate image of spots with identified clusters, and update "xcor.dat"
      if the user has manually added or deleted some clusters.

- For mode 1: spotmesh:
      Generate image of spots in each mesh of a given grid, where the
      spot intensity stands for the count of spots in that mesh.

- For mode 2: pixel:
      Generate image of sum of pixel values of the input image frames.
