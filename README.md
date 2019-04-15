This repository contains host software (Linux) for hackrf-hunt, a project to
produce a hackrf sweeping software by db level.  This can be used with a directional antenna to foxhunt known (e.g. lost drone) or unknown signals (wifi interference or targeting) 

I've been using a Vivaldi TSA-900 Ultra-Wideband (UWB) directional antenna with my hackrf one for successful signal hunting.

This code is preliminary. Use at your own risk.
                                  
## How to build the hackrf-hunt software on Linux:

### Prerequisites for Linux (Arch):
```
sudo pacman -S base-devel cmake libusb hackrf fftw
```

### Prerequisites for Linux (Debian/Ubuntu):
```
`sudo apt-get install build-essential hackrf cmake libusb-1.0-0-dev pkg-config libfftw3-dev`
```

### Build hackrf-hunt software on Linux:
```
mkdir build
cd build
cmake ..
make
sudo make install
```


## Clean CMake temporary files/dirs:
```
cd build
rm -rf *
```
This code was taken from the hackrf sweep code, so if you are familiar with the execution of hackrf_sweep this should be a breeze for you.  

Connect your hackrf and antenna (preferibly a directional antenna and to start, run the following:

```
$ hackrf-hunt -z "-20:-50"
call hackrf_sample_rate_set(20.000 MHz)
call hackrf_baseband_filter_bandwidth_set(15.000 MHz)
Sweeping from 0 MHz to 6000 MHz searching min DB: -50.000000 and max DB: -20.000000
Stop with Ctrl-C
1555364667, 65000000, 70000000, -43.94
1555364667, 85000000, 90000000, -48.63
1555364667, 105000000, 110000000, -23.71
1555364667, 205000000, 210000000, -48.25
1555364667, 520000000, 525000000, -48.24
1555364667, 565000000, 570000000, -45.23
1555364667, 65000000, 70000000, -40.76
^CCaught signal 2
```

From the above command and output you are receiving a sweep of the 1mhz-6Ghz spectrum for frequencies that fall between -50 on the lower end and -20 on the higher end.

You can specify a frequency range as well with something like the following command:

```
$ hackrf-hunt -z "-20:-50" -f 650:800
call hackrf_sample_rate_set(20.000 MHz)
call hackrf_baseband_filter_bandwidth_set(15.000 MHz)
Sweeping from 650 MHz to 810 MHz searching min DB: -50.000000 and max DB: -20.000000
Stop with Ctrl-C
1555364896, 735000000, 740000000, -41.65
1555364896, 735000000, 740000000, -44.83
1555364896, 735000000, 740000000, -44.35
1555364896, 735000000, 740000000, -46.18
1555364896, 735000000, 740000000, -43.78
1555364896, 735000000, 740000000, -38.40
1555364896, 735000000, 740000000, -41.26
^CCaught signal 2
```

From the above commend and output you will see that the output is familiar to the first example I provided, however you are now limiting the range of the scan to frequencies between 650mhz and 800mhz only.

The epoch time can be converted to recognizable time with utilities that can be found.  You can take this output and dump it into sqllite and then graph it using something like matplot.
