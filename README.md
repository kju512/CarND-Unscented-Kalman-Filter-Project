# Unscented Kalman Filter Project
This is a Self-Driving Car Engineer Nanodegree Program completed by Michael chen.This program implements an unscented kalman filter.

---


## Build and Run

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../data/output.txt`


## Specification
You can use Radar or Lidar dataset only as you wish.The method is:
Open the ukf.cpp file, set "use_laser_ = true;" or "use_laser_ = false;" at the Init() function,you can switch the lidar dataset utilization  ability.As the same method, set "use_radar_ = true;" or "use_radar_ = false;" at the Init() function,you can switch the radar dataset utilization ability.

Under my parameters setting, 
if only Radar data is used,RMSE output will be [0.228003,0.527052,0.284155,0.434147].
if only lidar data is used,RMSE output will be [0.182224,0.161809,0.690417,0.316651].
if both data are used,RMSE output will be [0.0830559,0.128115,0.32019,0.271567].


## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

