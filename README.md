# EventShape
Event shape calculator for MadAnalysis 5

Currently available functionalities:
  * Aplanarity
  * Sphericity
  * Transverse Sphericity
  * Thrust
  * Thrust Minor
  * Thrust Major
  * Thrust Axes

Simply add to the analysis as follows

In `user.h`:
```
private:
  EventShape * shape;
```

In `user.cpp`:
  * Initialization:
  ```
  shape = new EventShape();
  shape.Reset();
  ```
  
  * Execute:
  ```
  shape->calculateSphericity(my_jets);
  shape->calculateThrust(my_jets);
  ```
