# EventShape
Event shape calculator for MadAnalysis 5

Currently available functionalities:
  * [Aplanarity](https://www.sciencedirect.com/science/article/pii/0370269378900618)
  * [Sphericity](https://www.sciencedirect.com/science/article/pii/0370269378900618)
  * [Transverse Sphericity](https://www.sciencedirect.com/science/article/pii/0370269378900618)
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
  
  ## TODO:
  - [ ] Add broadening
