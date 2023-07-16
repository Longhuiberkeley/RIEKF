# RIEKF

- See [IEKF][pdf] to understand the theories
- See [RIEKF.m][main] to see the code;
    - use IMU to perform `prediction()`
    - use GPS's position to perform `correction`
        - use the adjoint to put it on the left, then we put it back to the right
    - use `odometry()` to combine wheel vehicle's non-holonomic constraint and encoder to perform correction
    - use `nonholonomic()` to correct the state by the non-holonomic assumption 

- not finished:
    - LIEKF.m 


[main]: filters/RIEKF.m
[pdf]: IEKF.pdf


This code is mostly modified from the GitHub repo, Invariant-ekf.