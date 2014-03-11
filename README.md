helm
====

![Controller block diagram](helm.png)

This PID controller, implemented wholly within [helm.h], features:
 * low pass filtering of the process derivative,
 * windup protection,
 * automatic reset on actuator saturation,
 * anti-kick on setpoint change using "derivative on measurement",
 * incremental output for bumpless manual-to-automatic transitions,
 * a unified controller gain parameter,
 * exposure of all independent physical time scales, and
 * the ability to accommodate varying sample rate.

References
----------

-- Åström, Karl J. and Richard M. Murray.  [Feedback systems: an introduction for scientists and engineers](http://www.worldcat.org/isbn/9781400828739).  Princeton University Press, April 2008.  [Available online](http://www.cds.caltech.edu/~murray/amwiki/index.php/Main_Page).
