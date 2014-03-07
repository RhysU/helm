helm
====

A header-only C99 proportional-integral-derivative (PID) controller.

The controller, implemented wholly in [helm.h](helm.h), features:
 * low pass filtering of the process derivative,
 * windup protection,
 * automatic reset on actuator saturation,
 * anti-kick on setpoint change using "derivative on measurement",
 * incremental output for bumpless manual-to-automatic transitions,
 * a unified controller gain parameter,
 * exposure of all independent physical time scales, and
 * the ability to accommodate varying sample rate.

<center><img src="figures/helm.png" width="95%"/></center>

The design and nomenclature is based largely on Figure 10.11 of
[Feedback Systems](http://www.worldcat.org/isbn/9781400828739) by
Astrom and Murray.

References
----------

TODO
