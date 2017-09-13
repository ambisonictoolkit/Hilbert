The Hilbert Transform Library : Read Me
========================
_Various implementations of the Hilbert Transform and related functions._


The [Hilbert Transform](https://en.wikipedia.org/wiki/Hilbert_transform) is
offered in various forms:

* __HilbertW__ - via Weaver's Third Method.
* __HilbertH__ - via Weaver's Second Method, known as Hartley Phasing.
* __HilbertPDN__ - via Hartley Phasing, expressed as a 12th-order Phase
Differencing Network.

These
[pseudo-UGens](http://doc.sccode.org/Guides/WritingUGens.html#Pseudo-UGens)
return
[phase-quadrature](https://ccrma.stanford.edu/~jos/st/In_Phase_Quadrature_Sinusoidal.html)
outputs in a form suitable for use as an
[analytic signal](https://ccrma.stanford.edu/~jos/st/Analytic_Signals_Hilbert_Transform.html).
Each form offers separate pseudo-UGens for returning just the _real_ or
_imaginary_ output independently by appending __Re__ or __Im__, respectively,
to the class name. E.g., __HilbertWRe__ & __HilbertWIm__.

Additionally, each class includes further Hilbert related transforms and
analyses:
* Phase Rotation
* [Single-Sideband Modulation (SSB)](https://en.wikipedia.org/wiki/Single-sideband_modulation)
* Instantaneous Amplitude Analysis
* Instantaneous Phase Analysis


&nbsp;

&nbsp;

Installing
==========

Distributed via the
[SC3 Hilbert Transform Quark Library](https://github.com/ambisonictoolkit/Hilbert).

Start by reviewing the Quark installation instructions
[found here](https://github.com/supercollider-quarks/quarks#installing). See
also [Using Quarks](http://doc.sccode.org/Guides/UsingQuarks.html).

With [git](https://git-scm.com/) installed, you can easily install the
[SC3 Hilbert Transform Quark Library](https://github.com/ambisonictoolkit/Hilbert)
directly by running the following line of code in SuperCollider:

    Quarks.install("https://github.com/ambisonictoolkit/Hilbert.git");



Feedback and Bug Reports
========================

Known issues are logged at
[GitHub](https://github.com/ambisonictoolkit/Hilbert/issues).

&nbsp;


List of Changes
---------------



Version 0.1.0

* First Public Release.


&nbsp;

&nbsp;

Credits
=======

&nbsp;

Copyright the ATK Community, Joseph Anderson, and Michael McCrea, 2017.

* J Anderson : [[e-mail]](mailto:joanders[at]uw.edu)
* M McCrea : [[e-mail]](mailto:mtm5[at]uw.edu)

&nbsp;

The development of the Hilbert Transform Library for SuperCollider3 is
supported by
[DXARTS, Center for Digital Arts and Experimental Media](https://dxarts.washington.edu/).


Contributors
------------

Version 0.1.0
*  Joseph Anderson (@joslloand)
*  Michael McCrea (@mtmccrea)
