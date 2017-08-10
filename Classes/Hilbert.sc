/*
	Copyright Joseph Anderson and Michael McCrea, 2017
		J Anderson	joanders@uw.edu
		M McCrea		mtm5@uw.edu

	This file is part of the Hilbert quark for SuperCollider 3 and is free software:
	you can redistribute it and/or modify it under the terms of the GNU General
	Public License as published by the Free Software Foundation, either version 3
	of the License, or (at your option) any later version.

	This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the GNU General Public License for more details.
	<http://www.gnu.org/licenses/>.
*/


//---------------------------------------------------------------------
//	Various implementations of the Hilbert transform.
//
// ::Classes::
// Weaver method:
// HilbertW
// HilbertWRe
// HilbertWIm
//
// Harley method:
// HilbertH
// HilbertHRe
// HilbertHIm
//
// Phase Difference Network (PDN) method:
// HilbertPDN
// HilbertPDNRe
// HilbertPDNIm
//
//---------------------------------------------------------------------



/*
Hilbert transform - Weaver quadrature method
Well-behaved phase response.
Uses PV_BrickWall and is therefore delayed by:
(size - s.options.blockSize) / s.sampleRate
Real and imaginary components can be accessed independently using
HilbertWRe and HilbertWIm, respectively.
*/

HilbertW {

	*ar { |in, size = 2048, mul = 1, add = 0.0|

		^[
			HilbertWRe.ar(in, size, mul, add),
			HilbertWIm.ar(in, size, mul, add),
		]

	}

	*arMag { |in, size = 2048, mul = 1, add = 0.0|
		var mag, real, imag;

		real = HilbertWRe.ar(in, size);
		imag = HilbertWIm.ar(in, size);

		mag = (real.squared + imag.squared).sqrt;

		^((mag * mul) + add)
	}

	*arPhase { |in, size = 2048, mul = 1, add = 0.0|
		var phase, real, imag;

		real = HilbertWRe.ar(in, size);
		imag = HilbertWIm.ar(in, size);

		phase = imag.atan2(real);

		^((phase * mul) + add)
	}

	*arRotate { |in, angle = 0.0, size = 2048, mul = 1, add = 0.0|
		var out, real, imag;

		real = HilbertWRe.ar(in, size);
		imag = HilbertWIm.ar(in, size);

		out = (angle.cos * real) + (angle.sin * imag);

		^((out * mul) + add)
	}

	*arSSB { |in, freq = 0.0, size = 2048, mul = 1, add = 0.0|
		var out, real, imag;
		var quadOsc, delay, phaseOffset;

		delay = size * SampleDur.ir - ControlDur.ir;
		phaseOffset = (-2pi * freq * delay).mod(2pi);

		quadOsc = SinOsc.ar(freq, [pi/2, 0] + phaseOffset);

		real = HilbertWRe.ar(in, size);
		imag = HilbertWIm.ar(in, size);

		out = (quadOsc.at(0) * real) - (quadOsc.at(1) * imag);

		^((out * mul) + add)
	}

}

HilbertWRe {

	*ar { |in, size = 2048, mul = 1, add = 0.0|
		var del = (size - BlockSize.ir) / SampleRate.ir;
		^DelayN.ar(in, del, del, mul, add)
	}
}

HilbertWIm {

	*ar { |in, size = 2048, mul = 1, add = 0.0|
		var nyqDiv2, cos, sin, inputDbl, sinBr, cosBr, out;

		// quadrature oscillators at nyquist/2
		nyqDiv2 = SampleRate.ir/4;
		cos = Impulse.ar(nyqDiv2);
		cos = Delay2.ar(cos, -1, cos);
		sin = Delay1.ar(cos);
		inputDbl = in * 2;

		// cosine and sine branches using brickwall for lowpass at nyquist/2
		cosBr = IFFT.ar(
			PV_BrickWall(
				FFT(
					LocalBuf(size),
					inputDbl * cos
				),
				-0.5
			)
		);

		sinBr = IFFT.ar(
			PV_BrickWall(
				FFT(
					LocalBuf(size),
					inputDbl * sin
				),
				-0.5
			)
		);

		// modulate and sum
		out = ([cosBr, sinBr] * [sin, cos.neg]).sum
		^((out * mul) + add)
	}
}




/*
Hilbert transform - Hartley method

Uses an FIR filter with Convolution2.ar and is therefore delayed by:
size/2 + size - s.options.blockSize

Kernel is rendered to sample delay.
*/

HilbertH {

	*ar { |in, size=2048, mul=1.0, add=0.0|
		var r, kernel_r, real;

		^[
			HilbertHRe.ar(in, size, mul, add),
			HilbertHIm.ar(in, size, mul, add)
		]

	}

	// calculate real coefficients as delayed impulse
	*calcRealCoeffs { |size|

		^HilbertHRe.calcCoeffs(size)
	}

	// calculate imag coefficients via (1-cos(t)) / t
	*calcImagCoeffs { |size|

		^HilbertHIm.calcCoeffs(size)
	}

	*arMag { |in, size = 2048, mul = 1, add = 0.0|
		var mag, real, imag;

		real = HilbertHRe.ar(in, size);
		imag = HilbertHIm.ar(in, size);

		mag = (real.squared + imag.squared).sqrt;

		^((mag * mul) + add)
	}

	*arPhase { |in, size = 2048, mul = 1, add = 0.0|
		var phase, real, imag;

		real = HilbertHRe.ar(in, size);
		imag = HilbertHIm.ar(in, size);

		phase = imag.atan2(real);

		^((phase * mul) + add)
	}

	*arRotate { |in, angle = 0.0, size = 2048, mul = 1, add = 0.0|
		var out, real, imag;

		real = HilbertHRe.ar(in, size);
		imag = HilbertHIm.ar(in, size);

		out = (angle.cos * real) + (angle.sin * imag);

		^((out * mul) + add)
	}

	*arSSB { |in, freq = 0.0, size = 2048, mul = 1, add = 0.0|
		var out, real, imag;
		var quadOsc, delay, phaseOffset;

		delay = (size + ((size/ 2).floor).asInt) * SampleDur.ir - ControlDur.ir;
		phaseOffset = (-2pi * freq * delay).mod(2pi);

		quadOsc = SinOsc.ar(freq, [pi/2, 0] + phaseOffset);

		real = HilbertHRe.ar(in, size);
		imag = HilbertHIm.ar(in, size);

		out = (quadOsc.at(0) * real) - (quadOsc.at(1) * imag);

		^((out * mul) + add)
	}
}

HilbertHRe {

	*ar { |in, size=2048, mul=1.0, add=0.0|
		var delay;

		delay = size/2 + size-BlockSize.ir;  // in samples
		delay = delay / SampleRate.ir;  // in seconds

		// delay aligned to sample, use DelayN
		^DelayN.ar(in, delay, delay, mul, add);
	}

	// if identical signal paths are preferred
	// between real and imaginary parts, at the expense of an
	// extra convolution
	*arConv { |in, size=2048, mul=1.0, add=0.0|
		var r, kernel_r;

		r = HilbertHRe.calcCoeffs(size);
		kernel_r = LocalBuf(size, 1).set(r);

		^Convolution2.ar(in, kernel_r, framesize: size, mul: mul, add: add);
	}

	// calculate real coefficients as delayed impulse
	*calcCoeffs { |size|
		var half_win, xReal;

		half_win = size/2;

		// real response
		xReal = Array.fill(size, { 0.0 });
		xReal.put(half_win, 1.0);

		^xReal
	}

	// // calculate real coefficients via sinc
	// *calcCoeffs { |size|
	// 	var half_win, window, xReal;
	//
	// 	half_win = size/2;
	// 	window = Signal.hanningWindow(size+1);
	//
	// 	// imaginary response
	// 	xReal = Array.series(size+1, half_win.neg, 1);
	// 	xReal = xReal.collect({|i| (i == 0).if({ 1 }, { sin(pi * i) / (pi * i) }) }) * window;
	//
	// 	^xReal.keep(size)
	// }
}


HilbertHIm {
	*ar {
		| in, size=2048, mul=1.0, add=0.0 |
		var i, kernel_i, image;

		i = HilbertHIm.calcCoeffs(size);
		kernel_i = LocalBuf(size, 1).set(i);

		^Convolution2.ar(in, kernel_i, framesize: size, mul: mul, add: add);
	}

	// calculate imag coefficients via (1-cos(t)) / t
	*calcCoeffs { |size|
		var half_win, window, xImag;

		half_win = size/2;
		window = Signal.hanningWindow(size+1);

		// imaginary response
		xImag = Array.series(size+1, half_win.neg, 1);
		xImag = xImag.collect({|i| (i == 0).if({ 0 }, { 1 - cos(pi * i) / (pi * i) }) }) * window;
		^xImag.keep(size);
	}
}


/*
Hilbert transform - Phase Difference Network (PDN), IIR
a 12 pole (6 per side) Hilbert IIR filter
based on Sean Costello and Bernie Hutchins
created by J Anderson - 7 jan 2001

See also:

B. Hutchins, "The Design of Wideband Analog 90° Phase Differencing Networks
without Large Spread of Capacitor Values", Electronotes, Special Issue G, No. 168
"http://electronotes.netfirms.com/EN168-90degreePDN.PDF"

B. Hutchins, "Calculation of the Poles of 90° Phase Difference Networks [after Weaver]"
in Musical Engineer's Handbook: Musical Engineering for Electronic Music. Ithaca, NY: Electronotes. 1975.
"http://electronotes.netfirms.com/MEHCh6aPart.PDF"

Weaver, D. “Design of RC Wide-Band 90-Degree Phase-Difference Network.”
Proceedings of the IRE, vol. 42, no. 4, 1954, pp. 671–676.
"http://ieeexplore.ieee.org/document/4051669/"

Other useful design references:

Ansari, R. “IIR Discrete-Time Hilbert Transformers.” Acoustics, Speech and Signal Processing,
IEEE Transactions On, vol. 35, no. 8, 1987, pp. 1116–1119.
"http://ieeexplore.ieee.org/abstract/document/1165250/"

Schüssler, H., and W. Steffen. “Halfband Filters and Hilbert Transformers.” Circuits, Systems and
Signal Processing, vol. 17, no. 2, 1998, pp. 137–164.
"https://link.springer.com/article/10.1007/BF01202851"


*/

HilbertPDN {

	*ar { |in, mul = 1.0, add = 0.0|
		^[
			HilbertPDNRe.ar(in, mul, add),
			HilbertPDNIm.ar(in, mul, add)
		]
	}

	*calcSOSCoefs { |poles|
		var gammas, coefs, b1, b2;

		gammas = (15.0 * pi / SampleRate.ir) * poles;

		coefs = [];
		gammas.do({ arg gamma;
			coefs = coefs.add((gamma-1)/(gamma+1))
		});

		// 1st order allpass filters coefs are grouped into coefs for 2nd order sections
		b1 = [];
		b2 = [];
		// gathers even and odd
		poles.size.div(2).do({ arg i;
			b1 = b1.add(coefs[2*i] + coefs[(2*i)+1]);
			b2 = b2.add(coefs[2*i] * coefs[(2*i)+1]);
		});

		^[b1, b2]
	}

	*arMag { |in, mul = 1, add = 0.0|
		var mag, real, imag;

		real = HilbertPDNRe.ar(in);
		imag = HilbertPDNIm.ar(in);

		mag = (real.squared + imag.squared).sqrt;

		^((mag * mul) + add)
	}

	*arPhase { |in, mul = 1, add = 0.0|
		var phase, real, imag;

		real = HilbertPDNRe.ar(in);
		imag = HilbertPDNIm.ar(in);

		phase = imag.atan2(real);

		^((phase * mul) + add)
	}

	*arRotate { |in, angle = 0.0, mul = 1, add = 0.0|
		var out, real, imag;

		real = HilbertPDNRe.ar(in);
		imag = HilbertPDNIm.ar(in);

		out = (angle.cos * real) + (angle.sin * imag);

		^((out * mul) + add)
	}

	*arSSB { |in, freq = 0.0, mul = 1, add = 0.0|
		var out, real, imag;
		var quadOsc;

		quadOsc = SinOsc.ar(freq, [pi/2, 0]);  // non-linear phase is uncompensated

		real = HilbertPDNRe.ar(in);
		imag = HilbertPDNIm.ar(in);

		out = (quadOsc.at(0) * real) - (quadOsc.at(1) * imag);

		^((out * mul) + add)
	}
}

HilbertPDNRe {
	*ar { |in, mul = 1.0, add = 0.0|
		var b1, b2, hilbertCos;
		var poles;

		// values taken from Bernie Hutchins, "Musical Engineer's Handbook"
		// also found in Electronotes #43
		// pole values are grouped in order as described by Hutchins
		// for optimal realization of the second order section coefficients
		poles = [1.2524, 2770.1114, 5.5671, 364.7914, 22.3423, 89.6271];  // cos: real

		#b1, b2 = HilbertPDN.calcSOSCoefs(poles);

		hilbertCos = in;
		3.do({ |i|
			hilbertCos = SOS.ar(hilbertCos, b2[i], b1[i], 1.0, b1[i].neg, b2[i].neg)
		});

		^((hilbertCos * mul) + add)
	}
}

HilbertPDNIm {
	*ar { |in, mul = 1.0, add = 0.0|
		var b1, b2, hilbertSin;
		var poles;

		// values taken from Bernie Hutchins, "Musical Engineer's Handbook"
		// also found in Electronotes #43
		// pole values are grouped in order as described by Hutchins
		// for optimal realization of the second order section coefficients
		poles = [0.3609, 798.4578, 2.7412, 179.6242, 11.1573, 44.7581];  // sin: imag

		#b1, b2 = HilbertPDN.calcSOSCoefs(poles);

		hilbertSin = in;
		3.do({ |i|
			hilbertSin = SOS.ar(hilbertSin, b2[i], b1[i], 1.0, b1[i].neg, b2[i].neg)
		});

		^((hilbertSin * mul) + add)
	}
}

