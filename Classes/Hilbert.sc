/*
	Copyright Joseph Anderson and Michael McCrea, 2017
		J Anderson	j.anderson[at]ambisonictoolkit.net
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
// Phase Delay Network (PDN) method:
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

	*ar { |in, size = 2048|

		^[
			HilbertWRe.ar(in, size),
			HilbertWIm.ar(in, size),
		]

	}
}

HilbertWRe {

	*ar { |in, size = 2048|
		var del = (size - BlockSize.ir) / SampleRate.ir;
		del.poll;
		^DelayN.ar(in, del, del)
	}
}

HilbertWIm {

	*ar { |in, size = 2048|
		var cos, sin, inputDbl, sinBr, cosBr;

		// quadrature oscillators at nyquist/2
		cos = SinOsc.ar(SampleRate.ir/4, pi/2);
		sin = SinOsc.ar(SampleRate.ir/4, 0);
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
		^([cosBr, sinBr] * [sin, cos.neg]).sum
	}
}




/*
Hilbert transform - Hartley method

Uses an FIR filter with Convolution2.ar and is therefore delayed by:
(size - 1) / 2 + size - s.options.blockSize
*/

HilbertH {

	// // all-in-one
	// *ar { arg in, size=2048, mul=1.0, add=0.0;
	// 	var kernel_r, kernel_i;
	// 	var hilbertCoeffs, r, i, real, imag;
	//
	// 	hilbertCoeffs =  { |size|
	// 		var xReal, xImag, reflect, window, half_win, arr;
	//
	// 		half_win = (size-1)/2;
	// 		reflect = [1.0, -1.0].dup(size/2).flat;
	// 		window = Signal.hanningWindow(size);
	//
	// 		// real response
	// 		xReal = Array.series(size, half_win.neg, 1);
	// 		xReal = xReal.collect({|i| (i == 0).if({ 1 }, { sin(pi * i) / (pi * i) }) }) * window;
	//
	// 		// imaginary response
	// 		xImag = xReal * reflect;
	// 		[xReal, xImag]
	// 	};
	//
	// 	#r, i = hilbertCoeffs.(size);
	//
	// 	kernel_r = LocalBuf(size, 1).set(r);
	// 	kernel_i = LocalBuf(size, 1).set(i);
	//
	// 	#real, imag = Convolution2.ar(in, [kernel_r, kernel_i], framesize: size, mul: mul, add: add);
	//
	// 	^[ real, imag ]
	// }

	// alternate
	*ar { |in, size=2048, mul=1.0, add=0.0|
		var r, kernel_r, real;

		^[
			HilbertHRe.ar(in, size, mul, add),
			HilbertHIm.ar(in, size, mul, add)
		]

	}

	*calcRealCoeffs { |size|
		var half_win, window, xReal;

		half_win = (size-1)/2;
		window = Signal.hanningWindow(size);

		// real response
		xReal = Array.series(size, half_win.neg, 1);
		^xReal.collect({|i| (i == 0).if({ 1 }, { sin(pi * i) / (pi * i) }) }) * window;
	}

	*calcImagCoeffs { |size, realCoeffs|
		var reflect;

		realCoeffs ?? {Error("realCoeffs not provided. You can generate them with HilbertH.calcRealCoeffs").throw};

		reflect = [1.0, -1.0].dup(size/2).flat;
		^realCoeffs * reflect;
	}
}

HilbertHRe {

	*ar { |in, size=2048, mul=1.0, add=0.0|
		var delay;

		delay = (size-1)/2 + size-BlockSize.ir;
		delay = delay / SampleRate.ir;

		// delay will be exactly between samples, so DelayL is OK
		// compared to generating the real component with
		// the convolution method, the difference in signal is ~120dB
		^DelayL.ar(in, delay, delay, mul, add);
	}

	// if identical signal paths are preferred
	// between real and imaginary parts, at the expense of an
	// extra convolution
	*arConv { |in, size=2048, mul=1.0, add=0.0|
		var r, kernel_r;

		r = HilbertH.calcRealCoeffs(size);
		kernel_r = LocalBuf(size, 1).set(r);

		^Convolution2.ar(in, kernel_r, framesize: size, mul: mul, add: add);
	}

}


HilbertHIm : HilbertH {
	*ar {
		| in, size=2048, mul=1.0, add=0.0 |
		var r, i, kernel_i, image;

		r = HilbertH.calcRealCoeffs(size);
		"realz: ".post; r.postln;
		i = HilbertH.calcImagCoeffs(size, r);
		kernel_i = LocalBuf(size, 1).set(i);

		^Convolution2.ar(in, kernel_i, framesize: size, mul: mul, add: add);
	}
}


/*
Hilbert transform - Phase Delay Network (PDN), IIR
*/

HilbertPDN {

}