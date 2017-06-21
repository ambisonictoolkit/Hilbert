/*
Hilbert transform - Weaver quadrature method

Well-behaved phase response.

Uses PV_BrickWall and is therefore delayed by:
(fftsize - s.options.blockSize) / s.sampleRate

Real and imaginary components can be accessed independently using
HilbertWRe and HilbertWIm, respectively.
*/

HilbertW { |in, fftsize = 2048|
	*ar {
		^[
			HilbertWRe.ar(in, fftsize),
			HilbertWIm.ar(in, fftsize),
		]
	}
}

HilbertWRe { |in, fftsize = 2048|
	*ar {
		var del = (fftsize - BlockSize.ir) / SampleRate.ir;

		^DelayN.ar(in, del, del)
	}
}

HilbertWIm { |in, fftsize = 2048|
	*ar {
		var cos, sin, inputDbl, sinBr, cosBr;

		// quadrature oscillators at nyquist/2
		cos = SinOsc.ar(SampleRate.ir/4, pi/2);
		sin = SinOsc.ar(SampleRate.ir/4, 0);
		inputDbl = in * 2;

		// cosine and sine branches using brickwall for lowpass at nyquist/2
		cosBr = IFFT.ar(
			PV_BrickWall(
				FFT(
					LocalBuf(fftsize),
					inputDbl * cos
				),
				-0.5
			)
		);

		sinBr = IFFT.ar(
			PV_BrickWall(
				FFT(
					LocalBuf(fftsize),
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
(kernelsize - 1) / 2 + kernelsize - s.options.blockSize
*/

HilbertH {

}


/*
Hilbert transform - Phase Delay Network (PDN), IIR
*/

HilbertPDN {

}