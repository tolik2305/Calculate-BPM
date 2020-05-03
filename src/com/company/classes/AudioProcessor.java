package com.company.classes;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;

import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;


/** Audio processing class (adapted from PerformanceMatcher). */
public class AudioProcessor {

    /** Input data stream for this performance (possibly in compressed format) */
    protected AudioInputStream rawInputStream;

    /** Uncompressed version of <code>rawInputStream</code>.
     *  In the (normal) case where the input is already PCM data,
     *  <code>rawInputStream == pcmInputStream</code> */
    protected AudioInputStream pcmInputStream;

    /** Line for audio output (not used, since output is done by AudioPlayer) */
//	protected SourceDataLine audioOut;

    /** Format of the audio data in <code>pcmInputStream</code> */
    protected AudioFormat audioFormat;

    /** Number of channels of audio in <code>audioFormat</code> */
    protected int channels;

    /** Sample rate of audio in <code>audioFormat</code> */
    protected float sampleRate;

    /** Source of input data.
     *  Could be extended to include live input from the sound card. */
    protected String audioFileName;

    /** Spacing of audio frames (determines the amount of overlap or skip
     *  between frames). This value is expressed in seconds and can be set by
     *  the command line option <b>-h hopTime</b>. (Default = 0.020s) */
    protected double hopTime;

    /** The approximate size of an FFT frame in seconds, as set by the command
     *  line option <b>-f FFTTime</b>. (Default = 0.04644s).  The value is
     *  adjusted so that <code>fftSize</code> is always power of 2. */
    protected double fftTime;

    /** Spacing of audio frames in samples (see <code>hopTime</code>) */
    protected int hopSize;

    /** The size of an FFT frame in samples (see <code>fftTime</code>) */
    protected int fftSize;

    /** The number of overlapping frames of audio data which have been read. */
    protected int frameCount;

    /** RMS amplitude of the current frame. */
    protected double frameRMS;

    /** Long term average frame energy (in frequency domain representation). */
    protected double ltAverage;

    /** Audio data is initially read in PCM format into this buffer. */
    protected byte[] inputBuffer;

    /** Audio data is scaled to the range [0,1] and averaged to one channel and
     *  stored in a circular buffer for reuse (if hopTime &lt; fftTime). */
    protected double[] circBuffer;

    /** The index of the next position to write in the circular buffer. */
    protected int cbIndex;

    /** The window function for the STFT, currently a Hamming window. */
    protected double[] window;

    /** The real part of the data for the in-place FFT computation.
     *  Since input data is real, this initially contains the input data. */
    protected double[] reBuffer;

    /** The imaginary part of the data for the in-place FFT computation.
     *  Since input data is real, this initially contains zeros. */
    protected double[] imBuffer;

    /** Phase of the previous frame, for calculating an onset function
     *  based on spectral phase deviation. */
    protected double[] prevPhase;

    /** Phase of the frame before the previous frame, for calculating an
     *  onset function based on spectral phase deviation. */
    protected double[] prevPrevPhase;

    /** Phase deviation onset detection function, indexed by frame. */
    protected double[] phaseDeviation;

    /** Spectral flux onset detection function, indexed by frame. */
    protected double[] spectralFlux;

    /** A mapping function for mapping FFT bins to final frequency bins.
     *  The mapping is linear (1-1) until the resolution reaches 2 points per
     *  semitone, then logarithmic with a semitone resolution.  e.g. for
     *  44.1kHz sampling rate and fftSize of 2048 (46ms), bin spacing is
     *  21.5Hz, which is mapped linearly for bins 0-34 (0 to 732Hz), and
     *  logarithmically for the remaining bins (midi notes 79 to 127, bins 35 to
     *  83), where all energy above note 127 is mapped into the final bin. */
    protected int[] freqMap;

    /** The number of entries in <code>freqMap</code>. Note that the length of
     *  the array is greater, because its size is not known at creation time. */
    protected int freqMapSize;

    /** The magnitude spectrum of the most recent frame.
     *  Used for calculating the spectral flux. */
    protected double[] prevFrame;

    /** The magnitude spectrum of the current frame. */
    protected double[] newFrame;

    /** The magnitude spectra of all frames, used for plotting the spectrogram. */
    protected double[][] frames;

    /** The RMS energy of all frames. */
    protected double[] energy;

    /** The estimated onset times from peak-picking the onset detection function(s). */
    protected double[] onsets;

    /** The estimated onset times and their saliences. */
    public EventList onsetList;

    /** GUI component which shows progress of audio processing. */
//	protected ProgressIndicator progressCallback;

    /** Total number of audio frames if known, or -1 for live or compressed input. */
    protected int totalFrames;

    /** Standard input for interactive prompts (for debugging). */
    BufferedReader stdIn;

//	/** Object for plotting output (for debugging / development) . */
//	Plot plot;

    /** Flag for enabling or disabling debugging output */
    public static boolean debug = false;

//	/** Flag for plotting onset detection function. */
//	public static boolean doOnsetPlot = false;

    /** Flag for suppressing all standard output messages except results. */
    protected static boolean silent = true;

    /** Flag for batch mode. */
    public static boolean batchMode = false;

    /** RMS frame energy below this value results in the frame being set to zero,
     *  so that normalisation does not have undesired side-effects. */
    public static double silenceThreshold = 0.0004;

    /** For dynamic range compression, this value is added to the log magnitude
     *  in each frequency bin and any remaining negative values are then set to zero.
     */
    public static double rangeThreshold = 10;

    /** Determines method of normalisation. Values can be:<ul>
     *  <li>0: no normalisation</li>
     *  <li>1: normalisation by current frame energy</li>
     *  <li>2: normalisation by exponential average of frame energy</li>
     *  </ul>
     */
    public static int normaliseMode = 2;

    /** Ratio between rate of sampling the signal energy (for the amplitude envelope) and the hop size */
    public static int energyOversampleFactor = 2;

    /** Audio buffer for live input. (Not used yet) */
    public static final int liveInputBufferSize = 32768; /* ~195ms buffer @CD */

    /** Maximum file length in seconds. Used for static allocation of arrays. */
    public static final int MAX_LENGTH = 3600;	// i.e. 1 hour


    /** Constructor: note that streams are not opened until the input file is set
     *  (see <code>setInputFile()</code>). */
    public AudioProcessor() {
        cbIndex = 0;
        frameRMS = 0;
        ltAverage = 0;
        frameCount = 0;
        hopSize = 0;
        fftSize = 0;
        hopTime = 0.010;	// DEFAULT, overridden with -h
        fftTime = 0.04644;	// DEFAULT, overridden with -f
//		progressCallback = null;
        stdIn = new BufferedReader(new InputStreamReader(System.in));
//		if (doOnsetPlot)
//			plot = new Plot();
    } // constructor

    /** For debugging, outputs information about the AudioProcessor to
     *  standard error.
     */
    public void print() {
        System.err.println(this);
    } // print()

    /** For interactive pause - wait for user to hit Enter */
    public String readLine() {
        try { return stdIn.readLine(); } catch (Exception e) { return null; }
    } // readLine()

    /** Gives some basic information about the audio being processed. */
    public String toString() {
        return "AudioProcessor\n" +
                String.format("\tFile: %s (%3.1f kHz, %1d channels)\n",
                        audioFileName, sampleRate/1000, channels) +
                String.format("\tHop / FFT sizes: %5.3f / %5.3f",
                        hopTime, hopTime * fftSize / hopSize);
    } // toString()

//	/** Adds a link to the GUI component which shows the progress of matching.
//	 *  @param c the AudioProcessor representing the other performance
//	 */
//	public void setProgressCallback(ProgressIndicator c) {
//		progressCallback = c;
//	} // setProgressCallback()

//	/** Sets up the streams and buffers for live audio input (CD quality).
//	 *  If any Exception is thrown within this method, it is caught, and any
//	 *  opened streams are closed, and <code>pcmInputStream</code> is set to
//	 *  <code>null</code>, indicating that the method did not complete
//	 *  successfully.
//	 */
//	public void setLiveInput() {
//		try {
//			channels = 2;
//			sampleRate = 44100;
//			AudioFormat desiredFormat = new AudioFormat(
//						AudioFormat.Encoding.PCM_SIGNED, sampleRate, 16,
//						channels, channels * 2, sampleRate, false);
//			TargetDataLine tdl = AudioSystem.getTargetDataLine(desiredFormat);
//			tdl.open(desiredFormat, liveInputBufferSize);
//			pcmInputStream = new AudioInputStream(tdl);
//			audioFormat = pcmInputStream.getFormat();
//			init();
//			tdl.start();
//		} catch (Exception e) {
//			e.printStackTrace();
//			closeStreams();	// make sure it exits in a consistent state
//		}
//	} // setLiveInput()

    /**
     * Sets the input stream directly.
     * Added by nw.
     *
     * @param ais AudioInputStream in 16bit, stereo, 44.1kHz, pcm signed
     */
    public void setInputStream(AudioInputStream ais) {

        rawInputStream = ais;
        pcmInputStream = ais;

        audioFormat = rawInputStream.getFormat();
        channels = audioFormat.getChannels();
        sampleRate = audioFormat.getSampleRate();

        audioFormat = new AudioFormat(
                AudioFormat.Encoding.PCM_SIGNED, sampleRate, 16,
                channels, channels * 2, sampleRate, false);

        //pcmInputStream = AudioSystem.getAudioInputStream(audioFormat, rawInputStream);

        init();
    }

    /** Allocates memory for arrays, based on parameter settings */
    protected void init() {
        hopSize = (int) Math.round(sampleRate * hopTime);
        fftSize = (int) Math.round(Math.pow(2,
                Math.round( Math.log(fftTime * sampleRate) / Math.log(2))));
        makeFreqMap(fftSize, sampleRate);
        int buffSize = hopSize * channels * 2;
        if ((inputBuffer == null) || (inputBuffer.length != buffSize))
            inputBuffer = new byte[buffSize];
        if ((circBuffer == null) || (circBuffer.length != fftSize)) {
            circBuffer = new double[fftSize];
            reBuffer = new double[fftSize];
            imBuffer = new double[fftSize];
            prevPhase = new double[fftSize];
            prevPrevPhase = new double[fftSize];
            prevFrame = new double[fftSize];
            window = FFT.makeWindow(FFT.HAMMING, fftSize, fftSize);
            for (int i=0; i < fftSize; i++)
                window[i] *= Math.sqrt(fftSize);
        }

        // Had to fix this. See repo for old version.
        totalFrames = (int)(pcmInputStream.getFrameLength() / hopSize);

        if(totalFrames == 0) {
            totalFrames = (int) (MAX_LENGTH / hopTime);
        }


        if ((newFrame == null) || (newFrame.length != freqMapSize)) {
            newFrame = new double[freqMapSize];
            frames = new double[totalFrames][freqMapSize];
        } else if (frames.length != totalFrames)
            frames = new double[totalFrames][freqMapSize];
        energy = new double[totalFrames*energyOversampleFactor];
        phaseDeviation = new double[totalFrames];
        spectralFlux = new double[totalFrames];
        frameCount = 0;
        cbIndex = 0;
        frameRMS = 0;
        ltAverage = 0;
//		progressCallback = null;
    } // init()

    /** Closes the input stream(s) associated with this object. */
    public void closeStreams() {
        if (pcmInputStream != null) {
            try {
                pcmInputStream.close();
                if (pcmInputStream != rawInputStream)
                    rawInputStream.close();
//				if (audioOut != null) {
//					audioOut.drain();
//					audioOut.close();
//				}
            } catch (Exception e) {}
            pcmInputStream = null;
//			audioOut = null;
        }
    } // closeStreams()

    /** Creates a map of FFT frequency bins to comparison bins.
     *  Where the spacing of FFT bins is less than 0.5 semitones, the mapping is
     *  one to one. Where the spacing is greater than 0.5 semitones, the FFT
     *  energy is mapped into semitone-wide bins. No scaling is performed; that
     *  is the energy is summed into the comparison bins. See also
     *  processFrame()
     */
    protected void makeFreqMap(int fftSize, float sampleRate) {
        freqMap = new int[fftSize/2+1];
        double binWidth = sampleRate / fftSize;
        int crossoverBin = (int)(2 / (Math.pow(2, 1/12.0) - 1));
        int crossoverMidi = (int)Math.round(Math.log(crossoverBin*binWidth/440)/
                Math.log(2) * 12 + 69);
        // freq = 440 * Math.pow(2, (midi-69)/12.0) / binWidth;
        int i = 0;
        while (i <= crossoverBin)
            freqMap[i++] = i;
        while (i <= fftSize/2) {
            double midi = Math.log(i*binWidth/440) / Math.log(2) * 12 + 69;
            if (midi > 127)
                midi = 127;
            freqMap[i++] = crossoverBin + (int)Math.round(midi) - crossoverMidi;
        }
        freqMapSize = freqMap[i-1] + 1;
    } // makeFreqMap()

    /** Calculates the weighted phase deviation onset detection function.
     *  Not used.
     *  TODO: Test the change to WPD fn */
    protected void weightedPhaseDeviation() {
        if (frameCount < 2)
            phaseDeviation[frameCount] = 0;
        else {
            for (int i = 0; i < fftSize; i++) {
                double pd = imBuffer[i] - 2 * prevPhase[i] + prevPrevPhase[i];
                double pd1 = Math.abs(Math.IEEEremainder(pd, 2 * Math.PI));
                phaseDeviation[frameCount] += pd1 * reBuffer[i];
                // System.err.printf("%7.3f   %7.3f\n", pd/Math.PI, pd1/Math.PI);
            }
        }
        phaseDeviation[frameCount] /= fftSize * Math.PI;
        double[] tmp = prevPrevPhase;
        prevPrevPhase = prevPhase;
        prevPhase = imBuffer;
        imBuffer = tmp;
    } // weightedPhaseDeviation()

    /** Reads a frame of input data, averages the channels to mono, scales
     *  to a maximum possible absolute value of 1, and stores the audio data
     *  in a circular input buffer.
     *  @return true if a frame (or part of a frame, if it is the final frame)
     *  is read. If a complete frame cannot be read, the InputStream is set
     *  to null.
     */
    public boolean getFrame() {
        if (pcmInputStream == null)
            return false;
        try {
            int bytesRead = (int) pcmInputStream.read(inputBuffer);
            if (bytesRead < inputBuffer.length) {
                closeStreams();
                return false;
            }
        } catch (IOException e) {
            e.printStackTrace();
            closeStreams();
            return false;
        }
        frameRMS = 0;
        double sample;
        switch(channels) {
            case 1:
                for (int i = 0; i < inputBuffer.length; i += 2) {
                    sample = ((inputBuffer[i+1]<<8) |
                            (inputBuffer[i]&0xff)) / 32768.0;
                    frameRMS += sample * sample;
                    circBuffer[cbIndex++] = sample;
                    if (cbIndex == fftSize)
                        cbIndex = 0;
                }
                break;
            case 2: // saves ~0.1% of RT (total input overhead ~0.4%) :)
                for (int i = 0; i < inputBuffer.length; i += 4) {
                    sample = (((inputBuffer[i+1]<<8) | (inputBuffer[i]&0xff)) +
                            ((inputBuffer[i+3]<<8) | (inputBuffer[i+2]&0xff)))
                            / 65536.0;
                    frameRMS += sample * sample;
                    circBuffer[cbIndex++] = sample;
                    if (cbIndex == fftSize)
                        cbIndex = 0;
                }
                break;
            default:
                for (int i = 0; i < inputBuffer.length; ) {
                    sample = 0;
                    for (int j = 0; j < channels; j++, i+=2)
                        sample += (inputBuffer[i+1]<<8) | (inputBuffer[i]&0xff);
                    sample /= 32768.0 * channels;
                    frameRMS += sample * sample;
                    circBuffer[cbIndex++] = sample;
                    if (cbIndex == fftSize)
                        cbIndex = 0;
                }
        }
        frameRMS = Math.sqrt(frameRMS / inputBuffer.length * 2 * channels);
        return true;
    } // getFrame()

    /** Processes a frame of audio data by first computing the STFT with a
     *  Hamming window, then mapping the frequency bins into a part-linear
     *  part-logarithmic array, then computing the spectral flux
     *  then (optionally) normalising and calculating onsets.
     */
    protected void processFrame() {
        if (getFrame()) {
            for (int i = 0; i < fftSize; i++) {
                reBuffer[i] = window[i] * circBuffer[cbIndex];
                if (++cbIndex == fftSize)
                    cbIndex = 0;
            }
            Arrays.fill(imBuffer, 0);
            FFT.magnitudePhaseFFT(reBuffer, imBuffer);
            Arrays.fill(newFrame, 0);
            double flux = 0;
            for (int i = 0; i <= fftSize/2; i++) {
                if (reBuffer[i] > prevFrame[i])
                    flux += reBuffer[i] - prevFrame[i];
                newFrame[freqMap[i]] += reBuffer[i];
            }
            spectralFlux[frameCount] = flux;
            for (int i = 0; i < freqMapSize; i++)
                frames[frameCount][i] = newFrame[i];
            int index = cbIndex - (fftSize - hopSize);
            if (index < 0)
                index += fftSize;
            int sz = (fftSize - hopSize) / energyOversampleFactor;
            for (int j = 0; j < energyOversampleFactor; j++) {
                double newEnergy = 0;
                for (int i = 0; i < sz; i++) {
                    newEnergy += circBuffer[index] * circBuffer[index];
                    if (++index == fftSize)
                        index = 0;
                }
                energy[frameCount * energyOversampleFactor + j] =
                        newEnergy / sz <= 1e-6? 0: Math.log(newEnergy / sz) + 13.816;
            }
            double decay = frameCount >= 200? 0.99:
                    (frameCount < 100? 0: (frameCount - 100) / 100.0);
            if (ltAverage == 0)
                ltAverage = frameRMS;
            else
                ltAverage = ltAverage * decay + frameRMS * (1.0 - decay);
            if (frameRMS <= silenceThreshold)
                for (int i = 0; i < freqMapSize; i++)
                    frames[frameCount][i] = 0;
            else {
                if (normaliseMode == 1)
                    for (int i = 0; i < freqMapSize; i++)
                        frames[frameCount][i] /= frameRMS;
                else if (normaliseMode == 2)
                    for (int i = 0; i < freqMapSize; i++)
                        frames[frameCount][i] /= ltAverage;
                for (int i = 0; i < freqMapSize; i++) {
                    frames[frameCount][i] = Math.log(frames[frameCount][i]) + rangeThreshold;
                    if (frames[frameCount][i] < 0)
                        frames[frameCount][i] = 0;
                }
            }
            double[] tmp = prevFrame;
            prevFrame = reBuffer;
            reBuffer = tmp;
            frameCount++;
        }
    } // processFrame()

    /** Processes a complete file of audio data. */
    public void processFile() {
        while (pcmInputStream != null) {
            // Profile.start(0);
            processFrame();
            // Profile.log(0);
            if (Thread.currentThread().isInterrupted()) {
                System.err.println("info: INTERRUPTED in processFile()");
                return;
            }
        }

        double hop = hopTime;
        Peaks.normalise(spectralFlux);
        LinkedList<Integer> peaks = Peaks.findPeaks(spectralFlux, (int)Math.round(0.06 / hop), 0.35, 0.84, true);
        onsets = new double[peaks.size()];
        double[] y2 = new double[onsets.length];
        Iterator<Integer> it = peaks.iterator();
        onsetList = new EventList();
        double minSalience = Peaks.min(spectralFlux);
        for (int i = 0; i < onsets.length; i++) {
            int index = it.next();
            onsets[i] = index * hop;
            y2[i] = spectralFlux[index];
            Event e = newBeat(onsets[i], 0);
            e.salience = spectralFlux[index] - minSalience;
            onsetList.add(e);
        }
    }

    /** Creates a new Event object representing a beat.
     *  @param time The time of the beat in seconds
     *  @param beatNum The index of the beat
     *  @return The Event object representing the beat
     */
    public static Event newBeat(double time, int beatNum) {
        return new Event(time,time, time, 56, 64, beatNum, 0, 1);
    } // newBeat()


    /** Perform beat tracking where the GUI is not active;
     *  there is no selected region.
     *  @param events The onsets or peaks in a feature list
     *  @param beats The initial beats which are given, if any
     *  @return The list of beats, or an empty list if beat tracking fails
     */
    public static EventList beatTrack(EventList events, EventList beats) {
        AgentList agents = null;
        int count = 0;
        double beatTime = -1;
        if (beats != null) {
            count = beats.size() - 1;
            beatTime = beats.l.getLast().keyDown;
        }
        if (count > 0) { // tempo given by mean of initial beats
            double ioi = (beatTime - beats.l.getFirst().keyDown) / count;
            agents = new AgentList(new Agent(ioi), null);
        } else									// tempo not given; use tempo induction
            agents = Induction.beatInduction(events);
        if (beats != null)
            for (AgentList ptr = agents; ptr.ag != null; ptr = ptr.next) {
                ptr.ag.beatTime = beatTime;
                ptr.ag.beatCount = count;
                ptr.ag.events = new EventList(beats);
            }

        agents.beatTrack(events, -1);

        Agent best = agents.bestAgent();

        if (best != null) {
            best.fillBeats(beatTime);
            return best.events;
        }

        return new EventList();
    }


    /** Finds the median tempo (as inter-beat interval) from an array of beat times
     *  @param d An array of beat times
     *  @return The median inter-beat interval
     */
    public static double getMedianIBI(double[] d) {
        if ((d == null) || (d.length < 2))
            return -1.0;
        double[] ibi = new double[d.length-1];
        for (int i = 1; i < d.length; i++)
            ibi[i-1] = d[i] - d[i-1];
        Arrays.sort(ibi);
        if (ibi.length % 2 == 0)
            return (ibi[ibi.length / 2] + ibi[ibi.length / 2 - 1]) / 2;
        else
            return ibi[ibi.length / 2];
    }

}