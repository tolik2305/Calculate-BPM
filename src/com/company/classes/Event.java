package com.company.classes;

public class Event implements Comparable<Object>, Cloneable, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    public double keyDown, keyUp, pedalUp, scoreBeat, scoreDuration, salience;
    public int midiPitch, midiVelocity, flags, midiCommand, midiChannel,
            midiTrack;
    //public String label;

    public Event(double onset, double offset, double eOffset, int pitch,
                 int velocity, double beat, double duration, int eventFlags,
                 int command, int channel, int track) {
        this(onset, offset, eOffset, pitch, velocity, beat,duration,eventFlags);
        midiCommand = command;
        midiChannel = channel;
        midiTrack = track;
    } // constructor

    public Event(double onset, double offset, double eOffset, int pitch,
                 int velocity, double beat, double duration, int eventFlags) {
        keyDown = onset;
        keyUp = offset;
        pedalUp = eOffset;
        midiPitch = pitch;
        midiVelocity = velocity;
        scoreBeat = beat;
        scoreDuration = duration;
        flags = eventFlags;
        midiCommand = javax.sound.midi.ShortMessage.NOTE_ON;
        midiChannel = 1;
        midiTrack = 0;
        salience = 0;
    } // constructor

    public Event clone() {
        return new Event(keyDown, keyUp, pedalUp, midiPitch, midiVelocity,
                scoreBeat, scoreDuration, flags, midiCommand, midiChannel,
                midiTrack);
    } // clone()

    // Interface Comparable
    public int compareTo(Object o) {
        Event e = (Event) o;
        return (int)Math.signum(keyDown - e.keyDown);
    } // compareTo()

    public String toString() {
        return "n=" + midiPitch + " v=" + midiVelocity + " t=" + keyDown +
                " to " + keyUp + " (" + pedalUp + ")";
    } // toString()

    public void print(Flags f) {
        System.out.printf("Event:\n");
        System.out.printf("\tkeyDown / Up / pedalUp: %5.3f / %5.3f /  %5.3f\n",
                keyDown, keyUp, pedalUp);
        //System.out.printf("\tkeyUp: %5.3f\n", keyUp);
        //System.out.printf("\tpedalUp: %5.3f\n", pedalUp);
        System.out.printf("\tmidiPitch: %d\n", midiPitch);
        System.out.printf("\tmidiVelocity: %d\n", midiVelocity);
        System.out.printf("\tmidiCommand: %02x\t", midiCommand | midiChannel);
        //System.out.printf("\tmidiChannel: %d\n", midiChannel);
        System.out.printf("\tmidiTrack: %d\n", midiTrack);
        System.out.printf("\tsalience: %5.3f\t", salience);
        System.out.printf("\tscoreBeat: %5.3f\t", scoreBeat);
        System.out.printf("\tscoreDuration: %5.3f\n", scoreDuration);
        System.out.printf("\tflags: %X", flags);
        if (f != null) {
            int ff = flags;
            for (int i=0; ff != 0; i++) {
                if (ff % 2 == 1)
                    System.out.print(" " + f.getLabel(i));
                ff >>>= 1;
            }
        }
        System.out.print("\n\n");
    } // print()

} // class Event