package com.company;

import com.company.classes.BeatRoot;

import javax.sound.sampled.*;
import java.io.File;
import java.io.IOException;

public class Main {

    public static void main(String[] args) throws IOException, UnsupportedAudioFileException {
	    File file = new File("src/com/company/resources/test.wav");
        AudioInputStream audioInputStream = AudioSystem.getAudioInputStream(file);
        BeatRoot beatRoot = new BeatRoot();
        Integer bpm = beatRoot.calculateBPM(audioInputStream);
        System.out.println("BPM: " + bpm);
    }
}
