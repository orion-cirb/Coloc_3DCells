package Tools;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


/**
 *
 * @author phm
 */
public class Nucleus {
    
    // nucleus index
    private float label;
    // nucleus volume
    private double vol;
    // is gfp surface
    private boolean gfp;
    // is cc1
    private boolean cc1; 
    // Mean intensity after dilataion in NG2 channel
    private double NG2MeanInt;
    // Corrected mean intensity (minus BG mean + bg Std)
    // after dilataion in NG2 channeNumber of gene2 dots in cell
    private double NG2MeanIntCor;  

   

   
	
	public Nucleus(float label, double vol, boolean gfp, boolean cc1, double NG2MeanInt, double NG2MeanIntCor) {
            this.label = label;
            this.vol = vol;
            this.gfp = gfp;
            this.cc1 = cc1;
            this.NG2MeanInt = NG2MeanInt;
            this.NG2MeanIntCor = NG2MeanIntCor;
	}
        
        public void setLabel(float label) {
            this.label = label;
	}
        
        public void setVol(double vol) {
            this.vol = vol;
	}
        
        public void setisGFP(boolean gfp) {
            this.gfp = gfp;
	}
        
        public void setisCC1(boolean cc1) {
            this.cc1 = cc1;
	}
        
        public void setNG2MeanInt(double NG2MeanInt) {
            this.NG2MeanInt = NG2MeanInt;
	}
        
        public void setNG2MeanIntCor(double NG2MeanIntCor) {
            this.NG2MeanIntCor = NG2MeanIntCor;
        }
        
        public float getLabel() {
            return label;
        }
        
        public double getVol() {
            return vol;
        }
        
        public boolean getisGFP() {
            return gfp;
        }
                
        public boolean getisCC1() {
            return cc1;
	}
        
        public double getNG2MeanInt() {
            return NG2MeanInt;
        }
        
        public double getNG2MeanIntCor() {
            return NG2MeanIntCor;
        }
                
}
