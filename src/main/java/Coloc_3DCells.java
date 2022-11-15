/*
 * Find dapi, cells in channels (PML) in outlined Cells 
 * Measure integrated intensity, nb of dots per cells detect if cell infected 
 * Author Philippe Mailly
 */


import Tools.Nucleus;
import ij.*;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.Region;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;
import Tools.Process;
import java.util.Arrays;



public class Coloc_3DCells implements PlugIn {
    
    Tools.Process proc = new Process();

    private final boolean canceled = false;
    private String imageDir = "";
    public static String outDirResults = "";
    private BufferedWriter nucleusResults;
    
    // StarDist
    private String stardistModel = "StandardFluo.zip";
    private double nucProbThr = 0.40;
    private double nucOver = 0.25;
    
    // CellPose
    int cellPoseDiameter = 30;
    double cellPoseMaskTh = 0;
    double cellPoseFlowTh = 0.4;
    double cellPoseStitchTh = 0.25;
    private boolean useGpu = true;
    public String cellGFPModel = "cyto2";
    public String cellCC1Model = "cyto2";
    
    // min max size for nucleus µm3
    public double minNucVol = 10;
    public double maxNucVol = 1000;
    // min max size for cell µm3
    public int minCellVol = 10;
    public int maxCellVol = 5000;
    
    // gene Intensity threshold
    public double gfpIntTh = 100;
    public double cc1IntTh = gfpIntTh;
            
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }
            // Find images files
            String fileExt = proc.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = proc.findImages(imageDir, fileExt);
            if (imageFiles == null) {
                IJ.showMessage("Error", "No images found with "+fileExt+" extension");
                return;
            }
            // create output folder
            outDirResults = imageDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            
             // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            // Find calibration
            reader.setId(imageFiles.get(0));
            proc.cal = proc.findImageCalib(meta);
            
            // Find channel names
            String[] imgChs = proc.findChannels(imageFiles.get(0), meta, reader);
            ArrayList<String> channels = new ArrayList<>(Arrays.asList(imgChs));
            channels.add("None");
            imgChs = channels.toArray(imgChs);
            String[] channelsName = {"Nucleus", "GFP", "CC1", "NG2"};

            // Channels dialog
            String[] chs = proc.dialog(stardistModel, imgChs, channelsName, minNucVol, maxNucVol, minCellVol, maxCellVol, gfpIntTh, cc1IntTh );
                
            if (chs == null) {
                    IJ.showStatus("Plugin cancelled");
                    return;
            }
                
            // Write headers results for results file
            FileWriter fileResultsNumber = new FileWriter(outDirResults + "CellsNumberResults.xls", false);
            nucleusResults = new BufferedWriter(fileResultsNumber);
            try {
                nucleusResults.write("ImageName\t#Nucleus\tVol\tIs GFP\tIs CC1\tNG2 mean intensity\tNG2 mean corrected intensity\n");
                nucleusResults.flush();
            }
            catch (IOException ex) {
                Logger.getLogger(Coloc_3DCells.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            
            // Read images
            for (String f : imageFiles) {
                reader.setId(f);
                String rootName = FilenameUtils.getBaseName(f);
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                
                ArrayList<Nucleus> nuclei = new ArrayList<>();
                // find rois
                String roiFile = imageDir+rootName+".roi";
                if (new File(roiFile).exists()) {
                    Roi roi = new Opener().openRoi(roiFile);
                    options.setCrop(true);
                    Region reg = new Region(roi.getBounds().x, roi.getBounds().y, roi.getBounds().width, roi.getBounds().height);
                    options.setCropRegion(0, reg);
                    options.doCrop();
                }
                // find dapi and cells channel images
                
                // Open DAPI Channel and detect nucleus
                int indexCh = ArrayUtils.indexOf(imgChs, chs[0]);
                System.out.println("Opening nucleus channel " + chs[0] +" ...");
                ImagePlus imgNuc = BF.openImagePlus(options)[indexCh];
                Objects3DIntPopulation nucPop = proc.stardistObjectsPop(imgNuc, 0.5f, true, 5, stardistModel, nucProbThr, nucOver);
                System.out.println("Total nucleus found = "+nucPop.getNbObjects());
                proc.popFilterSize(nucPop, minNucVol, maxNucVol);
                System.out.println("Total nucleus found aflter size filter "+nucPop.getNbObjects());
                // Initialize infos in nucleus class
                for (Object3DInt nucObj : nucPop.getObjects3DInt()) {
                    double vol = new MeasureVolume(nucObj).getVolumeUnit();
                    nuclei.add(new Nucleus(nucObj.getLabel(), vol, false, false, 0, 0));
                }
                
                // Open gfp image
                ImagePlus imgGFP;
                Objects3DIntPopulation gfpPop = new Objects3DIntPopulation();
                indexCh = ArrayUtils.indexOf(imgChs, chs[1]);
                if (!chs[1].equals("None")) {
                    System.out.println("Opening GFP channel " + chs[1] +" ...");
                    imgGFP = BF.openImagePlus(options)[indexCh];
                    gfpPop = proc.cellposeDetection(imgGFP, cellGFPModel, cellPoseDiameter, cellPoseMaskTh, cellPoseFlowTh, cellPoseStitchTh, 0.5f, true, useGpu);
                    System.out.println("Total GFP "+gfpPop.getNbObjects());
                    proc.popFilterSize(gfpPop, minCellVol, maxCellVol);
                    System.out.println("Total GFP after size filter "+gfpPop.getNbObjects());
                    proc.intensityFilter(gfpPop, imgGFP, gfpIntTh);
                    System.out.println("Total GFP after intensity filter "+gfpPop.getNbObjects());
                    proc.flush_close(imgGFP);
                }
                
                // Open cc1 image
                ImagePlus imgCC1;
                Objects3DIntPopulation cc1Pop = new Objects3DIntPopulation();
                indexCh = ArrayUtils.indexOf(imgChs, chs[2]);
                if (!chs[2].equals("None")) {
                    System.out.println("Opening CC1 channel " + chs[2] +" ...");
                    imgCC1 = BF.openImagePlus(options)[indexCh];
                    cc1Pop = proc.cellposeDetection(imgCC1, cellCC1Model, cellPoseDiameter, cellPoseMaskTh, cellPoseFlowTh, cellPoseStitchTh, 0.5f, true, useGpu);
                    System.out.println("Total CC1 "+cc1Pop.getNbObjects());
                    proc.popFilterSize(cc1Pop, minCellVol, maxCellVol);
                    System.out.println("Total CC1 after size filter "+cc1Pop.getNbObjects());
                    proc.intensityFilter(cc1Pop, imgCC1, cc1IntTh);
                    System.out.println("Total CC1 after intensity filter "+cc1Pop.getNbObjects());
                    proc.flush_close(imgCC1);
                }
                
                // Open NG2 image
                ImagePlus imgNG2 = null;
                indexCh = ArrayUtils.indexOf(imgChs, chs[3]);
                double bgNG2 = 0;
                if (!chs[3].equals("None")) {
                    System.out.println("Opening NG2 channel " + chs[3] +" ...");
                    imgNG2 = BF.openImagePlus(options)[indexCh];
                    // find background of NG2
                    bgNG2 = proc.findMedianBackground(imgNG2);
                }
                    
                
                // Find colocalized cells

                // dapi/gfp
                System.out.println("Finding dapi/GFP colocalization ...");
                Objects3DIntPopulation dapiGFPPop = proc.findColoc(nucPop, gfpPop, channelsName[1], nuclei);
                System.out.println(dapiGFPPop.getNbObjects() + " dapi with GFP found");
                
                // dapi/cc1
                Objects3DIntPopulation dapiCC1Pop = proc.findColoc(nucPop, cc1Pop, channelsName[2], nuclei);
                System.out.println(dapiCC1Pop.getNbObjects() + " dapi with CC1 found");
                               
                // Measure intensity of dapi, coloc dapi with gfp and cc1 after dapi dilatation
                Objects3DIntPopulation nucDilPop = proc.measureIntensity(nucPop, imgNG2, nuclei, bgNG2);
                
                // Save objects image
                ImageHandler imhNuc = ImageHandler.wrap(imgNuc).createSameDimensions();
                ImageHandler imhNucDil = imhNuc.createSameDimensions();
                ImageHandler imhGFP = imhNuc.createSameDimensions();
                ImageHandler imhCC1 = imhNuc.createSameDimensions();
                // draw objects pop
                nucPop.drawInImage(imhNuc);
                nucDilPop.drawInImage(imhNucDil);
                if (gfpPop.getNbObjects() > 0)
                    gfpPop.drawInImage(imhGFP);
                if (cc1Pop.getNbObjects() > 0)
                    cc1Pop.drawInImage(imhCC1);
                ImagePlus[] imgColors =  {imhCC1.getImagePlus(), imhGFP.getImagePlus(), imhNuc.getImagePlus(), imgNG2, imhNucDil.getImagePlus()};
                ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                imgObjects.setCalibration(proc.cal);
                FileSaver imgObjectsFile = new FileSaver(imgObjects);
                imgObjectsFile.saveAsTiff(outDirResults+rootName+"_Objects.tif");
                proc.flush_close(imgNuc);
                proc.flush_close(imgObjects);
                proc.flush_close(imgNG2);
                
                // write results
                IJ.showStatus("Writing results ...");
                
                for (Nucleus nuc : nuclei) {
                   nucleusResults.write(rootName+"\t"+nuc.getLabel()+"\t"+nuc.getVol()+"\t"+nuc.getisGFP()+"\t"+nuc.getisCC1()+"\t"+
                           nuc.getNG2MeanInt()+"\t"+nuc.getNG2MeanIntCor()+"\n");
                    nucleusResults.flush();
                }
                   
            }
            
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Coloc_3DCells.class.getName()).log(Level.SEVERE, null, ex);
        }
        try {
            nucleusResults.close();
        } catch (IOException ex) {
            Logger.getLogger(Coloc_3DCells.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done");
    }
}
