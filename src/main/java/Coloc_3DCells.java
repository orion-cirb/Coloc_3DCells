/*
 * Find dapi, cells in channels (PML) in outlined Cells 
 * Measure integrated intensity, nb of dots per cells detect if cell infected 
 * Author Philippe Mailly
 */


import ij.*;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
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
import java.util.List;
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


public class Coloc_3DCells implements PlugIn {

    private final boolean canceled = false;
    private String imageDir = "";
    public static String outDirResults = "";
    private BufferedWriter nucleusResults;
    
    private Cells_Processing proc = new Cells_Processing();

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
            ArrayList<String> channels = proc.findChannels(imageFiles.get(0));
            channels.add("None");
            List<String> chs = new ArrayList();
            List<String> channelsName = new ArrayList();
            channelsName.add("Nucleus");
            channelsName.add("GFP");
            channelsName.add("CC1");
            channelsName.add("NG2");

            // Channels dialog
            chs = proc.dialog(channels, channelsName);
            
            if (chs == null || proc.canceled) {
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
                if (roiFile != null) {
                    Roi roi = new Opener().openRoi(roiFile);
                    options.setCrop(true);
                    Region reg = new Region(roi.getBounds().x, roi.getBounds().y, roi.getBounds().width, roi.getBounds().height);
                    options.setCropRegion(0, reg);
                    options.doCrop();
                }
                // find dapi and cells channel images
                
                // Open DAPI Channel and detect nucleus
                int indexCh = channels.indexOf(chs.get(0));
                System.out.println("Opening nucleus channel " + chs.get(0) +" ...");
                ImagePlus imgNuc = BF.openImagePlus(options)[indexCh];
                Objects3DIntPopulation nucPop = proc.stardistCellsPop(imgNuc);
                System.out.println("Total nucleus found = "+nucPop.getNbObjects());
                
                // Initialize infos in nucleus class
                for (Object3DInt nucObj : nucPop.getObjects3DInt()) {
                    double vol = new MeasureVolume(nucObj).getVolumeUnit();
                    nuclei.add(new Nucleus(nucObj.getLabel(), vol, false, false, 0, 0));
                }
                
                // Open gfp image
                ImagePlus imgGFP;
                Objects3DIntPopulation gfpFilterPop = new Objects3DIntPopulation();
                indexCh = channels.indexOf(chs.get(1));
                if (!chs.get(1).equals("None")) {
                    System.out.println("Opening GFP channel " + chs.get(1) +" ...");
                    imgGFP = BF.openImagePlus(options)[indexCh];
                    Objects3DIntPopulation gfpPop = proc.cellPoseCellsPop(imgGFP, proc.gfpCellModel, proc.gfpCellPoseDiameter);
                    System.out.println("Total GFP = "+gfpPop.getNbObjects());
                    gfpFilterPop = proc.intensityFilter(gfpPop, imgGFP, proc.gfpIntTh);
                    System.out.println("Total GFP after intensity filter = "+gfpFilterPop.getNbObjects());
                    proc.closeImages(imgGFP);
                }
                
                // Open cc1 image
                ImagePlus imgCC1;
                Objects3DIntPopulation cc1FilterPop = new Objects3DIntPopulation();
                indexCh = channels.indexOf(chs.get(2));
                if (!chs.get(2).equals("None")) {
                    System.out.println("Opening CC1 channel " + channels.get(2) +" ...");
                    imgCC1 = BF.openImagePlus(options)[indexCh];
                    Objects3DIntPopulation cc1Pop = proc.cellPoseCellsPop(imgCC1, proc.cc1CellModel, proc.cc1CellPoseDiameter);
                    System.out.println("Total CC1 = "+cc1Pop.getNbObjects());
                    cc1FilterPop = proc.intensityFilter(cc1Pop, imgCC1, proc.cc1IntTh);
                    System.out.println("Total CC1 after intensity filter = "+cc1FilterPop.getNbObjects());
                    proc.closeImages(imgCC1);
                }
                
                // Open NG2 image
                ImagePlus imgNG2 = null;
                indexCh = channels.indexOf(chs.get(3));
                double bgNG2 = 0;
                if (!chs.get(3).equals("None")) {
                    System.out.println("Opening NG2 channel " + channels.get(3) +" ...");
                    imgNG2 = BF.openImagePlus(options)[indexCh];
                    // find background of NG2 as mean + std of min projection
                    double[] bg = proc.find_background(imgNG2);
                    bgNG2 = bg[0]+ bg[1];
                }
                    
                
                // Find colocalized cells

                // dapi/gfp
                System.out.println("Finding dapi/GFP colocalization ...");
                Objects3DIntPopulation dapiGFPPop = proc.findColoc(nucPop, gfpFilterPop, channelsName.get(1), nuclei);
                System.out.println(dapiGFPPop.getNbObjects() + " dapi with GFP found");
                
                // dapi/cc1
                Objects3DIntPopulation dapiCC1Pop = proc.findColoc(nucPop, cc1FilterPop, channelsName.get(2), nuclei);
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
                if (gfpFilterPop.getNbObjects() > 0)
                    gfpFilterPop.drawInImage(imhGFP);
                if (cc1FilterPop.getNbObjects() > 0)
                    cc1FilterPop.drawInImage(imhCC1);
                ImagePlus[] imgColors =  {imhCC1.getImagePlus(), imhGFP.getImagePlus(), imhNuc.getImagePlus(), imgNG2, imhNucDil.getImagePlus()};
                ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                imgObjects.setCalibration(proc.cal);
                FileSaver imgObjectsFile = new FileSaver(imgObjects);
                imgObjectsFile.saveAsTiff(outDirResults+rootName+"_Objects.tif");
                proc.closeImages(imgNuc);
                proc.closeImages(imgObjects);
                proc.closeImages(imgNG2);
                
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
